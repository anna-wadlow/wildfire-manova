# ============================================================
# Build the wildfire dataset for MATH 257 MANOVA project
# California fires 2000-2024 (FPA FOD) + gridMET weather + 3DEP terrain
# Output: fires_ca_weather_terrain_joint_0024.csv  (read by the R script)
# ============================================================
import math
import os
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import requests
import xarray as xr
from scipy.ndimage import uniform_filter

# ---- Settings ----------------------------------------------
FOD_DB = "FPA_FOD_20260615.sqlite"      # https://www.fs.usda.gov/rds/archive/catalog/RDS-2013-0009.7
OUT_CSV = "fires_ca_weather_terrain_joint_0024.csv"

GRIDMET_URL = "http://www.northwestknowledge.net/metdata/data/"
GRIDMET_DIR = Path("data_raw/gridmet");  GRIDMET_DIR.mkdir(parents=True, exist_ok=True)
GRIDMET_VARS = {"tmmx": "max_temp_c", "rmin": "min_rh", "vs": "wind_speed", "erc": "erc"}
LAT_MIN, LAT_MAX, LON_MIN, LON_MAX = 32.0, 42.5, -125.0, -114.0   # California box

DEM_URL = "https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/1/TIFF/current/{t}/USGS_1_{t}.tif"
DEM_DIR = Path("dem_tiles");  DEM_DIR.mkdir(exist_ok=True)
WIN = 51 # slope smoothing window: 51 cells x ~30 m = ~1.5 km


# ---- 1. California fires, 2000-2024 ------------------------
query = """
SELECT FOD_ID, FIRE_SIZE, FIRE_SIZE_CLASS, LATITUDE, LONGITUDE, DISCOVERY_DATE, FIRE_YEAR
FROM Fires
WHERE STATE = 'CA' AND FIRE_YEAR >= 2000
"""
con = sqlite3.connect(FOD_DB)
fires = pd.read_sql(query, con)
con.close()

fires = fires.drop_duplicates(subset="FOD_ID")
fires["discovery_date"] = pd.to_datetime(fires["DISCOVERY_DATE"] - 2440587.5, unit="D").dt.normalize()  # julian -> date
fires["weather_date"]   = fires["discovery_date"] - pd.Timedelta(days=1)   # weather from the day before discovery

# size groups: small A-C (<100 acres), medium D-E (100-999), large F-G (1000+)
fires["size_group3"] = pd.cut(fires["FIRE_SIZE"], bins=[0, 100, 1000, np.inf], right=False,
                              labels=["small", "medium", "large"])
print(pd.crosstab(fires["FIRE_SIZE_CLASS"], fires["size_group3"]))   # check: A-C small, D-E medium, F-G large


# ---- 2. Weather: gridMET (4 km grid, nearest cell, day before) ----
def download_gridmet(var, year):
    path = GRIDMET_DIR / f"{var}_{year}.nc"
    if not path.exists():
        with requests.get(f"{GRIDMET_URL}{var}_{year}.nc", stream=True, timeout=120) as r:
            r.raise_for_status()
            with open(path, "wb") as f:
                for chunk in r.iter_content(chunk_size=1 << 20):
                    f.write(chunk)
    return path

def gridmet_values(var, year, rows):
    """gridMET value on weather_date at the cell nearest each fire in rows."""
    ds = xr.open_dataset(download_gridmet(var, year))
    da = ds[[v for v in ds.data_vars if ds[v].ndim == 3][0]]
    lat_slice = slice(LAT_MAX, LAT_MIN) if float(da.lat[0]) > float(da.lat[-1]) else slice(LAT_MIN, LAT_MAX)
    da = da.sel(lat=lat_slice, lon=slice(LON_MIN, LON_MAX)).load()
    tdim = [d for d in da.dims if d not in ("lat", "lon")][0]
    pts = da.sel(lat=xr.DataArray(rows["LATITUDE"].values, dims="pt"),
                 lon=xr.DataArray(rows["LONGITUDE"].values, dims="pt"),
                 **{tdim: xr.DataArray(rows["weather_date"].values, dims="pt")},
                 method="nearest")
    ds.close()
    return pd.Series(pts.values, index=rows.index)

fires["weather_year"] = fires["weather_date"].dt.year
for var, col in GRIDMET_VARS.items():
    fires[col] = np.nan
    for year, rows in fires.groupby("weather_year"):
        fires.loc[rows.index, col] = gridmet_values(var, year, rows)
        print(f"done: {var} {year}")

fires["max_temp_c"] -= 273.15  # kelvin -> celsius
fires = fires.dropna(subset=list(GRIDMET_VARS.values()))


# ---- 3. Terrain: USGS 3DEP 1 arc-second DEM (~30 m) ---------
def tile_name(lat, lon):
    """USGS tiles are named by their northwest corner, e.g. n38w121 covers lat 37-38, lon -121 to -120."""
    return f"n{math.ceil(lat):02d}w{math.ceil(abs(lon)):03d}"

def download_tile(t, tries=3):
    dest, part = DEM_DIR / f"USGS_1_{t}.tif", DEM_DIR / f"USGS_1_{t}.part"
    for attempt in range(tries):
        try:
            with requests.get(DEM_URL.format(t=t), stream=True, timeout=120) as r:
                r.raise_for_status()
                with open(part, "wb") as f:
                    for chunk in r.iter_content(chunk_size=1 << 20):
                        f.write(chunk)
            part.replace(dest)          # only rename once the download finished
            return True
        except Exception as e:
            print(f"{t}: attempt {attempt + 1} failed -> {e}")
    return False

def slope_degrees(z, lat_center, win=WIN):
    """Slope for every cell of a tile, after smoothing elevation over a win x win window."""
    valid = np.isfinite(z)
    smooth = uniform_filter(np.where(valid, z, 0).astype("float32"), size=win, mode="nearest")
    share_valid = uniform_filter(valid.astype("float32"), size=win, mode="nearest")
    z = np.where(share_valid > 1 - 0.5 / win**2, smooth, np.nan).astype("float32")  # NaN if any cell in window is missing
    dy = 30.87 # meters per arc-second of latitude
    dx = 30.87 * np.cos(np.radians(lat_center)) # longitude shrinks with latitude
    dz_dn = -np.gradient(z, dy, axis=0) # rows run north to south
    dz_de = np.gradient(z, dx, axis=1)
    return np.degrees(np.arctan(np.hypot(dz_de, dz_dn)))

fires["tile"] = [tile_name(la, lo) for la, lo in zip(fires["LATITUDE"], fires["LONGITUDE"])]

for t in sorted(fires["tile"].unique()):
    if not (DEM_DIR / f"USGS_1_{t}.tif").exists():
        download_tile(t) # a few tiles are open ocean and don't exist; those fires get NaN

fires[["elev_m", "slope_deg"]] = np.nan
fires["near_tile_edge"] = False
half = WIN // 2

for t, grp in fires.groupby("tile"):
    path = DEM_DIR / f"USGS_1_{t}.tif"
    if not path.exists():
        print(f"missing tile {t}: {len(grp)} fires skipped")
        continue
    with rasterio.open(path) as src:
        z = src.read(1).astype("float32")
        if src.nodata is not None:
            z[z == src.nodata] = np.nan
        z[z < -1000] = np.nan
        slope = slope_degrees(z, grp["LATITUDE"].mean())
        r, c = rasterio.transform.rowcol(src.transform, grp["LONGITUDE"].values, grp["LATITUDE"].values)
        r = np.clip(np.array(r), 0, z.shape[0] - 1)
        c = np.clip(np.array(c), 0, z.shape[1] - 1)
        fires.loc[grp.index, "elev_m"] = z[r, c] # raw 30 m cell at the fire
        fires.loc[grp.index, "slope_deg"] = slope[r, c] # slope from the ~1.5 km smoothed surface
        fires.loc[grp.index, "near_tile_edge"] = ((r < half) | (r >= z.shape[0] - half) |
                                                  (c < half) | (c >= z.shape[1] - half))
    print(f"{t}: {len(grp)} fires")


# ---- 4. Save ------------------------------------------------
fires = fires.drop(columns="tile")
print(fires["size_group3"].value_counts().reindex(["small", "medium", "large"]))
print(fires[["max_temp_c", "min_rh", "wind_speed", "erc", "elev_m", "slope_deg"]].describe().round(2))
print("missing elevation or slope:", int(fires[["elev_m", "slope_deg"]].isna().any(axis=1).sum()))

fires.to_csv(OUT_CSV, index=False)
