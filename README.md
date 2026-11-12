# Small Fires, Big Data: MANOVA vs. Six ANOVAs on California Wildfires

## Author: Anna Wadlow
## Submission Date: 11/11/2026

## Data

The project uses California wildfires discovered from 2000 to 2024, with weather and terrain joined to each fire. Three sources are used:

- **Fire records:** USDA Forest Service Fire Program Analysis fire-occurrence database (FPA FOD), release FPA_FOD_20260615. This is US federal government data and is in the public domain. Short, Karen C. 2026. "Spatial Wildfire Occurrence Data for the United States, 1992–2024 [FPA_FOD_20260615]." Forest Service Research Data Archive. https://doi.org/10.2737/RDS-2013-0009.7
- **Weather:** gridMET daily surface meteorological data (~4 km grid): maximum temperature, minimum relative humidity, wind speed, and Energy Release Component. The data is in the public domain, and the authors ask that it be cited. Data acquired 11/01/2026. Abatzoglou, J. T. 2013. "Development of Gridded Surface Meteorological Data for Ecological Applications and Modelling." *International Journal of Climatology* 33 (1): 121–31. https://doi.org/10.1002/joc.3413. https://www.climatologylab.org/gridmet.html
- **Terrain:** USGS 3D Elevation Program (3DEP) 1 arc-second (~30 m) elevation tiles, used for elevation and slope. USGS 3DEP products are public domain. U.S. Geological Survey. 2026. "About 3DEP Products & Services." https://www.usgs.gov/3d-elevation-program/about-3dep-products-services

The raw files (the FPA FOD database, the gridMET files, and the DEM tiles) are too large to store in GitHub but can be downloaded from the links above. The build script below downloads the gridMET and DEM files automatically. The FPA FOD database must be downloaded by hand.

## Folder Description

- `report/` contains the Quarto and PDF documents for the report.
- `wildfire_data_build.py` extracts the California fires, joins the weather and terrain data, and writes `fires_ca_weather_terrain_joint_0024.csv`.
- `wildfire_manova_clean.R` contains the R code for the analysis: Table 1 and Figure 1, the independence and assumption checks, and the MANOVA and six ANOVAs.

To reproduce the results, run the Python script first to build the data file, then run the R script.

## External Sources
External sources (including LLMs) were used for debugging, plot and table formatting, and to assist in writing Python code for `wildfire_data_build.py`. 
