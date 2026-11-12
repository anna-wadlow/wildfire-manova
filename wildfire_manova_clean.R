# ============================================================
# Wildfire size groups: one-way MANOVA vs six one-way ANOVAs
# Groups: small (A-C), medium (D-E), large (F-G)
# Responses: slope, ERC, min RH, max temp, elevation, wind speed
# ============================================================

# ---- 0. Setup ----------------------------------------------
library(tidyverse)
set.seed(257)

DATA_FILE <- "fires_ca_weather_terrain_joint_0024.csv"

V     <- c("slope_deg", "erc", "min_rh", "max_temp_c", "elev_m", "wind_speed")
sizes <- c("small", "medium", "large")

# ---- 1. Data -----------------------------------------------
# raw = original units, fires = transformed copy used for all the tests
raw <- read.csv(DATA_FILE) |>
  select(size_group3, LATITUDE, LONGITUDE, weather_date, all_of(V)) |>
  mutate(size = factor(size_group3, levels = sizes)) |>
  drop_na(all_of(V))

fires <- raw |>
  mutate(slope_deg = log1p(slope_deg),
         wind_speed = log1p(wind_speed),
         min_rh = sqrt(min_rh),
         elev_m = sqrt(pmax(elev_m, 0))) # a few elevations are below 0

table(fires$size)

# Table 1: group means in original units
table1 <- raw |>
  group_by(Group = size) |>
  summarize(n = n(), across(all_of(V), mean))
table1

# Figure 1: distributions after transforming
strip_labels <- c(slope_deg = "log(1 + slope)", erc = "ERC", min_rh = "sqrt(min RH)",
                  max_temp_c = "Max temp (°C)", elev_m = "sqrt(elevation)",
                  wind_speed = "log(1 + wind)")
n_by <- table(raw$size)
size_labels <- paste0(c("Small", "Medium", "Large"), " (n = ",
                      formatC(as.integer(n_by), big.mark = ",", format = "d"), ")")

dens_long <- fires |>
  mutate(size = factor(size, levels = sizes, labels = size_labels)) |>
  pivot_longer(all_of(V), names_to = "measure", values_to = "values") |>
  mutate(measure = factor(measure, levels = V, labels = strip_labels[V]))

cols <- c("#6baed6", "#2171b5", "#08306b")
ggplot(dens_long, aes(values, colour = size, fill = size)) +
  geom_density(alpha = 0.15, linewidth = 0.8) +
  scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
  scale_x_continuous(n.breaks = 5) +
  facet_wrap(~ measure, scales = "free", ncol = 3) +
  labs(x = NULL, y = NULL, fill = NULL, colour = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", axis.text.y = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        panel.spacing.x = unit(1.5, "lines"),
        strip.text = element_text(face = "bold", hjust = 0))

# ---- 2. Independence Checks  ---------------------------------------
# fires in the same gridMET cell (1/24 degree) on the same date share identical weather
fires <- fires |>
  mutate(cell = paste(round(LATITUDE * 24), round(LONGITUDE * 24)),
         key  = paste(cell, weather_date))

fires$shared <- duplicated(fires$key) | duplicated(fires$key, fromLast = TRUE)
mean(fires$shared) # overall share
round(100 * tapply(fires$shared, fires$size, mean), 1) # by group

# cell-dates that contain more than one size group
n_sizes <- tapply(as.character(fires$size), fires$key, n_distinct)
sum(n_sizes > 1)
mixed_keys <- names(n_sizes)[n_sizes > 1]
round(100 * tapply(fires$key %in% mixed_keys, fires$size, mean), 1) # % of each group's fires

# fires with exactly the same coordinates (same terrain values)
xy <- paste(fires$LATITUDE, fires$LONGITUDE)
same_xy <- duplicated(xy) | duplicated(xy, fromLast = TRUE)
round(100 * mean(same_xy), 1)
round(100 * tapply(same_xy, fires$size, mean), 1)

# ---- 3. Assumption checks (on the transformed data) --------
p <- length(V)

# equal variance / covariance
sds <- aggregate(fires[, V], by = list(size = fires$size), FUN = sd)
round(apply(sds[, V], 2, max) / apply(sds[, V], 2, min), 2) # SD ratios, worry above ~2

cors <- lapply(split(fires[, V], fires$size), cor)
round(max(abs(cors$small - cors$medium),
          abs(cors$small - cors$large),
          abs(cors$medium - cors$large)), 2) # biggest correlation gap between groups

# normality: skewness, then squared Mahalanobis distances vs the chi-square cutoff
skew <- function(x) mean((x - mean(x))^3) / sd(x)^3
skew_tab <- t(sapply(split(fires[, V], fires$size), function(d) sapply(d, skew)))
round(skew_tab, 2)
round(apply(abs(skew_tab), 2, max), 2)

cutoff <- qchisq(0.999, df = p)
outlier_rate <- function(d) {
  sapply(split(d[, V], d$size), function(x) {
    x <- as.matrix(x)
    100 * mean(mahalanobis(x, colMeans(x), cov(x)) > cutoff)  # % beyond the 99.9% cutoff, expect 0.1
  })
}
round(outlier_rate(fires), 2) # after transforming
round(outlier_rate(raw), 2) # before transforming

# multicollinearity
R <- cor(fires[, V])
round(R, 2)
pairs <- which(upper.tri(R), arr.ind = TRUE)
pair_tab <- data.frame(var1 = V[pairs[, 1]], var2 = V[pairs[, 2]], r = round(R[pairs], 2))
pair_tab[order(-abs(pair_tab$r)), ] # strongest pairs first
sum(abs(pair_tab$r) >= 0.2) # pairs with |r| >= 0.2
round(max(diag(solve(R))), 2) # biggest VIF
round(kappa(R, exact = TRUE), 1) # condition number

# ---- 4. Full-data MANOVA and the six ANOVAs ----------------
fit <- manova(as.matrix(fires[, V]) ~ size, data = fires)

fit_sum <- summary(fit, test = "Pillai")
fit_sum
summary(fit, test = "Wilks")

s <- min(p, nlevels(fires$size) - 1)
fit_sum$stats[1, "Pillai"] / s # V/s, average share explained

# the six ANOVAs: F, p, and share of variance explained (eta squared)
aov_list <- summary.aov(fit)
names(aov_list) <- V
anova_results <- map_dfr(aov_list, .id = "variable", \(a)
  tibble(F = a[1, "F value"],
         p = a[1, "Pr(>F)"],
         eta2 = a[1, "Sum Sq"] / sum(a[, "Sum Sq"])))
anova_results

all(anova_results$p < 0.05 / p) # all pass the Bonferroni adjusted cutoff?

# B and W matrices: diagonals are the six ANOVA sums of squares, off-diagonals are cross-products
B <- summary(fit)$SS$size
W <- summary(fit)$SS$Residuals
round(diag(B)); round(diag(W))
round(cov2cor(W), 2) # within-group correlations
