
# In this script, we edit environmental variables in two ways: 
# 1. we set precipitation to zero when subzero temperatures and 
# 2. we scale variables. 

# We also calculate daily average biomass values. 
# Finally, we check correlation between snow, temperature and biomass to impute 
# samples with biomass = 0 

library(tidyverse)
library(scales)

rm(list = ls())

source("code/utils.R")

predictors_raw <- readRDS("data/env_predictors_raw.rds")
dat_raw <- readRDS("data/dat.rds")

env_variables <- c("avg_t2m",
                   "GDD5",
                   "tot_tp",
                   "tot_evavt",
                   "avg_sde", 
                   "avg_swvl1",
                   "HFI", 
                   "HII", 
                   "elevation")

env_labs <- c(
  GDD5      = "'GDD5'~'('~degree*C~')'",
  avg_t2m   = "'Temperature'~'('~degree*C~')'",
  avg_sde   = "'Snow depth'~'(m)'",
  avg_swvl1 = "'Soil moisture'~'(m'^3~m^-3*')'",
  tot_evavt = "'Evapotranspiration'~'(mwe)'",
  tot_tp    = "'Precipitation'~'(m)'"
)

# EDIT ENVIRONMENTAL VARIABLES -------------------------------------------------

predictors <- 
  predictors_raw |> 
  filter(!is.na(avg_sde)) |> 
  # GDD5 can not be lower than 0 
  mutate(GDD5 = if_else(GDD5 < 0, 0, GDD5), 
         jdate = map_dbl(date, ~lubridate::yday(.x)))

# CLIMATE DATA: 10-YEAR AVERAGES -----------------------------------------------

clim.data <-
  predictors |> 
  group_by(SITE, SITE_LATITUDE, SITE_LONGITUDE, jdate) |> 
  summarise(across(intersect(env_variables, colnames(predictors_raw)), mean)) |> 
  ungroup() |> 
  mutate(site_sde = mean(avg_sde), .by = SITE)

saveRDS(clim.data, "data/clim_data.rds")

p_climate <- 
  clim.data |>
  # Shift southern hemisphere 
  mutate(jdate = map2_dbl(jdate, SITE_LATITUDE, ~if_else(.y < 0, shift_s_date(.x), .x))) |> 
  pivot_longer(all_of(intersect(colnames(clim.data), env_variables[-which(env_variables == "GDD5")])), 
               names_to = "var", values_to = "value") |>
  mutate(var = factor(var, levels = env_variables)) |> 
  ggplot() +
  aes(x = jdate, y = value, color = SITE_LATITUDE, group = SITE) +
  geom_line(alpha = 0.5) +
  facet_wrap(~var, scales = "free_y", labeller = as_labeller(env_labs, label_parsed)) +
  scale_color_gradient2(
    midpoint = 0,
    low  = brewer_pal(palette = "PiYG")(11)[1],
    mid  = brewer_pal(palette = "PiYG")(11)[6],
    high = brewer_pal(palette = "PiYG")(11)[11]
  ) +
  labs(color = "Latitude", 
       x = "Day of the year") + 
  theme_bw() + 
  theme(axis.title.y = element_blank(), 
        panel.grid.minor = element_blank())

ggsave(p_climate, 
       filename = "plots/climate.pdf", 
       width = 183, 
       height = 100, 
       units = "mm")

# DAILY BIOMASS VALUES ---------------------------------------------------------

# Some additional filtering 
dat <- 
  dat_raw |>
  # Remove NEGreenland - snow data is missing 
  filter(SITE != "NEGreenland") |> 
  # Remove samples that are the only sample of a site-year (consider start and end)
  mutate(year_start = format(START_DATE, "%Y"), 
         year_end = format(COLL_DATE, "%Y")) |> 
  mutate(alone_start = n() == 1, .by = c(SITE, year_start)) |> 
  mutate(alone_end = n() == 1, .by = c(SITE, year_end)) |> 
  filter(!(alone_start & alone_end)) |> 
  dplyr::select(-starts_with("year_"), -starts_with("alone_")) |> 
  # Remove sites with less than four samples
  filter(n() > 3, .by = SITE) |> 
  # Calculate daily biomass
  mutate(n_days = as.numeric(COLL_DATE - START_DATE)) |> 
  mutate(DB_sample = Weight/n_days) |> 
  rename(SAMPLE_WEIGHT = Weight) |>
  mutate(date = map2(START_DATE, COLL_DATE,
                     ~lubridate::as_date(.x:(.y-1)))) |>
  unnest(date)

# SEARCH FOR THRESHOLD FOR ZERO BIOMASS ----------------------------------------

biomass <- readRDS("data/PREDICTED_BIOMASS.rds")

# Plot to decide threshold for temp and snow 
p_snow <- 
  dat |> 
  left_join(biomass) |> 
  left_join(predictors) |> 
  summarise(across(intersect(env_variables, colnames(predictors)), mean), 
            .by = c(SAMPLE_CODE, 
                    SAMPLE_WEIGHT, 
                    BIOMASS_PostMean
                    )) |>
  mutate(cold = if_else(avg_t2m < 5, T, F)) |> 
  ggplot() + 
  geom_point(aes(x = avg_sde, y = BIOMASS_PostMean, color = cold), alpha = 0.3, shape = 1) + 
  labs(color = "< 5°C", 
       x = "Average snow depth", 
       y = "Sample weight (g)") + 
  geom_vline(aes(xintercept = 0.1), 
             linetype = "dashed") +
  theme_bw(base_size = 7) + 
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.2), 
        legend.position = "inside", 
        legend.position.inside = c(0.8,0.8), 
        legend.box.background = element_rect(color = "black")) + 
  scale_y_log10()

ggsave(plot = p_snow,
       filename = "plots/zero_imputation.pdf", 
       width = 57, 
       heigh = 60, 
       units = "mm")

# IMPUTE ZERO BIOMASS ----------------------------------------------------------

# Impute zeros when weekly average sde > 0.1 and weekly average t2m < 5

imputed_zeros <- 
  predictors |>
  # Join with SITE_TYPE data 
  left_join(dat |> distinct(SITE, SITE_TYPE)) |> 
  # Remove already sampled days 
  anti_join(dat) |> 
  # Limit to active sampling year per site 
  mutate(year = format(date, "%Y")) |>
  inner_join(dat |>
               mutate(year = format(COLL_DATE, "%Y")) |>
               distinct(SITE, year),
             by = c("SITE", "year")) |> 
  group_by(SITE) |> 
  arrange(date) |> 
  # Group unsampled days into sample weeks, days of the sample week must be consecutive 
  mutate(gap = c(0, diff(date)), 
         group_id = cumsum(gap > 1)) |> 
  group_by(SITE, group_id) |> 
  mutate(sub_group = (row_number() - 1) %/% 7) |> 
  ungroup() |> 
  mutate(SAMPLE_CODE = paste(SITE, group_id, sub_group, sep = "_")) |> 
  group_by(SAMPLE_CODE, SITE, SITE_LONGITUDE, SITE_LATITUDE) |> 
  # Use only periods of at least 6 consecutive days 
  mutate(n_days = n()) |> 
  filter(n_days >= 6) |> 
  mutate(START_DATE = first(date), 
         COLL_DATE = last(date)) |> 
  group_by(SAMPLE_CODE,
           SITE,
           SITE_LONGITUDE,
           SITE_LATITUDE,
           SITE_TYPE,
           START_DATE,
           COLL_DATE,
           n_days) |>
  summarise(avg_t2m = mean(avg_t2m), 
            avg_sde = mean(avg_sde)) |>
  filter(!is.na(avg_sde) & avg_sde > 0.1 & avg_t2m < 5) |> 
  mutate(DB_sample = 0) |> 
  # Sync with dat 
  mutate(n_days = as.numeric(COLL_DATE - START_DATE)) |> 
  mutate(date = map2(START_DATE, COLL_DATE,
                     ~lubridate::as_date(.x:(.y-1)))) |>
  unnest(date)
  
dat_w_zeros <- bind_rows(dat, imputed_zeros |> dplyr::select(intersect(colnames(dat), colnames(imputed_zeros))))

# WEATHER DATA: ANOMALIES ------------------------------------------------------

weather.data <- 
  dat_w_zeros |> 
  mutate(jdate = map_dbl(date, ~lubridate::yday(.x))) |> 
  left_join(predictors_raw, by = c("SITE", "SITE_LATITUDE", "SITE_LONGITUDE", "date")) |> 
  left_join(clim.data, by = c("SITE", "SITE_LATITUDE", "SITE_LONGITUDE", "jdate"), 
            suffix = c("", ".clim")) |> 
  # Per day weather anomalies 
  mutate(avg_t2m.anom = avg_t2m - avg_t2m.clim,
         GDD5.anom = GDD5 - GDD5.clim, 
         avg_sde.anom = avg_sde - avg_sde.clim, 
         avg_swvl1.anom = avg_swvl1 - avg_swvl1.clim, 
         tot_evavt.anom = tot_evavt - tot_evavt.clim, 
         tot_tp.anom = tot_tp - tot_tp.clim) |> 
  dplyr::select(-avg_t2m, -avg_sde, -avg_swvl1, -tot_evavt, -tot_tp, -GDD5) |> 
  # Calculate means across sample period
  summarise(across(c(avg_t2m.anom, avg_t2m.clim,
                   GDD5.anom, GDD5.clim, 
                   avg_sde.anom, avg_sde.clim, 
                   avg_swvl1.anom, avg_swvl1.clim, 
                   tot_evavt.anom, tot_evavt.clim, 
                   tot_tp.anom, tot_tp.clim), mean), 
            jdate = median(jdate), 
            .by = c(SAMPLE_CODE, SAMPLE_WEIGHT, SITE, SITE_LATITUDE, SITE_LONGITUDE, SITE_TYPE, DB_sample, n_days, START_DATE, COLL_DATE, site_sde)) 


# ADD HUMAN FOOTPRINT INDEX AND ELEVATION --------------------------------------

hfi <- 
  read_csv("data/raw_data/LIFEPLAN_covariates_PN.csv") |> 
  dplyr::select(NAME, human_footprint_index, human_influence_index) |> 
  rename(SITE = NAME, 
         HFI = human_footprint_index, 
         HII = human_influence_index)

sites_w_elevation <- 
  readRDS("data/site_with_elevation.rds") |> 
  rename(elevation = altitude)

sample.data <- 
  weather.data |> 
  left_join(hfi) |> 
  left_join(sites_w_elevation) 

# CHECK VALUES OF ENVIRONMENTAL PREDICTORS -------------------------------------
# 
# # Plot predictors
# predictors_raw |>
#   pivot_longer(cols = all_of(env_variables), 
#                names_to = "variable", 
#                values_to = "value") |> 
#   ggplot() + 
#   geom_histogram(aes(x = value)) + 
#   facet_wrap(~variable, scales = "free") + 
#   scale_y_log10()
# 
# CHECK CORRELATION BETWEEN PREDICTORS -----------------------------------------
# 
# # Apply Z-score standardization: (X - mean) / SD
# predictors_scaled <- predictors
# predictors_scaled[env_variables] <- scale(predictors_scaled[env_variables])
# 
# # Check scaling 
# apply(predictors_scaled[env_variables], 2, function(x) mean(x, na.rm = TRUE))  # Should be close to 0
# apply(predictors_scaled[env_variables], 2, function(x) sd(x, na.rm = TRUE))    # Should be close to 1
# 
# # Check correlation 
# predictors_scaled |>
#   dplyr::select(all_of(env_variables), SITE_LATITUDE) |> 
#   mutate(SITE_LATITUDE = abs(SITE_LATITUDE)) |> 
#   slice_sample(n = 2000) |> 
#   psych::pairs.panels() 
  
# SAVE DATASETS WITH DAILY VALUES AND ENV VARIABLES ----------------------------

#min(sample.data |> filter(DB_sample > 0) |> pull(avg_t2m.clim))
tempcutoff_t2m <- -7

sample.data |> 
  # Truncate temperature at minimum observed value to not drive model with only zeros 
  # When truncating, set anomaly to 0
  # Motivation: colder than -7 should not make a difference, it is perceived as cold  
  mutate(avg_t2m.anom = if_else(avg_t2m.clim < tempcutoff_t2m, 0, avg_t2m.anom),
         avg_t2m.clim = if_else(avg_t2m.clim < tempcutoff_t2m, tempcutoff_t2m, avg_t2m.clim)) |> 
  saveRDS("data/sample_data_unscaled.rds")

sample.data |> 
  # Truncate temperature at minimum observed value to not drive model with only zeros 
  # When truncating, set anomaly to 0
  # Motivation: colder than -7 should not make a difference, it is perceived as cold  
  mutate(avg_t2m.anom = if_else(avg_t2m.clim < tempcutoff_t2m, 0, avg_t2m.anom),
         avg_t2m.clim = if_else(avg_t2m.clim < tempcutoff_t2m, tempcutoff_t2m, avg_t2m.clim)) |>
  mutate(across(c(avg_t2m.anom, 
                  avg_t2m.clim, 
                  GDD5.anom,
                  GDD5.clim,
                  avg_sde.anom,
                  avg_sde.clim,
                  avg_swvl1.anom,
                  avg_swvl1.clim,
                  tot_evavt.anom,
                  tot_evavt.clim,
                  tot_tp.anom,
                  tot_tp.clim,
                  HFI,
                  HII,
                  elevation, 
                  site_sde), 
                ~as.numeric(scale(.x)))) |> 
  saveRDS("data/sample_data_scaled.rds")
