
rm(list = ls())

# Load libraries ---------------------------------------------------------------

source("code/utils.R")

library(tidyverse)
library(lme4)
library(MuMIn)
library(ggpubr)
library(spdep)

# Data preparation -------------------------------------------------------------

## Read data ----
dat_raw <- readRDS("data/sample_data_scaled.rds")           # Biomass+site data 
clim_data <- readRDS("data/clim_data.rds")                  # Climate
climate.residuals <- readRDS("data/climate_residuals.rds")  # Climate deviations
biomass <- readRDS("data/PREDICTED_BIOMASS.rds")            # Predicted biomass

dat <- 
  left_join(dat_raw, biomass) |> 
  filter(n_days >= 3) |> 
  mutate(BIOMASS_low_included = map2_dbl(DB_sample, BIOMASS_PostMean,
                                        \(x, y) if(x == 0) runif(1, 0.01, min(biomass$BIOMASS_PostMean)) else y)) |>  # Sample low values instead of imputed zeros
  mutate(DB_post_sample = BIOMASS_low_included/n_days)  # Calculate daily biomass 

  
## Sample-level data for seasonal analyses ----
dat_sample <- list()

dat_sample$all <- 
  dat |> 
 # filter(!is.na(BIOMASS_PostMean)) |> 
  mutate(logDB_post_sample = log(DB_post_sample), 
         year = factor(format(COLL_DATE, "%Y")),
         SITE = factor(SITE), 
         lat = abs(SITE_LATITUDE/100)) |> 
  # Shift jdate for Southern hemisphere 
  mutate(jdate = map2_dbl(jdate, SITE_LATITUDE, 
                          ~if_else(.y < 0, round(shift_s_date(.x)), round(.x)))) |> 
  mutate(jdate = if_else(jdate == 0, 1, jdate)) |> 
  # Join with climate residual data 
  left_join(climate.residuals, 
            by = c("SITE", "SITE_LATITUDE", "SITE_LONGITUDE", "jdate")) |> 
  # Add seasonality terms 
  mutate(cos1 = cos(2 * pi * jdate / 365), 
         sin1 = sin(2 * pi * jdate / 365), 
         cos2 = cos(4 * pi * jdate / 365),
         sin2 = sin(4 * pi * jdate / 365)) |> 
  # Scale residuals
  mutate(across(c(temp.res,
                  gdd5.res, 
                  prec.res, 
                  et.res, 
                  snow.res, 
                  soil.res), 
                ~as.numeric(scale(.x)))) 

# Annual mean climate per site 
siteclim <- 
  clim_data |>
  group_by(SITE, SITE_LATITUDE, SITE_LONGITUDE) |> 
  summarise(across(c(avg_t2m, GDD5, tot_tp, tot_evavt, avg_sde, avg_swvl1), mean)) |> 
  ungroup() |> 
  mutate(across(c(avg_t2m, GDD5, tot_tp, tot_evavt, avg_sde, avg_swvl1), 
                ~as.numeric(scale(.x)))) |> 
  rename(avg_t2m.siteclim = avg_t2m,
         GDD5.siteclim = GDD5,
         tot_tp.siteclim = tot_tp,
         tot_evavt.siteclim = tot_evavt,
         avg_sde.siteclim = avg_sde,
         avg_swvl1.siteclim = avg_swvl1)

## Site-level data for spatial analyses ----
dat_site <- list()

dat_site$all <- 
  dat |> 
  mutate(year = factor(format(COLL_DATE, "%Y")),
         SITE = factor(SITE),
         lat = abs(SITE_LATITUDE/100)) |> 
  summarise(logDB_site = log(sum(BIOMASS_low_included, na.rm = T) / sum(n_days)),
            #logDB_site_old = log(sum(SAMPLE_WEIGHT, na.rm = T) / sum(n_days)), 
            .by = c(SITE, SITE_TYPE, HFI, HII, elevation, lat)) |>
  left_join(siteclim, by = "SITE") |> 
  mutate(country = case_when(
    str_detect(SITE_TYPE, "madagascar") ~ "MAD",
    str_detect(SITE_TYPE, "Nordic") ~ "NORD",
    TRUE ~ str_remove(SITE, "_[A-Za-z]+$")
    )) |>
  mutate(cluster = case_when(
    SITE_TYPE == "madagascar_hierarchical" ~ "MAD_hier",
    SITE_TYPE == "Nordic Hierarchical" ~ "NORD_hier",
    str_detect(SITE_TYPE, "madagascar") ~ "MAD",
    str_detect(SITE_TYPE, "Nordic") ~ "NORD",
    TRUE ~ str_remove(SITE, "_[A-Za-z]+$")
    ))

## Data subsets ----

# Random selection of one hierarchical site/cluster 
# (run this once and use the result)
keep_hierarchical <-
  dat_sample$all |>
  filter(SITE_TYPE == "madagascar_hierarchical" |
           SITE_TYPE == "Nordic Hierarchical") |>
  add_count(SITE, name = "n_samples") |>
  filter(n_samples > 100) |>
  distinct(SITE, SITE_TYPE) |>
  group_by(SITE_TYPE) |>
  sample_n(1) |>
  pull(SITE)

# Here 
#keep_hierarchical <- c("SH5", "MRT6")

dat_sample$no_hierarchical <- 
  dat_sample$all |> 
  filter(SITE_TYPE != "madagascar_hierarchical" &
           SITE_TYPE != "Nordic Hierarchical" | 
           SITE %in% keep_hierarchical) 

dat_site$no_hierarchical <- 
  dat_site$all |> 
  filter(SITE_TYPE != "madagascar_hierarchical" &
           SITE_TYPE != "Nordic Hierarchical" | 
           SITE %in% keep_hierarchical) 

# SEASONALITY ANALYSIS ---------------------------------------------------------

## I: Latitude-dependent seasonality ----

lmer_re <- "+ (1|SITE) + (1|year)"          # Random effect

### Null model -- only site variation ----

null_formula <- as.formula(paste0("logDB_post_sample ~ 1", lmer_re))

# Null model for all data sets 
m_null <- lapply(dat_sample, \(x) lmer(null_formula, data = x, REML = F))

### M1 -- Latitude-dependent seasonality ----

# Main effect of seasonality, and its interaction with latitude
m1_fe <- "logDB_post_sample ~ (cos1 + sin1) + lat:(cos1 + sin1)"
m1_formula <- as.formula(paste(m1_fe, lmer_re))

m1 <- lapply(dat_sample, \(x) lmer(m1_formula, data = x, REML = F))

### M2 -- Second periodic function ----

m2_fe <- "logDB_post_sample ~ (cos1 + sin1 + cos2 + sin2) + lat:(cos1 + sin1 + cos2 + sin2)"
m2_formula <- as.formula(paste(m2_fe, lmer_re))

m2 <- lapply(dat_sample, \(x) lmer(m2_formula, data = x, REML = F))

### Check model choice ----

lapply(names(dat_sample), \(x) AIC(m_null[[x]], m1[[x]], m2[[x]])) #AIC 

#R2
data.frame(
  null_R2m = sapply(m_null, function(m) r.squaredGLMM(m)[1]),
  null_R2c = sapply(m_null, function(m) r.squaredGLMM(m)[2]),
  m1_R2m = sapply(m1, function(m) r.squaredGLMM(m)[1]),
  m1_R2c = sapply(m1, function(m) r.squaredGLMM(m)[2]),
  m2_R2m = sapply(m2, function(m) r.squaredGLMM(m)[1]),
  m2_R2c = sapply(m2, function(m) r.squaredGLMM(m)[2])
)

## II: Effect of environmental predictors on seasonality -----------------------

### Climate ----

season <- "(cos1 + sin1 + cos2 + sin2)"

# Candidate variables (in order)
clim_cand_se <- 
  list(
    temp = "temp.res",                                     # Temperature
    temp.season = paste0("temp.res:", season),
    temp_nl = "I(temp.res ^ 2)",
    temp_nl.season = paste0("I(temp.res ^ 2):", season),
    temp.latseason = paste0("temp.res:lat:", season),
    prec = "prec.res",                                     # Precipitation
    prec.season = paste0("prec.res:", season),
    prec_nl = "I(prec.res ^ 2)",
    prec_nl.season = paste0("I(prec.res ^ 2):", season),
    prec.latseason = paste0("prec.res:lat:", season),
    et = "et.res",                                         # Evapotranspiration
    et.season = paste0("et.res:", season),
    et.latseason = paste0("et.res:lat:", season),
    snow = "snow.res",                                     # Snow depth
    snow.season = paste0("snow.res:", season),
    snow.latseason = paste0("snow.res:lat:", season),
    soil = "soil.res",                                     # Soil moisture
    soil.season = paste0("soil.res:", season),
    soil_nl = "I(soil.res ^ 2)",
    soil_nl.season = paste0("I(soil.res ^ 2):", season),
    soil.latseason = paste0("soil.res:lat:", season),
    elevation = "elevation",                               # Elevation
    elevation.season = paste0("elevation:", season),
    elevation.latseason = paste0("elevation:lat:", season)
  )


# Run forward selection across datasets
result_se_clim <- lapply(names(dat_sample), 
                         \(x) forward_selection(start_model = m2[[x]], 
                                                candidates = clim_cand_se, 
                                                aic_criterion = 2, 
                                                dat = dat_sample[[x]])
                         )
names(result_se_clim) <- names(dat_sample)

# Are the same predictors selected? 
all(result_se_clim$all$selected_vars == result_se_clim$no_hierarchical$selected_vars)

# Which variables were not retained? 
setdiff(names(clim_cand_se), result_se_clim$all$selected_vars)

# Save model for predictions 
saveRDS(result_se_clim$all$model, "data/best_clim_model.rds")

### Weather ----

# Candidates (weather, in order)
weather_cand <- list(
  avg_t2m.anom = "avg_t2m.anom",                          # Temperature
  avg_t2m.season = paste0("avg_t2m.anom:", season), 
  tot_tp.anom = "tot_tp.anom",                            # Precipitation
  tot_tp.season = paste0("tot_tp.anom:", season),
  tot_evavt.anom = "tot_evavt.anom",                      # Evapotranspiration
  tot_evavt.season = paste0("tot_evavt.anom:", season),
  avg_sde.anom = "avg_sde.anom",                          # Snow depth
  avg_sde.season = paste0("avg_sde.anom:", season),
  avg_swvl1.anom = "avg_swvl1.anom",                      # Soil moisture 
  avg_swvl1.season = paste0("avg_swvl1.anom:", season)
)

# Forward selection across datasets
result_se_weather <- 
  lapply(names(dat_sample), 
         \(x) forward_selection(start_model = result_se_clim[[x]]$model, 
                                candidates = weather_cand,
                                aic_criterion = 2, 
                                dat = dat_sample[[x]])
         )
names(result_se_weather) <- names(dat_sample)

all(result_se_weather$all$selected_vars == result_se_weather$no_hierarchical$selected_vars)

### Human footprint ----

human_cand_se <- c(hfi = "HFI", 
                hfi.season = paste0("HFI:", season), 
                hfi.latseason = paste0("HFI:lat:", season))

# Forward selection across datasets
result_se_human <- 
  lapply(names(dat_sample), 
         \(x) forward_selection(start_model = result_se_weather[[x]]$model, 
                                candidates = human_cand_se,
                                aic_criterion = 2, 
                                dat = dat_sample[[x]])
  )
names(result_se_human) <- names(dat_sample)

lapply(result_se_human, \(x) x$selected_vars)

## Re-fit final model with REML ----
se_reml <- lapply(names(result_se_human),
                  \(x) update(result_se_human[[x]]$model, 
                              data = dat_sample[[x]],
                              REML = T))
names(se_reml) <- names(result_se_human)

## Save seasonal result ----
seasonal <- list(data = dat_sample, 
                 se_lat = m2, 
                 se_clim = result_se_clim, 
                 se_weather = result_se_weather, 
                 se_human = result_se_human, 
                 se_reml = se_reml)


# saveRDS(list(data = dat_sample, 
#              se_lat = m2, 
#              se_clim = result_se_clim, 
#              se_weather = result_se_weather, 
#              se_human = result_se_human, 
#              se_reml = se_reml), 
#         file = "results/se_results.rds",
#         compress = FALSE)

# SPATIAL ANALYSIS -------------------------------------------------------------

## Get biomass peak data ----

# Create prediction grid
preddat <- 
  expand.grid(jdate = seq(1, 365), 
              SITE = unique(dat_site$all$SITE)) |> 
  left_join(dat_site$all |> distinct(SITE, lat)) |> 
  mutate(cos1 = cos(2 * pi * jdate / 365),
         sin1 = sin(2 * pi * jdate / 365),
         cos2 = cos(4 * pi * jdate / 365),
         sin2 = sin(4 * pi * jdate / 365), 
         year = 2023)

# Predict biomass 
preddat$pred <- predict(m2$all, newdata = preddat, re.form = ~(1|SITE))

# Find three-week period with maximum biomass
peak_period <- 
  preddat |>
  group_by(SITE) |>
  arrange(jdate) |> 
  mutate(
    roll_mean = zoo::rollapply(pred, 
                               width = 21, 
                               FUN = mean,
                               align = "left",
                               fill = NA)
  ) |>
  slice_max(roll_mean, n = 1, with_ties = FALSE) |>
  mutate(
    start_day = jdate,
    end_day = jdate + 21 - 1,
    end_day = ifelse(end_day > 365, 365, end_day)
  ) |>
  dplyr::select(SITE, start_day, end_day, roll_mean)

# Only years and site where the peak has been sampled should be included 
peak_sampled_year <- 
  dat_sample$all |> 
  dplyr::select(SITE, year, jdate) |> 
  left_join(# Add one week before and after
    peak_period |> 
      mutate(jdate = map2(start_day, end_day, ~(.x-14):(.y+14))) |> 
      unnest(jdate) |> 
      dplyr::select(SITE, jdate) |> 
      mutate(jdate = if_else(jdate < 1, jdate + 365, 
                             if_else(jdate > 365, jdate - 365, jdate)), 
             peak_sampled = T)) |> 
  filter(peak_sampled) |> 
  distinct(SITE, year) 

dat_days <-
  dat_sample$all |> 
  # Filter so that only years where the peak was feasibly sampled is included 
  inner_join(peak_sampled_year) |> 
  mutate(sampled_days = map2(START_DATE, COLL_DATE, ~(lubridate::yday(.x) + 1):lubridate::yday(.y))) |> 
  unnest(sampled_days) |> 
  dplyr::select(SITE, DB_post_sample, year, sampled_days) |> 
  group_by(SITE, year) |> 
  arrange(sampled_days, .by_group = T)

max_window_length <- 21
min_required_days <- 19

peaks <- 
  dat_days |>
  arrange(SITE, year, sampled_days) |>
  group_by(SITE, year) |>
  # For each site-year, slide a 21-day window along sampled_days
  mutate(
    total_3wk = slider::slide_dbl(
      .x = seq_along(sampled_days),
      .f = ~{
        start_day <- sampled_days[.x[1]]
        end_day   <- start_day + max_window_length - 1
        # subset days in this range
        window_data <- DB_post_sample[sampled_days >= start_day & sampled_days <= end_day]
        n_days <- length(window_data)
        if (n_days >= min_required_days) sum(window_data) else NA_real_
      },
      .before = 0, .after = 0  # not rolling by index, custom inside
    )
  ) |> 
  # Identify the window with maximum biomass
  slice_max(total_3wk, n = 1, with_ties = FALSE) |>
  mutate(
    start_day = sampled_days,
    end_day   = start_day + max_window_length - 1
  ) |>
  dplyr::select(SITE, year, start_day, end_day, total_3wk) 

dat_peak <- 
  lapply(dat_site, \(x) inner_join(x,
                                   peaks |>
                                     ungroup() |> 
                                     summarise(log_DB_peak = mean(log(total_3wk/21),
                                                                  na.rm = T),
                                               .by = SITE)))
## Latitudinal gradient ----

# Mean daily biomass 
lm_lat <- lapply(dat_site, \(x) lm(logDB_site ~ lat, data = x))
lapply(lm_lat, summary)

lm_lat_peak <- lapply(dat_peak, \(x) lm(log_DB_peak ~ lat, data = x))
lapply(lm_lat_peak, summary)

## Effect of climate -----------------------------------------------------------

### Spatial structure ----

sites_sf <- st_as_sf(dat_site$all,
                     coords = c("SITE_LONGITUDE", "SITE_LATITUDE"),
                     crs = 4326)
sites_sf_m <- st_transform(sites_sf, 3857)  # Project for m-based dist 
coords <- st_coordinates(sites_sf_m)

dat_site$all$in_nord <- as.numeric(dat_site$all$cluster == "NORD_hier") 
dat_site$all$in_mad <- as.numeric(dat_site$all$cluster == "MAD_hier")
dat_site$all$X <- coords[,1]
dat_site$all$Y <- coords[,2]

dat_peak$all <- 
  inner_join(dat_peak$all,
             dat_site$all)

### Climate ----

clim_cand_sp <- 
  list(
    avg_t2m.siteclim = "avg_t2m.siteclim",
    nl_avg_t2m.siteclim = "I(avg_t2m.siteclim^2)",
    tot_tp.siteclim = "tot_tp.siteclim", 
    nl_tot_tp.siteclim = "I(tot_tp.siteclim^2)", 
    tot_evavt.siteclim = "tot_evavt.siteclim", 
    avg_sde.siteclim = "avg_sde.siteclim", 
    avg_swvl1.siteclim = "avg_swvl1.siteclim",
    nl_avg_swvl1.siteclim = "I(avg_swvl1.siteclim^2)", 
    elevation = "elevation",
    latitude = "lat"
)

#### Mean daily biomass -----

m_mean <- 
  mgcv::gam(
    logDB_site ~ 1 +
      s(X, Y, bs = "gp", by = in_nord, k = 15) +   # smooth only within that cluster
      s(X, Y, by = in_mad, k = 28),
    family = gaussian(),
    data = dat_site$all,
    method = "ML"
)

result_sp_clim <- 
  forward_selection_gam(start_model = m_mean, 
                        candidates = clim_cand_sp, 
                        aic_criterion = 2, 
                        dat = dat_site$all)

# Re-fit with REML 
spat_m <- mgcv::gam(
  formula(result_sp_clim$model),
  family = gaussian(),
  data = dat_site$all,
  method = "REML"
)

##### Check residual spatial correlation ----

spat_res <- residuals(spat_m)

moran.mc(spat_res[dat_site$all$country == "NORD"],   # Nordic sites 
           build_neighbour(sf = sites_sf_m, 
                           dist = 10000, 
                           clustcol = "country",
                           clustval = "NORD"), 
         nsim = 999)

moran.mc(spat_res[dat_site$all$country == "MAD"],   # Madagascar sites 
           build_neighbour(sf = sites_sf_m, 
                           dist = 10000, 
                           clustcol = "country",
                           clustval = "MAD"),
           nsim = 999)

moran.mc(spat_res[dat_site$all$cluster == "NORD_hier"],   # Nordic hier 
           build_neighbour(sf = sites_sf_m, 
                           dist = 10000, 
                           clustcol = "cluster",
                           clustval = "NORD_hier"),
           nsim = 999)

moran.mc(spat_res[dat_site$all$cluster == "MAD_hier"],   # Madagascar hier 
           build_neighbour(sf = sites_sf_m, 
                           dist = 10000, 
                           clustcol = "cluster",
                           clustval = "MAD_hier"),
           nsim = 999)

#### Peak biomass -----

m_peak <- 
  mgcv::gam(
    log_DB_peak ~ 1 +
      s(X, Y, bs = "gp", by = in_nord, k = 15) +   # smooth only within that cluster
      s(X, Y, by = in_mad, k = 28),
    family = gaussian(),
    data = dat_peak$all,
    method = "ML"
  )

result_sp_peak <- 
  forward_selection_gam(start_model = m_peak, 
                        candidates = clim_cand_sp, 
                        aic_criterion = 2, 
                        dat = dat_peak$all)

# Re-fit with REML 
spat_p <- mgcv::gam(
  formula(result_sp_peak$model),
  family = gaussian(),
  data = dat_peak$all,
  method = "REML"
)

##### Check residual spatial correlation ----

peak_res <- residuals(spat_p)

moran.mc(peak_res[dat_peak$all$country == "NORD"],   # Nordic sites 
           build_neighbour(sf = sites_sf_m |> inner_join(dat_peak$all), 
                           dist = 10000, 
                           clustcol = "country",
                           clustval = "NORD"), 
         nsim = 999)

moran.mc(peak_res[dat_peak$all$country == "MAD"],   # Madagascar sites 
           build_neighbour(sf = sites_sf_m |> inner_join(dat_peak$all), 
                           dist = 10000, 
                           clustcol = "country",
                           clustval = "MAD"), 
           nsim = 999)

moran.mc(peak_res[dat_peak$all$cluster == "NORD_hier"],   # Nordic hier 
           build_neighbour(sf = sites_sf_m |> inner_join(dat_peak$all), 
                           dist = 10000, 
                           clustcol = "cluster",
                           clustval = "NORD_hier"), 
         nsim = 999)

moran.mc(peak_res[dat_peak$all$cluster == "MAD_hier"],   # Madagascar hier 
           build_neighbour(sf = sites_sf_m |> inner_join(dat_peak$all), 
                           dist = 10000, 
                           clustcol = "cluster",
                           clustval = "MAD_hier"), 
         nsim = 999)

### Human footprint ----


human_cand_sp <- c(hfi = "HFI", 
                   hfilat = "lat*HFI")

result_sp_human <- 
  forward_selection_gam(start_model = result_sp_clim$model, 
                        candidates = human_cand_sp, 
                        aic_criterion = 2, 
                        dat = dat_site$all)

result_sp_human_peak <- 
  forward_selection_gam(start_model = result_sp_peak$model, 
                        candidates = human_cand_sp, 
                        aic_criterion = 2, 
                        dat = dat_peak$all)


spat_m_human <-
  mgcv::gam(
    formula(result_sp_human$model),
    family = gaussian(),
    data = dat_site$all,
    method = "REML"
  )

### Natural vs. Urban ----

naturb_data <- 
  dat_site$all |> 
  rename(log_DB_site = logDB_site) |> 
  left_join(dat_peak$all |> dplyr::select(SITE, log_DB_peak)) |>
  separate(SITE, into = c("SITEPAIR", "NATURB"), sep = "_") |> 
  filter(!is.na(NATURB)) |> 
  filter(n() == 2, .by = SITEPAIR) |> 
  dplyr::select(SITEPAIR, NATURB, log_DB_site, log_DB_peak) |> 
  pivot_wider(names_from = NATURB, values_from = c(log_DB_site, log_DB_peak), 
              names_glue = "{sub('log_DB_', '', .value)}_{NATURB}") 

t.test(naturb_data$site_Natural, naturb_data$site_Urban, paired = TRUE)
t.test(naturb_data$peak_Natural, naturb_data$peak_Urban, paired = TRUE)

# Save spatial result ---- 

spatial <- list(data = dat_site, 
                data_peak = dat_peak, 
                data_naturb = naturb_data, 
                mean_lat = lm_lat, 
                peak_lat = lm_lat_peak, 
                mean_clim = result_sp_clim, 
                mean_human = result_sp_human, 
                peak_clim = result_sp_peak, 
                peak_human = result_sp_human_peak, 
                sp_reml = spat_m_human)

# saveRDS(list(data = dat_site, 
#              data_peak = dat_peak, 
#              data_naturb = naturb_data, 
#              mean_lat = lm_lat, 
#              peak_lat = lm_lat_peak, 
#              mean_clim = result_sp_clim, 
#              mean_human = result_sp_human, 
#              peak_clim = result_sp_peak, 
#              peak_human = result_sp_human_peak, 
#              sp_reml = spat_m_human), 
#         "results/sp_results.rds", 
#         compress = FALSE)
