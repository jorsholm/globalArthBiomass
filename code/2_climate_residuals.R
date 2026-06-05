
library(tidyverse)

rm(list = ls())

source("code/utils.R")

climdat_raw <- readRDS("data/env_predictors_raw.rds")
cleandat <- readRDS("data/sample_data_scaled.rds")

climdat <- 
  climdat_raw |> 
  filter(SITE %in% unique(cleandat$SITE)) |> 
  mutate(jdate = lubridate::yday(date)) |> 
  filter(jdate != 366) |> 
  summarise(across(c("GDD5",
                     "avg_t2m",
                     "avg_sde",
                     "avg_swvl1",
                     "tot_evavt",
                     "tot_tp"), 
            \(x) mean(x, na.rm = T)), 
            .by = c(SITE, SITE_LONGITUDE, SITE_LATITUDE, jdate)) |> 
  # Shift southern hemisphere 
  mutate(jdate = map2_dbl(jdate, SITE_LATITUDE, ~if_else(.y < 0, shift_s_date(.x), .x))) |> 
  mutate(lat = abs(SITE_LATITUDE / 100),
         cos1 = cos(2 * pi * jdate / 365), 
         sin1 = sin(2 * pi * jdate / 365), 
         cos2 = cos(4 * pi * jdate / 365),
         sin2 = sin(4 * pi * jdate / 365))

# Build seasonal models for each variable --------------------------------------

season <- "(cos1 + sin1 + cos2 + sin2)"

## Temperature ---- 

temp <- lm(as.formula(paste0("avg_t2m ~ lat * ", season)), data = climdat)

## Prec ----

prec <- lm(as.formula(paste0("tot_tp ~ lat * ", season)), data = climdat)

## Snow ----

snow <- lm(as.formula(paste0("avg_sde ~ lat * ", season)), data = climdat) 

## ET ----

et <- lm(as.formula(paste0("tot_evavt ~ lat * ", season)), data = climdat) 

## Soil ----

soil <- lm(as.formula(paste0("avg_swvl1 ~ lat * ", season)), data = climdat) 

## GDD5 ---- 

gdd5 <- lm(as.formula(paste0("GDD5 ~ lat * ", season)), data = climdat)

## Save models for predictions 

saveRDS(list("temp" = temp, 
             "prec" = prec, 
             "snow" = snow, 
             "et" = et, 
             "soil" = soil), 
        "data/climate_seasonality_models.rds")

# Extract the residuals --------------------------------------------------------

climate.residuals <- 
  climdat |> 
  # Add predictions of variables capped at zero
  mutate(prec.pred = predict(prec), 
         et.pred = predict(et), 
         snow.pred = predict(snow), 
         soil.pred = predict(soil),
         # Add residuals 
         temp.res = resid(temp), 
         gdd5.res = resid(gdd5), 
         prec.res.raw = resid(prec), 
         et.res.raw = resid(et), 
         snow.res.raw = resid(snow), 
         soil.res.raw = resid(soil)) |> 
  # Modify residuals as if predictions were capped at zero 
  mutate(prec.res = if_else(prec.pred < 0, (prec.res.raw + prec.pred), prec.res.raw),
         et.res = if_else(et.pred < 0, (et.res.raw + et.pred), et.res.raw),
         snow.res = if_else(snow.pred < 0, (snow.res.raw + snow.pred), snow.res.raw),
         soil.res = if_else(soil.pred < 0, (soil.res.raw + soil.pred), soil.res.raw)) |> 
  dplyr::select(starts_with("SITE"), 
         jdate, 
         ends_with(".res"))


# Get temperature and snow cutoff 
sampled_clim <- 
  climate.residuals |> 
  inner_join(cleandat |> 
               filter(DB_sample != 0) |> 
               mutate(sampleday = map2(START_DATE, COLL_DATE, ~as.Date(.x:.y))) |> 
               unnest(sampleday) |> 
               mutate(jdate = map_dbl(sampleday, ~lubridate::yday(.x))) |> 
               mutate(jdate = map2_dbl(jdate, SITE_LATITUDE, ~if_else(.y < 0, shift_s_date(.x), .x))) |> 
               distinct(SITE, jdate))

sampled_clim$temp.res |> min()
sampled_clim$snow.res |> max()

min.temp <- -13.5
max.snow <- 0.95

climate.residuals.limited <- 
  climate.residuals |> 
  mutate(temp.res = if_else(temp.res < min.temp, min.temp, temp.res), 
         snow.res = if_else(snow.res > max.snow, max.snow, snow.res)) 

saveRDS(climate.residuals.limited, 
        "data/climate_residuals.rds")

# Predictions ------------------------------------------------------------------

jdate_pred <- seq(1, 365, by = 1)

CZ_sites <- readRDS("data/CZ_sites.rds") |> rename(SITE = Site)
CZ_lat <- 
  cleandat |> 
  left_join(CZ_sites) |> 
  group_by(CZ) |> 
  summarise(lat = mean(abs(SITE_LATITUDE)/100)) 

#lat_groups <- c(0.6, 0.3, 0)

# Create a data frame for predictions
pred_lat <- 
  expand.grid(jdate = jdate_pred,
              lat = CZ_lat$lat
  ) |> 
  mutate(cos1 = cos(2 * pi * jdate / 365),
         sin1 = sin(2 * pi * jdate / 365),
         cos2 = cos(4 * pi * jdate / 365),
         sin2 = sin(4 * pi * jdate / 365)) |> 
  left_join(CZ_lat)

pretty_names <- c(
  GDD5      = "'GDD5'~'('*degree*C*')'",
  avg_t2m   = "'Temperature'~'('*degree*C*')'",
  avg_sde   = "'Snow depth'~'(m)'",
  avg_swvl1 = "atop('Soil moisture','(m'^3~m^-3*')')",
  tot_evavt = "atop('Evapotranspiration','(mwe)')",
  tot_tp    = "'Precipitation'~'(m)'", 
  "0" = "0", 
  "0.3" = "30", 
  "0.6" = "60", 
  "Polar-Continental" = "Polar-Continental", 
  "Temperate" = "Temperate",
  "Tropical-Subtropical" = "Tropical-Subtropical"
)

p_climmodels <- 
  pred_lat |> 
  mutate(avg_t2m = predict(temp, newdata = pred_lat, allow.new.levels = T), 
         tot_tp = predict(prec, newdata = pred_lat, allow.new.levels = T), 
         tot_evavt = predict(et, newdata = pred_lat, allow.new.levels = T), 
         avg_swvl1 = predict(soil, newdata = pred_lat, allow.new.levels = T), 
         avg_sde = predict(snow, newdata = pred_lat, allow.new.levels = T), 
         #GDD5 = predict(gdd5, newdata = pred_lat, allow.new.levels = T)
         ) |> 
  pivot_longer(cols = c(#"GDD5",
                        "avg_t2m",
                        "avg_sde",
                        "avg_swvl1",
                        "tot_evavt",
                        "tot_tp"), 
               names_to = "var", 
               values_to = "val") |> 
  ggplot() +
  geom_point(data = climdat[sample(1:nrow(climdat), 20000),] |> 
               left_join(CZ_sites) |> 
               #mutate(lat_group = as.character(sapply(lat, function(x) lat_groups[which.min(abs(x - lat_groups))]))) |> 
               pivot_longer(cols = c(#"GDD5",
                                     "avg_t2m",
                                     "avg_sde",
                                     "avg_swvl1",
                                     "tot_evavt",
                                     "tot_tp"), 
                            names_to = "var", 
                            values_to = "val"),  
             aes(x = jdate,
                 y = val),
             color = "gray20",
             alpha = 0.1, 
             size = 0.2)  +
  geom_line(aes(x = jdate, 
                y = val#,
                #color = CZ
            ), 
            linewidth = 1, 
            color = "#d55e00") +
  theme_bw(base_size = 7) + 
  theme(legend.title = element_blank(), 
        panel.grid = element_blank(), 
        panel.spacing = unit(0, "lines"), 
        legend.position = "None") +
    # scale_color_manual(
    #   values = c(
    #     "Tropical-Subtropical" = "#009E73",  #"#009E73", 
    #     "Temperate"            = "#E69F00",  #F0E442", 
    #     "Polar-Continental"    = "#CC79A7"  #"#D55E00"  
    #   )
    # ) + 
  # scale_color_manual(
  #   values = c(
  #     "0" = "#009E73",  # green-teal (south)
  #     "0.3"    = "#F0E442",  # yellow (equator)
  #     # "0.3"  = "#E69F00",  # orange (mid-north)
  #     "0.6"  = "#D55E00"   # dark orange-red (far north)
  #   ),
  #   name = "Latitude"#, 
  #   #labels = c("0", "30", "60")
  # ) + 
  facet_grid(var~CZ, scales = "free", labeller = as_labeller(pretty_names, label_parsed)) + 
  labs(y = "Climate variable value", 
       x = "Day of the year")

ggsave(plot = p_climmodels, 
       "plots/climate_models.pdf", 
       width = 121, 
       height = 121, 
       units = "mm")


# Conceptual ----- 


p_climres <- 
  pred_lat |> 
  mutate(avg_t2m = predict(temp, newdata = pred_lat, allow.new.levels = T)) |> 
  filter(CZ == "Temperate") |> 
  select(jdate, predtemp = avg_t2m) |> 
  left_join(climdat |> 
              filter(SITE == "Perth_Natural"  & jdate %% 25 == 0) |> 
              select(jdate, mestemp = avg_t2m)) |> 
  ggplot() + 
  geom_segment(aes(xend = jdate, 
                   yend = predtemp, 
                   x = jdate, 
                   y = mestemp), 
               color = "firebrick") + 
  geom_line(aes(x = jdate, 
                y = predtemp), 
            linetype = "dashed") + 
  geom_point(aes(x = jdate, 
                 y = mestemp), 
             shape = 4) + 
  theme_classic() + 
  labs(x = "DOY", 
       y = "Climate") + 
  scale_x_continuous(breaks = c(1, 365)) + 
  theme(axis.text.y = element_blank(), 
        axis.ticks = element_blank())
  
ggsave(p_climres, 
       filename = "plots/conc_climres.pdf", 
       width = 40, 
       height = 30, 
       units = "mm")

p_wres <- 
  climdat |> 
  filter(SITE == "Perth_Natural" & jdate %% 5 == 0) |> 
  select(jdate, mestemp = avg_t2m) |> 
  left_join(climdat_raw |> filter(SITE == "Perth_Natural" & format(date, "%Y") == 2022) |> 
              mutate(jdate = lubridate::yday(date)) |> 
              filter(jdate %% 25 == 0) |> 
              select(jdate, avg_t2m)) |> 
  ggplot() + 
  geom_segment(aes(xend = jdate, 
                   yend = avg_t2m, 
                   x = jdate, 
                   y = mestemp), 
               color = "steelblue") + 
  geom_line(aes(x = jdate, 
                y = mestemp)) + 
  geom_point(aes(x = jdate, 
                 y = avg_t2m)) + 
  theme_classic() + 
  labs(x = "DOY", 
       y = "Weather") + 
  scale_x_continuous(breaks = c(1, 365)) + 
  theme(axis.text.y = element_blank(), 
        axis.ticks = element_blank())

ggsave(p_wres, 
       filename = "plots/conc_wres.pdf", 
       width = 40, 
       height = 30, 
       units = "mm")
