
library(tidyverse)
library(kgc)

dat_raw <- readRDS("data/sample_data_scaled.rds")

dat <- 
  dat_raw |> 
  # mutate(rndCoord.lon = RoundCoordinates(SITE_LONGITUDE), 
  #        rndCoord.lat = RoundCoordinates(SITE_LATITUDE)) |> 
  dplyr::rename(Longitude = SITE_LONGITUDE, 
         Latitude = SITE_LATITUDE, 
         Site = SITE) |> 
  select(Site, Longitude, Latitude) |> 
  distinct()

dat$CZ <- as.character(kgc::LookupCZ(data = dat, rc = TRUE))

dat <- 
  dat |> 
  dplyr::mutate(
    CZ = dplyr::case_when(
      grepl("A.", CZ) ~ "Tropical-Subtropical", # tropical
      grepl("B.h", CZ) ~ "Tropical-Subtropical", # warm deserts
      grepl("B.k", CZ) ~ "Temperate", # cold deserts
      grepl("C[fw]a", CZ) ~ "Tropical-Subtropical", # subtropical
      grepl("C..", CZ) ~ "Temperate", # Mediterannean, oceanic, highland
      grepl("[DE]..?", CZ) ~ "Polar-Continental"
      #    Latitude > 50 ~ "Polar-Continental",
      #    abs(Latitude) < 30 ~ "Tropical-Subtropical",
      #    TRUE ~ "Temperate"
    )
  )

# Check site pairs with distinct CZ 
unmatched <- 
  dat |> 
  mutate(sitepair = map_chr(Site, ~str_split(.x, "_") |>
                              unlist() |> head(1))) |>
  distinct(sitepair, CZ) |>
  group_by(sitepair) |>
  filter(n() > 1) |>
  pull(sitepair) |> 
  unique()

dat |> 
  filter(str_detect(Site, paste(unmatched, collapse = "|")))

saveRDS(dat, 
        "data/CZ_sites.rds")
