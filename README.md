# biomass

Data and code for analyzing global patterns of insect biomass using Lifeplan data.

## Folder Structure

```
biomass/
├── data/
│   ├── raw_data/                  # Raw data of biomass, snow depth for sites near glaciers, and HFI
│   ├── dat.rds                    # Cleaned biomass data
│   ├── PREDICTED_BIOMASS.rds      # Estimated biomass data, removing ethanol weight
│   ├── env_predictors_raw.rds     # Extracted environmental covariates from ERA5
│   └── site_with_elevation.rds    # Extracted elevation from Copernicus DEM
│
├── code/
│   ├── utils.R                    # Helper functions
│   ├── 1_compile_data.R           # Compile climate, weather, elevation, and HFI data
│   ├── 2_climate_residuals.R      # Model seasonal variation in climate and calculate residual values
│   ├── 3_CZ                       # Extract climate zone for each site
│   └── 4_biomass_models.R         # Spatial and seasonal analyses
│
└── result/
    ├── season_model_eval.csv      # AIC, R2 for each step of the model of seasonal variation 
    └── spatial_model_eval.csv     # R2 for each step of the model of spatial variation
```

> Each script in `code/` builds on the previous, creating files necessary for the next step.
