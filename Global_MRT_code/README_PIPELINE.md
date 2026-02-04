# Global MRT Pipeline - Quick Reference

## Directory Structure

```
~/
├── Datasets/Global_SOC/              # Raw observation data (one level up)
└── Global_MRT_code/                  # ← Pipeline code and outputs
    ├── 00_config_expanded.R          # Configuration (paths, params, packages)
    ├── 01_utils_expanded.R           # Utility functions (point + spatial)
    ├── verify_setup.R                # ← Run this first!
    ├── 03_extract_climate.R          # (to be created)
    ├── run_pipeline.R                # (to be created)
    │
    ├── outputs/                      # Final point dataset (tabular .rds)
    ├── spatialized_layers/           # Final raster layers (0.5° GeoTIFFs)
    │   ├── climate/
    │   ├── soilgrids/
    │   ├── topography/
    │   ├── landcover/
    │   ├── soilclass/
    │   └── productivity/
    │
    ├── temp/                         # Intermediates (safe to delete)
    ├── era5_data/                    # Cached ERA5 downloads
    ├── soilgrids/                    # Cached SoilGrids downloads
    ├── koppen_data/                  # Cached Köppen data
    └── hwsd_data/                    # Cached HWSD data
```

## First Steps

1. **Copy files to Global_MRT_code/**
   ```bash
   cd ~/Global_MRT_code
   # You should have:
   # - 00_config_expanded.R
   # - 01_utils_expanded.R  
   # - verify_setup.R
   ```

2. **Verify everything works**
   ```r
   source("verify_setup.R")
   ```
   This will:
   - Load config and utilities
   - Check that input data exists (../Datasets/Global_SOC/)
   - Create all output directories
   - Test that the 0.5° global grid can be created

3. **Once verified**, we create extraction modules:
   - `03_extract_climate.R` — Climate extraction + rasterization
   - (Then soil variables, landcover, etc.)
   - `run_pipeline.R` — Main orchestrator

## Key Features

**Point Workflow** (existing):
- Extract variables at observation point locations
- Store in tabular format (.rds)
- Used for model calibration

**Spatial Workflow** (new):
- All variables rasterized to 0.5° global grid
- Stored as compressed GeoTIFFs in `spatialized_layers/`
- Used for global prediction with trained model

**Both workflows run in parallel** from the same data sources.

## File Naming Convention

**Point data:**
- `01_data_filtered.rds`
- `02_with_climate.rds`
- `03_with_soilgrids.rds`
- etc.

**Raster layers:**
- `spatialized_layers/climate/mean_annual_temp_0.5deg.tif`
- `spatialized_layers/soilgrids/clay_0_20cm_0.5deg.tif`
- etc.

## Configuration Parameters

Key settings in `00_config_expanded.R`:

```r
PARAMS$spatial_resolution = 0.5          # degrees (lon/lat)
PARAMS$raster_datatype = "FLT4S"         # 32-bit float
PARAMS$target_depth = 20                 # cm for soil variables
PARAMS$climate_years = 1999:2020
```

## Next Steps

Once `verify_setup.R` passes:

1. Create climate extraction module
2. Test it runs for both point + spatial
3. Follow same pattern for other modules
4. Create `run_pipeline.R` to orchestrate all steps
