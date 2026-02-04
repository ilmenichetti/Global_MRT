# Step 3: Climate Extraction Module

## Overview

This module extracts ERA5-Land climate data at two levels:

1. **Point-based extraction** → Add columns to observation dataset
2. **Spatial rasterization** → 0.5° global GeoTIFFs for global prediction

## Variables Extracted

| Category | Variable | Aggregation | Output |
|----------|----------|-------------|--------|
| **Temperature** | 2m air temperature | Annual mean | `temperature_2m_mean_K` |
| | | Seasonality (std dev) | `temperature_seasonality_K` |
| **Precipitation** | Total precipitation | Annual mean | `precipitation_annual_m` |
| | | Seasonality (std dev) | `precipitation_seasonality_m` |
| | Solid precipitation (snow) | Annual mean | `solid_precipitation_m` |
| **Soil (0-20cm)** | Soil temperature | Annual mean (depth-weighted) | `soil_temperature_0_20cm_K` |
| | Soil moisture | Annual mean (depth-weighted) | `soil_moisture_0_20cm` |
| **Radiation** | Downward shortwave | Annual mean | `shortwave_radiation_J_m2` |
| | Downward longwave | Annual mean | `longwave_radiation_J_m2` |
| **Wind** | 10m U & V wind | Annual mean | `wind_10m_u_m_s`, `wind_10m_v_m_s` |
| | Wind speed magnitude | Annual mean (derived) | `wind_10m_speed_m_s` |
| **Evapotranspiration** | Total actual ET | Annual mean | `total_evaporation_m` |
| **Vegetation** | LAI (high veg) | Annual mean | `lai_high_veg` |
| | LAI (low veg) | Annual mean | `lai_low_veg` |

**Total: 15 variables (13 primary + 2 seasonality)**

## Data Processing Steps

### 1. Download ERA5-Land Monthly Data

**Source:** Copernicus Climate Data Store (CDS)
**Dataset:** `reanalysis-era5-land-monthly-means`
**Period:** 1999-2020 (22 years × 12 months = 264 files)
**Resolution:** 0.1° (~11 km)
**Size:** ~15 GB compressed

**Required:**
- CDS API key (already set in config)
- `ecmwfr` R package
- ~20 GB disk space in `era5_data/`

### 2. Calculate Annual Aggregates

For each variable, calculate from 12 monthly layers:

```
annual_mean = mean(Jan, Feb, ..., Dec)
seasonality = std_dev(Jan, Feb, ..., Dec)
```

**For accumulated variables (precipitation, snow):**
```
annual_total = sum(all monthly accumulated values)
```

### 3. Depth-Weight Soil Layers

ERA5-Land provides soil variables at 4 depths:
- Layer 1: 0-7 cm
- Layer 2: 7-28 cm
- Layer 3: 28-100 cm
- Layer 4: 100-255 cm

To standardize to **0-20 cm** (matching your OC standardization):

```
value_0_20cm = (layer1_value × 7/20) + (layer2_value × 13/20)
```

Where:
- Layer 1 contributes full 7 cm thickness
- Layer 2 contributes 13 cm (from 7 to 20 cm)
- Total target thickness: 20 cm

### 4. Extract at Point Locations

For each observation point (1.46M records):
- Extract annual mean values
- Add as new columns to point dataset

**Result:** Point dataset grows from 44 → 59 columns

### 5. Rasterize to 0.5° Global Grid

For each variable:
1. Create 0.5° template grid (720 × 360 cells)
2. Aggregate point observations to grid cells → `mean(points_in_cell)`
3. Save as compressed GeoTIFF

**Result:** 15 GeoTIFFs in `spatialized_layers/climate/`
**Size:** ~6 MB total (highly compressed)

## Implementation Status

### ✓ Code Structure Created
- `03_extract_climate.R` — Main module (skeleton)
- `03_extract_climate_helper.R` — Reference implementation

### ⏳ To Complete Implementation

You'll need to:

1. **Set up CDS API** (one time)
   ```r
   # Install if needed
   install.packages("ecmwfr")
   
   # Register key (run once)
   library(ecmwfr)
   wf_set_key(user = "cds", 
              key = Sys.getenv("CDS_API_KEY"))
   ```

2. **Download monthly data** (~4-6 hours runtime)
   ```r
   source("03_extract_climate_helper.R")
   # Uncomment download lines and run
   ```

3. **Complete the aggregation functions** in `03_extract_climate.R`
   - Load monthly NetCDF files
   - Calculate annualized statistics
   - Depth-weight soil layers
   - Extract at points
   - Rasterize

4. **Test and validate**
   - Check extracted values make sense (e.g., temperature in K, ~250-310)
   - Verify point coverage
   - Inspect raster outputs

## File Structure After Completion

```
Global_MRT_code/
├── outputs/
│   ├── 02_mineral_soils_filtered.rds
│   └── 03_with_climate.rds              ← Point data + climate
│
└── spatialized_layers/climate/
    ├── temperature_2m_mean_K_0.5deg.tif
    ├── temperature_seasonality_K_0.5deg.tif
    ├── precipitation_annual_m_0.5deg.tif
    ├── precipitation_seasonality_m_0.5deg.tif
    ├── solid_precipitation_m_0.5deg.tif
    ├── soil_temperature_0_20cm_K_0.5deg.tif
    ├── soil_moisture_0_20cm_0.5deg.tif
    ├── shortwave_radiation_J_m2_0.5deg.tif
    ├── longwave_radiation_J_m2_0.5deg.tif
    ├── wind_10m_u_m_s_0.5deg.tif
    ├── wind_10m_v_m_s_0.5deg.tif
    ├── wind_10m_speed_m_s_0.5deg.tif
    ├── total_evaporation_m_0.5deg.tif
    ├── lai_high_veg_0.5deg.tif
    └── lai_low_veg_0.5deg.tif
```

## Key Decisions Made

1. **Time period:** 1999-2020 (22 years of consistent, high-quality data)
2. **Aggregation:** All annual means (except seasonality = std dev)
3. **Depth standardization:** 0-20 cm (matches OC standardization)
4. **Spatial resolution:** 0.5° (matches CMIP6, manageable file sizes)
5. **Format:** Compressed float32 GeoTIFFs (efficient storage)

## Next Steps

Once climate extraction is complete:

1. **Step 4:** Productivity extraction (DMP from Copernicus)
2. **Step 5:** Depth standardization (OC to 0-20 cm)
3. **Step 6:** SoilGrids variables (clay, pH, CEC, nitrogen)
4. **Step 7:** Topography (slope, aspect, elevation)
5. **Step 8:** Soil classification (WRB)
6. **Step 9:** Land cover (ESA WorldCover)
7. **Step 10:** Bulk density gap-filling
8. **Step 11:** Random Forest modeling

## Estimated Runtimes

- Download ERA5 (1999-2020): 4-6 hours (network dependent)
- Aggregate to annual: 1-2 hours
- Extract at points: 30 min
- Rasterize: 5-10 min
- **Total: ~6-9 hours**

## Resources

- ERA5-Land documentation: https://confluence.ecmwf.int/display/CKB/ERA5-Land:+data+documentation
- CDS API guide: https://cds.climate.copernicus.eu/how-to-api
- ecmwfr R package: https://github.com/khufkens/ecmwfr
