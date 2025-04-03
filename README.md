# Masterthesis_BA_Portugal
This repository contains scripts and data for a Master's Thesis analyzing burned areas (BAs) in protected vs. unprotected regions in Portugal.

## **Overview**
The study analyzes the influence of non-climate and climate predictors on burned area in mainland Portugal, while focusing of protected and non-protected areas. Three models are applied to analyze the influence of predictors on BA - Quantile Regression, Generalized Additive Model, and Random Forest.
The pre-processing and visualizations are done in Python , the modeling is performed in R version 4.2.


## **Data**
### **Included in this Repository**
- **Portugal shapefile**: Used for clipping raster data.
- **WDPA shapefile and raster (1 km resolution)**: Represents protected areas.

### **Required External Data (Not Included)**
You need to obtain the following datasets manually:
- **Burned Area Data** provided by European Space Agency, find [here](https://catalogue.ceda.ac.uk/uuid/3628cb2fdba443588155e15dee8e5352)
- **Land Cover Data** provided by Copernicus Climate Data Store, find [here](https://cds.climate.copernicus.eu/datasets/satellite-land-cover?tab=download)
- **ERA5-Land Hourly Data** provided by Copernicus Climate Data Store, find [here](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land?tab=download)


## **Scripts**
### **Pre-processing**
- Scripts handle clipping and resampling of each dataset for the analysis.
Burned Area:
- Two scripts to first clip BA files to the mainland Portugal shape, then resample to 1 km with percentage of BA as new pixel value
WDPA:
- Script to rasterize protected areas andd save as netCDF
Land cover:
- Creates new (broader) classification system for land cover (LC) classes, resamples to 1 km by keeping the most abundant LC value as new pixel value
ERA5-Land:
- Pre-processing of climate variables (cumulative precipitation, calculate relative humidity and total wind speed)
- Calculation of FWI components in `6_FWI_calculation.py`, requires `fwi_function.py` for computation
- Script also includes resampling from 9 km to 1 km
Data table:
- Script `7_Create_table.ipynb` for sampling and generation of table used for modeling analysis

### **Modeling**
- R script for statistical modeling:
  - Quantile Regression
  - Generalized Additive Model
  - Random Forest
- Includes code for generating visualizations used in the thesis
- Computation in R

### **Visualization**
- Scripts for producing figures related to FWI and Land Cover.

## **Requirements**
Ensure you have the following installed:
- Python
- R (version 4.2)
- Packages are listed in each script




