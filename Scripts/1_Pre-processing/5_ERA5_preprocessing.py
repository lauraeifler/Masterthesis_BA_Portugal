## Import packages
import xarray as xr
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
import math
import pandas as pd
import sys

# Loop through years
year = int(sys.argv[1])

# Set path
path_precip = 'path-to-precipitation-data'
path_clim = 'path-to-other-ERA5-data' #all others in one nc file per year
outpath_data = '/work/lh88wasu-BA_Portugal/Data/FWI/Preprocessed/'

# Load data of other variables (noon values)
data = xr.open_dataset(path_clim + f'ERA5_Land_climvars_{year}.nc')

################################## CUMULATIVE PRECIPITATION NOON TO NOON #######################################
# import dataset
data_precip = xr.open_dataset(path_precip + f'ERA5_Land_precip_{year}.nc')

precip = data_precip["tp"] * 1000 # from m to mm

cumulative_precip = precip.resample(valid_time='24h', offset='12h').sum(skipna=True)

# Adjust time to represent the midpoint (12:00) of each period
cumulative_precip['valid_time'] = cumulative_precip['valid_time'] + pd.Timedelta(hours=24)

cumulative_precip = cumulative_precip.isel(valid_time=slice(None, -1))

# all ocean pixels to nan
ocean_mask = np.isnan(data['t2m'])
cumulative_precip = cumulative_precip.where(~ocean_mask, other=np.nan)

########################## Pre-process other variables ###########################

# add cumulative precipitation values to data of all other variables
data['tp_cum'] = cumulative_precip

## Temperature is needed as celsius
temp_kelvin = data.t2m
temp_celsius = temp_kelvin - 273.15
data['t2m_celsius'] = temp_celsius

dew_kelvin = data.d2m
dew_celsius = dew_kelvin - 273.15
data['d2m_celsius'] = dew_celsius

## Calculate the total wind based on the eastward (u) and northward (v) wind component
u_wind = data["u10"]  
v_wind = data["v10"] 

# Calculate total wind speed
total_wind = np.sqrt(u_wind**2 + v_wind**2)

# Add the new variable to the dataset
data["ws"] = (total_wind)  # Adjust dimensions as needed
data["ws"] = data["ws"] * 3.6
data["ws"].attrs = {"units": "km/h", "long_name": "Total Wind Speed"}


## Calculate the relative humidity based on the temperature and dew point temperature
t2m = data['t2m_celsius']
d2m = data['d2m_celsius']

# Constants for the saturation vapor pressure equation
a = 17.62
b = 243.12
c = 6.112

# Calculate saturation vapor pressure
e_temp = c * np.exp((a * t2m) / (b + t2m))
e_dew = c * np.exp((a * d2m) / (b + d2m))

# Calculate relative humidity
relative_humidity = 100 * (e_dew / e_temp)

data['rh'] = relative_humidity
data['rh'].attrs = {"units": "%", "long_name": "Relative Humidity"}

# save to netCDF
data.to_netcdf(outpath_data + f'ERA5_Land_allvars_{year}.nc')








