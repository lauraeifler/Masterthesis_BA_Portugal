## Import packages
import xarray as xr
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
import math
import pandas as pd
import sys
from fwi_functions import calculate_ffmc, calculate_dmc, calculate_dc, calculate_isi, calculate_bui, calculate_fwi

# Loop through years
year = sys.argv[1]

# Set path
path_data = '/work/lh88wasu-BA_Portugal/Data/FWI/Preprocessed/'
outpath_data = '/work/lh88wasu-BA_Portugal/Data/FWI/FWI/'

# lLad data
data = xr.open_dataset(path_data + f'ERA5_Land_allvars_{year}.nc')

### FFMC ###
# Initialize the previous FFMC value
ffmc_prev = np.full((data.dims['latitude'], data.dims['longitude']), 85.0)  # Initial FFMC array
ffmc_values = []

# Iterate over the time points
for i, time_point in enumerate(data['valid_time'].values):
    # Extract temperature, relative humidity, and precipitation for the current day
    temp_day = data['t2m_celsius'].isel(valid_time=i).values  
    rh_day = data['rh'].isel(valid_time=i).values            
    wind_day = data['ws'].isel(valid_time=i).values 
    precip_day = data['tp_cum'].isel(valid_time=i).values         
        
    # Initialize an empty array to store the FFMC values for each pixel (same shape as temp_day)
    ffmc_today = np.zeros_like(temp_day)  # Same shape as the input 
    
    # Loop through each pixel (i, j) to calculate FFMC
    for row in range(temp_day.shape[0]):  
        for col in range(temp_day.shape[1]):  
            temp_pixel = temp_day[row, col]
            rh_pixel = rh_day[row, col]
            wind_pixel = wind_day[row, col]
            precip_pixel = precip_day[row, col]
            ffmc_prev_pixel = ffmc_prev[row, col] 
            
            # Calculate FFMC for the current pixel and store it in the corresponding position in ffmc_today
            ffmc_today[row, col] = calculate_ffmc(temp_pixel, rh_pixel, precip_pixel, wind_pixel, ffmc_prev_pixel)
    
    # Append the current day's FFMC values for all pixels to the list
    ffmc_values.append(ffmc_today)
    
    # Update the ffmc_prev for the next iteration (next day)
    ffmc_prev = ffmc_today.copy()   # array for the next day

ffmc_array = xr.DataArray(
    ffmc_values,  
    coords=[data['valid_time'].values, data['latitude'].values, data['longitude'].values], 
    dims=["valid_time", "latitude", "longitude"]
)

# Assign the FFMC DataArray to the dataset
data["ffmc"] = ffmc_array

# Add attributes to the new FFMC variable
data["ffmc"].attrs = {
    "units": "dimensionless",
    "long_name": "Fine Fuel Moisture Code"
}


### DMC ####
# Initialize the previous DMC value
dmc_prev = np.full((data.dims['latitude'], data.dims['longitude']), 6.0)  # Initial DMC array
dmc_values = []

# Iterate over the time points
for i, time_point in enumerate(data['valid_time'].values):
    # Extract temperature, relative humidity, and precipitation for the current day
    temp_day = data['t2m_celsius'].isel(valid_time=i).values  
    rh_day = data['rh'].isel(valid_time=i).values            
    precip_day = data['tp_cum'].isel(valid_time=i).values         
    month_day = data['valid_time'].isel(valid_time=i).dt.month.values.item()  
    
    # Initialize an empty array to store the DMC values for each pixel (same shape as temp_day)
    dmc_today = np.zeros_like(temp_day)  # Same shape as the input (504, 324)

    # Loop through each pixel (i, j) to calculate DMC
    for row in range(temp_day.shape[0]):  # 504 rows
        for col in range(temp_day.shape[1]):  # 324 columns
            temp_pixel = temp_day[row, col]
            rh_pixel = rh_day[row, col]
            precip_pixel = precip_day[row, col]
            dmc_prev_pixel = dmc_prev[row, col] 
            
            # Calculate DMC for the current pixel and store it in the corresponding position in dmc_today
            dmc_today[row, col] = calculate_dmc(temp_pixel, rh_pixel, precip_pixel, dmc_prev_pixel, month_day)
    
    # Append the current day's DMC values for all pixels to the list
    dmc_values.append(dmc_today)
    
    # Update the dmc_prev for the next iteration (next day)
    dmc_prev = dmc_today.copy()   # array for the next day

dmc_array = xr.DataArray(
    dmc_values,
    coords=[data['valid_time'].values, data['latitude'].values, data['longitude'].values], 
    dims=["valid_time", "latitude", "longitude"]
)

# Assign the DMC DataArray to the dataset
data["dmc"] = dmc_array

# Add attributes to the new DMC variable
data["dmc"].attrs = {
    "units": "dimensionless",
    "long_name": "Duff Moisture Code"
}


######## DC ########
# Initialize the previous DC value
dc_prev = np.full((data.dims['latitude'], data.dims['longitude']), 15.0)  # Initial DC array
dc_values = []

# Iterate over the time points
for i, time_point in enumerate(data['valid_time'].values):
    # Extract temperature, relative humidity, and precipitation for the current day
    temp_day = data['t2m_celsius'].isel(valid_time=i).values  
    precip_day = data['tp_cum'].isel(valid_time=i).values         
    month_day = data['valid_time'].isel(valid_time=i).dt.month.values.item()  
    
    # Initialize an empty array to store the DC values for each pixel (same shape as temp_day)
    dc_today = np.zeros_like(temp_day)  # Same shape as the input (504, 324)

    # Loop through each pixel (i, j) to calculate DC
    for row in range(temp_day.shape[0]):  # 504 rows
        for col in range(temp_day.shape[1]):  # 324 columns
            temp_pixel = temp_day[row, col]
            precip_pixel = precip_day[row, col]
            dc_prev_pixel = dc_prev[row, col] 
            
            # Calculate DC for the current pixel and store it in the corresponding position in dc_today
            dc_today[row, col] = calculate_dc(temp_pixel, precip_pixel, dc_prev_pixel, month_day)
    
    # Append the current day's DC values for all pixels to the list
    dc_values.append(dc_today)
    
    # Update the dmc_prev for the next iteration (next day)
    dc_prev = dc_today.copy()   # array for the next day

dc_array = xr.DataArray(
    dc_values,  
    coords=[data['valid_time'].values, data['latitude'].values, data['longitude'].values], 
    dims=["valid_time", "latitude", "longitude"]
)

# Assign the DC DataArray to the dataset
data["dc"] = dc_array

# Add attributes to the new DC variable
data["dc"].attrs = {
    "units": "dimensionless",
    "long_name": "Drought Code"
}



#### ISI ######
ffmc = data.ffmc
wind = data.ws

# Apply IS calculation
isi = xr.apply_ufunc(
    calculate_isi,
    wind, ffmc,
    vectorize=True,  # Apply to each pixel
    dask="parallelized",  # Enable parallel computation if using dask
    output_dtypes=[float]
)

# Add IS to the dataset
data["isi"] = isi
data["isi"].attrs = {
    "units": "dimensionless",
    "long_name": "Initial Spread Index"
}



######### BUI #########
dmc = data.dmc
dc = data.dc

# Apply BUI calculation
bui = xr.apply_ufunc(
    calculate_bui,
    dmc, dc,
    vectorize=True,  # Apply to each pixel
    dask="parallelized",  # Enable parallel computation if using dask
    output_dtypes=[float]
)

# Add BUI to the dataset
data["bui"] = bui
data["bui"].attrs = {
    "units": "dimensionless",
    "long_name": "Buildup Index"
}



####### FWI #############
bui = data.bui
isi = data.isi

# Apply FWI calculation
fwi = xr.apply_ufunc(
    calculate_fwi,
    bui, isi,
    vectorize=True,  # Apply to each pixel
    dask="parallelized",  # Enable parallel computation if using dask
    output_dtypes=[float]
)

# Add FWI to the dataset
data["fwi"] = fwi
data["fwi"].attrs = {
    "units": "dimensionless",
    "long_name": "Fire Weather Index"
}



################################ RESAMPLING ####################################
# Define a function to expand the data for each variable
def upscale_to_higher_resolution(data_array, upscale_factor):
    """
    Upscales the resolution of a 2D data array by repeating each element.
    Parameters:
        data_array (np.ndarray): Input 2D array.
        upscale_factor (int): Factor to upscale the resolution.
    Returns:
        np.ndarray: Upscaled 2D array.
    """
    return np.repeat(np.repeat(data_array, upscale_factor, axis=1), upscale_factor, axis=2)

# Upscale factor: 9 (from 9km to 1km)
upscale_factor = 9

# Create new coordinates for the 1km resolution
lat_9km = data['latitude'].values
lon_9km = data['longitude'].values
lat_1km = np.linspace(lat_9km[0], lat_9km[-1], len(lat_9km) * upscale_factor)
lon_1km = np.linspace(lon_9km[0], lon_9km[-1], len(lon_9km) * upscale_factor)

# Upscale each variable in the dataset
resampled_vars = {}
for var in ['t2m_celsius', 'rh', 'ws', 'tp', 'tp_cum', 'ffmc', 'dmc', 'dc', 'isi', 'bui', 'fwi']:
    if var in data:
        original_data = data[var].values  # Extract the data as a NumPy arraya
        resampled_data = upscale_to_higher_resolution(original_data, upscale_factor)
        resampled_vars[var] = (["valid_time", "latitude", "longitude"], resampled_data)  # Preserve dimensions


# Create a new Dataset with the resampled data
resampled_ds = xr.Dataset(
    resampled_vars,
    coords={
        "valid_time": data["valid_time"].values,
        "latitude": lat_1km,
        "longitude": lon_1km
    }
)

resampled_ds.to_netcdf(outpath_data + f'ERA5_Land_FWI_1km_{year}.nc')

