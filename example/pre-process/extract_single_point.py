# python3

"""
This code extracts the data for a single point from the FDM timeseries data files and saves it to the example directory.
The number of dimensions, etc is maintained.
Outputs to example/input/timeseries; can then be run from `rundir/launch_example_job.sc` but if you change the point, domain, or forcing, you'll need to update this in `source/model_settings.f90`

The values for rlat_i, rlon_i, and r_lon_ts are extracted from logfiles from a run done on the ecmwf. 
Their values are 1 fewer than those reported in the logfiles because of 0 vs 1 indexing.
Fnumb corresponds to the file part number.

Elizabeth Case
2025-06-27
"""

import xarray as xr
import numpy as np

point = "1"
fnumb = "24"

ts_data_dir = "/ec/res4/scratch/nld4814/FGRN055_era055/input/timeseries/"
ts_save_dir = "/home/nld4814/perm/code/IMAU-FDM/example/input/timeseries/"
ts_data_suffix = f"_FGRN055_era055_1939-2023_p{fnumb}.nc"
ts_save_suffix = f"_FGRN055_era055_1939-2023_point_{point}.nc"

#all these values are 1 less than those reported in the logfiles because of 0 vs 1 indexing
rlat_i = 54
#rlon_i = 140
rlon_i_ts = 2

def get_single_point_ts(var,rlat_i,rlon_i,fnumb,point,data_dir=ts_data_dir,save_dir=ts_save_dir,save=False,do_return=False):
    """
    Get the value of a variable at a single point (rlat_i, rlon_i)
    from the timeseries data files.
    """
    ts_file = f"{data_dir}{var}{ts_data_suffix}"
    ds = xr.open_dataset(ts_file)

    # if you give single points as a list, xarray maintains all the dimensions
    ds_i = ds.isel(rlat=[rlat_i], rlon=[rlon_i_ts], drop=False)
    
    if save:
        ds_i.to_netcdf(save_dir+f"{var}{ts_save_suffix}", mode='w')

    if do_return:
        return ds_i
    
vars = ["evap","ff10m","precip","sndiv","snowfall","snowmelt","tskin"]

for var in vars:

    print(f"Processing {var} for timeseries data...")
    get_single_point_ts(var, rlat_i, rlon_i_ts, fnumb, point, save=True, do_return=False)


# run this to check whether remove value for Ptot is correct

# ds_i = ds.isel(rlat=[rlat_i], rlon=[rlon_i_ts], drop=False)

# dtobs = 10800 # 3 hourly observations in seconds
# ts_minimum =  1.e-04 

# remove_Val = 0.
# PreVal = ds_i[var].values * dtobs


# for ind_t in range(len(ds_i.time)-1):
        
#     if PreVal[ind_t] < ts_minimum:
#         remove_Val = remove_Val + PreVal[ind_t]
#         PreVal[ind_t] = 0.

# print("lat: ", ds_i.lat.values[0][0])
# print("lon: ", ds_i.lon.values[0][0])
# print("remove_Psol: ", remove_Val[0][0])