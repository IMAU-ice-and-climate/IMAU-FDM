# import modules 
import os
import glob
import pandas as pd
import xarray as xr
import numpy as np
import math
import datetime as dt
from datetime import timedelta, datetime

# set working directory
os.chdir('/perm/nld3562/code/IMAU-FDM/post-process/make_1d_files/')
variab = "FirnAir"           # FirnAir (can be any output from 1D files)
dir1 = "/perm/nld3562/code/IMAU-FDM/reference/FGRN055/"                 # directory for input pointlist  
dir2 = "/perm/nld3562/code/IMAU-FDM/reference/FGRN055/"          # directory for mask 
dir3 = "/ec/res4/scratch/nld3562/daily/processed_output/"  # directory for 1D output files 
dir4 = "/ec/res4/scratch/nld3562/daily/processed_output/"                   # directory save gridded file 

pre_fname1 = "FGRN055_era055"	  		# file name 1D output files 
fname1 = dir1 + "IN_ll_FGRN055.txt"    # pointlist 
Grid1 = pd.read_csv(fname1,header=None)			# open pointlist 
nof = len(Grid1) 					# number of files 

fout = dir4 +"FDMv12A_" + variab + "_1939-2023_ERA5_FGRN055_mean.nc"   # include not_detrended, when FirnAir or zs are not detrended

# open mask 
fname_lsm = dir2 + "FGRN055_Masks.nc"                     # Antarctica mask
model_lsm = xr.open_dataset(fname_lsm)			# open mask
model_IceMask = model_lsm['LSM_GR']			# IceMask variable
model_Lat = model_lsm['lat']				# Latitude variable
model_Lon = model_lsm['lon']				# Longitude variable
n_lat = len(model_IceMask)				# Dimension length latitudes
n_lon = len(model_IceMask[0])				# Dimension length longitudes 
n_time = 1
TimeY = np.empty([n_time])

# Create matrix for gridded output 
Matrix_1 = np.empty([n_time,n_lat,n_lon])           # Empty 3D array  
Matrix_1[:,:,:] = np.nan			   # Fill with nan	

# Create gridded data 
for nn in range(0,nof):						  	# Loop over grid points  
        numb = nn + 1						  	# This value corresponds to the number of the 1D output file 
        print(numb)						  	# Keep track of the file being processed 
        fname = dir3 + pre_fname1 + "_1D_" + str(numb) + ".nc"    	# 1D output file 
       
        if os.path.exists(fname):				  	# If this path exists, then open the file 	
                ec_out = xr.open_dataset(fname)				
                ds_full = np.array(ec_out[variab][:])  # Load full timeseries

                xx = Grid1.iloc[nn,5].astype(int)		  	# Obtain x location		
                yy = Grid1.iloc[nn,6].astype(int)		  	# Obtain y location 	
                print(xx)
                print(yy)
                Matrix_1[:,xx,yy] = ds_full[:]	

#save gridded data
output_file = xr.Dataset(data_vars=dict(FAC=(["time","rlat","rlon"], Matrix_1.data), lat=(["rlat","rlon"],model_Lat.data), lon=(["rlat","rlon"],model_Lon.data)), coords=dict(time=TimeY)) #, dims=["time","rlat","rlon"])
output_file.to_netcdf(fout)