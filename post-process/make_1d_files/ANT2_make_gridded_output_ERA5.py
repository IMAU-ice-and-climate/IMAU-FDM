#Things to change
#set directory
#location of directories
#reading input (pre_fname1 and fname1)
#IMAU mask
#create arrays

#import modules
import os
import glob
import pandas as pd
import xarray as xr
import numpy as np
import math
import datetime as dt
from datetime import timedelta, datetime


#set directory
os.chdir('/perm/nld3562/code/IMAU-FDM/post-process/make_1d_files/')

#settings
variab = "FirnAir"     # Firn air content (m) 
dt = 0             # (0=no detrending, 1=detrending) use for contemporary simulation (1979-2020) (FAC and zs)
atzero = 1  	   # (0=start at zero) use for surface elevation change (zs)

#location of directories 
dir1 = "/perm/nld3562/code/IMAU-FDM/reference/ANT11_larsenc/"                 # directory for input pointlist  
dir2 = "/perm/nld3562/code/IMAU-FDM/reference/ANT11_larsenc/"          # directory for mask 
dir3 = "/ec/res4/scratch/nld3562/RACMO_era05_larsenc_3hourly_3K/output/"  # directory for 1D output files 
dir4 = "/ec/res4/scratch/nld3562/RACMO_era05_larsenc_3hourly_3K/output/"                   # directory save gridded file  

#reading input 
pre_fname1 = "ANT11_larsenc_ANT11_era05_larsenc" # name of raw data
fname1 = dir1 + "IN_ll_ANT11_larsenc.txt" 		 # grid list 
Grid1 = pd.read_csv(fname1,header=None) 			 # open grid list 
numpoints1 = Grid1.shape[0]
num_in_grid = Grid1.shape[1]
print(numpoints1)
print(num_in_grid)

#IMAU-FDM mask 
fname_lsm = dir2 + "PXANT11_masks_larsenc.nc"    # IMAU-FDM mask 
model_lsm = xr.open_dataset(fname_lsm)                               # open IMAU-FDM mask 
LSM = model_lsm['LSM']
Lat = model_lsm['lat']
Lon = model_lsm['lon']
nlat = len(LSM)
nlon = len(LSM[0])
print(LSM)
print(nlat)
print(nlon)

#open example file 
fname = dir3 + pre_fname1 + "_1D_607.nc"
forT = xr.open_dataset(fname)
Time_zs = forT['FirnAir'][:]   # time variable #changed from 1524 to 1719
ntime = len(Time_zs)           # length of time array 

#create arrays 
Matrix_27 = np.empty([ntime,nlat,nlon])      	# create output array 
Matrix_27[:,:,:] = np.nan			# fill with nan values 
TimeY = np.empty([ntime])    			# create time array 

for t in range(0, ntime):
        TimeY[t] = 1979.02738 + (t/36.5)    # assign time, 10 day resolution, without leap years, in fractional years 

#process the raw data  
for nnn in range(0,(numpoints1)):        # loop over grid list 
        numb = nnn + 1
        print(numb)
        fname = dir3 + pre_fname1 + "_1D_" + str(numb) + ".nc"   

        if os.path.exists(fname):
                ec_out = xr.open_dataset(fname)			# open raw data file 
                D1D = np.array(ec_out[variab][:])           # retrieve data 
                xx = Grid1.iloc[nnn,5].astype(int)              # obtain x coordinate 
                yy = Grid1.iloc[nnn,6].astype(int)		# obtain y coordinate 
                print(xx)
                print(yy)
                if atzero == 0:
                    D1D[:] = D1D[:] - D1D[0]      # start at zero 
                if dt == 0:
                    Matrix_27[:,yy,xx] = D1D[:]   # write to output array  #changed xx,yy to yy,xx
                if dt == 1:
                    rc0 = (D1D[ntime]-D1D[0])/ntime        	 # detrend timeseries 
                    for tt in range(0,ntime):
                        Matrix_27[tt,yy,xx] = D1D[tt] - rc0*tt  # write to output array #changed xx,yy to yy,xx
        else:
                print("-------> this file is not available:", numb)

fout = dir4 + "FDMv1.2AD_RACMO2.3p2_" + variab + "_1979-2021_ERA5_ANT2.nc" # output file 

#save gridded data
output_file = xr.Dataset(data_vars=dict(FAC=(["time","rlat","rlon"], Matrix_27.data), lat=(["rlat","rlon"],Lat.data), lon=(["rlat","rlon"],Lon.data)), coords=dict(time=TimeY)) 
output_file.to_netcdf(fout)

