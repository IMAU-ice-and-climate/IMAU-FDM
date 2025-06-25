"""
This script makes gridded output (3D: lon,lat,time) of 1D  files (time) from IMAU-FDMv12A 
Can be used for all variables from the 1D output files 

Notes: 
- Only FirnAir and zs should be detrended over the spin up period 
- If using with RACMO2.3p2-CESM >> Remove leap years in L56 and L60 by removing .25
- check attributes at the end of the script

Made by: Sanne Veldhuijsen (nov-2024) 
"""

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

# settings 
variab = "FirnAir"           # FirnAir (can be any output from 1D files)
dt = 1  		# 0=no detrending, 1=detrending --> Should be done for "FirnAir" and "zs" (include in file name if not)
atzero = 1  		# 0=start at zero,1=do not start at zero --> (Actually I realized I do not use this anymore, because first output is after 10 days...) But could be used for "zs".)
nyears = 84 		# amount of years 
ts_su = 10950 		# timestep spin up: 1533 for Antarctica 1979-2020 - check Greenland 

# directories and files  
dir1 = "/perm/nld3562/code/IMAU-FDM/reference/FGRN055/"                 # directory for input pointlist  
dir2 = "/perm/nld3562/code/IMAU-FDM/reference/FGRN055/"          # directory for mask 
dir3 = "/ec/res4/scratch/nld3562/daily/output/"  # directory for 1D output files 
dir4 = "/ec/res4/scratch/nld3562/daily/output/"                   # directory save gridded file  

pre_fname1 = "FGRN055_era055"	  		# file name 1D output files 
fname1 = dir1 + "IN_ll_FGRN055.txt"    # pointlist 
Grid1 = pd.read_csv(fname1,header=None)			# open pointlist 
nof = len(Grid1) 					# number of files 

# output file name  
fout = dir4 +"FDMv12A_" + variab + "_1939-2023_ERA5_FGRN055_detrended2.nc"   # include not_detrended, when FirnAir or zs are not detrended

# open mask 
fname_lsm = dir2 + "FGRN055_Masks.nc"                     # Antarctica mask
model_lsm = xr.open_dataset(fname_lsm)			# open mask
model_IceMask = model_lsm['LSM_GR']			# IceMask variable
model_Lat = model_lsm['lat']				# Latitude variable
model_Lon = model_lsm['lon']				# Longitude variable
n_lat = len(model_IceMask)				# Dimension length latitudes
n_lon = len(model_IceMask[0])				# Dimension length longitudes 

full_nt = np.int32(nyears * 365.25 + 1)   # Full number of time steps (daily)
step = 10                                 # Take every 10th timestep
ntime = full_nt // step                   # Reduced number of time steps
TimeY = np.empty([ntime])

for t in range(0, ntime):
    #TimeY[t] = 1939.02738 + ((t * step) / 36.5250) 
    TimeY[t] = 1939.02738 + ((t * step)) 

# Create matrix for gridded output 
Matrix_1 = np.empty([ntime,n_lat,n_lon])           # Empty 3D array  
Matrix_1[:,:,:] = np.nan			   # Fill with nan	

# Create gridded data 
for nn in range(0,nof):						  	# Loop over grid points  
        numb = nn + 1						  	# This value corresponds to the number of the 1D output file 
        print(numb)						  	# Keep track of the file being processed 
        fname = dir3 + pre_fname1 + "_1D_" + str(numb) + ".nc"    	# 1D output file 
       
        if os.path.exists(fname):				  	# If this path exists, then open the file 	
                ec_out = xr.open_dataset(fname)				
                ds_full = np.array(ec_out[variab][:full_nt])  # Load full timeseries
                ds = ds_full[::step]          
                xx = Grid1.iloc[nn,5].astype(int)		  	# Obtain x location		
                yy = Grid1.iloc[nn,6].astype(int)		  	# Obtain y location 	
                print(xx)
                print(yy)
                if atzero == 0:					  	# Start at zero
                    ds[:] = ds[:] - ds[0]
                if dt == 1:			                  	# Detrending			
                    ts_su_ds = ts_su // step
                    rc0 = (ds[ts_su_ds]-ds[0]) / ts_su
                    for ss in range(0, ts_su_ds):
                        Matrix_1[ss, xx, yy] = ds[ss] - rc0 * (ss * step)       	# Correct for this trend during spin up 
                    jump = (Matrix_1[ts_su_ds-1,xx,yy]-ds[ts_su_ds-1])
                    for ss in range(ts_su_ds, ntime): 
                        Matrix_1[ss,xx,yy] = Matrix_1[ss,xx,yy] + jump
                else:
                    Matrix_1[:,xx,yy] = ds[:]				# no detrending
 
# create gridded data matrix
Matrix_2 = xr.Dataset(data_vars=dict(variab=(["time","rlat","rlon"], Matrix_1.data), lat=(["rlat","rlon"],model_Lat.data), lon=(["rlat","rlon"],model_Lon.data)), coords=dict(time=TimeY)) #, dims=["time","rlat","rlon"])

Matrix_2['%s'%variab] = Matrix_2['variab']
Matrix_2 = Matrix_2.drop_vars(['variab'])

print(Matrix_2)
print(Matrix_2.rlat)
print(Matrix_2.rlon)

#assign attributes, dimensions etc.
Matrix_2 = Matrix_2.assign_attrs(Conventions = 'CF-1.4')
Matrix_2 = Matrix_2.assign_attrs(source = 'IMAU-FDM v1.2A')
Matrix_2 = Matrix_2.assign_attrs(Domain = 'ANT27')
Matrix_2 = Matrix_2.assign_attrs(Experiment = 'ERA5-3H_RACMO2.3p2')
Matrix_2 = Matrix_2.assign_attrs(Institution = 'IMAU (Sanne Veldhuijsen)')
Matrix_2 = Matrix_2.assign_attrs(CreationDate = '13/06/2024')
Matrix_2 = Matrix_2.assign_coords(rlat = model_lsm.rlat)
Matrix_2 = Matrix_2.assign_coords(rlon = model_lsm.rlon)
Matrix_2['lat'] = Matrix_2.lat.assign_attrs(longname = 'latitude')
Matrix_2['lat'] = Matrix_2.lat.assign_attrs(units = 'degrees_north')
Matrix_2['lon'] = Matrix_2.lon.assign_attrs(longname = 'longitude')
Matrix_2['lon'] = Matrix_2.lon.assign_attrs(units = 'degrees_east')
#Matrix_2['rotated_pole'] = model_lsm['rotated_pole']
#Matrix_2 = Matrix_2.sel(time=slice("1979-01-01", "2024-01-01"))

if variab == 'zs':
	Matrix_2['zs'] = Matrix_2.zs.assign_attrs(units = 'm')
	Matrix_2['zs'] = Matrix_2.zs.assign_attrs(long_name = 'Surface elevation')
	Matrix_2['zs'] = Matrix_2.zs.assign_attrs(standard_name = 'zs')
	Matrix_2 = Matrix_2.assign_attrs(title = 'Surface elevation')

if variab == 'FirnAir':
	Matrix_2['FirnAir'] = Matrix_2.FirnAir.assign_attrs(units = 'm')
	Matrix_2['FirnAir'] = Matrix_2.FirnAir.assign_attrs(long_name = 'Firn air content')
	Matrix_2['FirnAir'] = Matrix_2.FirnAir.assign_attrs(standard_name = 'FAC')
	Matrix_2 = Matrix_2.assign_attrs(title = 'Firn air content')

if variab == 'vfc':
	Matrix_2['vfc'] = Matrix_2.vfc.assign_attrs(units = 'm')
	Matrix_2['vfc'] = Matrix_2.vfc.assign_attrs(long_name = 'Vertical velocity due to firn compaction')
	Matrix_2['vfc'] = Matrix_2.vfc.assign_attrs(standard_name = 'vfc')
	Matrix_2 = Matrix_2.assign_attrs(title = 'Vertical velocity due to firn compaction')

	
Matrix_2.to_netcdf(fout)  # save output file


