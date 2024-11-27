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
os.chdir('/perm/rusv/IMAU_FDM/scripts/scripts/')

# settings 
variab = "dens"         # dens, temp (can be any output from 2D files)
variab2 = "depth" 	# depth (or age)
nyears = 45 		# amount of years 
dens_t = 830 		# critical density threshold 
variab_output = "z830"  # z830: depth of critical density level z830

# directories and files  
dir1 = "/perm/rusv/IMAU_FDM/code_CESM2/future/DATA/"                 # directory for input pointlist  
dir2 = "/scratch/rusv/data/input/ant27_era_files/averages/"          # directory for mask 
dir3 = "/ec/res4/scratch/rusv/backup_HPC/output/ant27_exp_MO_2023/"  # directory for 1D output files 
dir4 = "/scratch/rusv/data/output/NCL/ANT27/ERA5/"                   # directory save gridded file  

pre_fname1 = "ECMWF_ANT27_1979-2023_CB"	  		# file name 1D output files 
fname1 = dir1 + "IN_ll_ANT27_CB1979-2020_nor_MO.txt"    # pointlist 
Grid1 = pd.read_csv(fname1,header=None)			# open pointlist 
nof = len(Grid1)					# number of files 
print(nof)

# output file name  
fout = dir4 +"FDMv12A_" + variab_output + "_1979-2023_ERA5_ANT27.nc"

# Open mask 
fname_lsm = dir2 + "ANT27_Masks.nc"                     # Antarctica mask
model_lsm = xr.open_dataset(fname_lsm)			# open mask
model_IceMask = model_lsm['IceMask']			# IceMask variable
model_Lat = model_lsm['lat']				# Latitude variable
model_Lon = model_lsm['lon']				# Longitude variable
n_lat = len(model_IceMask)				# Dimension length latitudes
n_lon = len(model_IceMask[0])				# Dimension length longitudes 

# Create timestep dimension
ntime = np.int32(nyears*365.25/30+1)			# Define amount of timesteps
TimeY = np.empty([ntime])		                # Empty time array 		

for t in range(0, ntime):
        TimeY[t] = 1979.08214 + (t/36.5250)*3           # Fill with fractional years (including leap years)

# create matrix for gridded output 
Matrix_1 = np.empty([ntime,n_lat,n_lon])           # Empty 3D array  
Matrix_1[:,:,:] = np.nan			   # Fill with nan	

# Create gridded data 
for nn in range(0,nof):					  		# Loop over grid points  
        numb = nn + 1						  	# This value corresponds to the number of the 1D output file 
        print(numb)						  	# Keep track of the file being processed 
        fname = dir3 + pre_fname1 + "_2D_" + str(numb) + ".nc"          # 1D output file 
       
        if os.path.exists(fname):				  	# If this path exists, then open the file 	
                ec_out = xr.open_dataset(fname)				
                ds = np.array(ec_out[variab][:,:ntime])			# Obtain variable
                dz = np.array(ec_out[variab2][:,:ntime])		# Obtain variabl2 (depth/age)
                dz = np.where(ds[:,:]>dens_t,np.nan,dz)			# Fill with np.nan if density threshold is exceeded 
                dz = np.nanmax(dz,axis=0)				# Maximum depth with density threshold
                xx = Grid1.iloc[nn,5].astype(int)		  	# Obtain x location		
                yy = Grid1.iloc[nn,6].astype(int)		  	# Obtain y location 	
                print(xx)
                print(yy)
                Matrix_1[:,xx,yy] = dz[:]				
  
print(variab)
# create gridded data matrix
Matrix_2 = xr.Dataset(data_vars=dict(variab=(["time","rlat","rlon"], Matrix_1.data), lat=(["rlat","rlon"],model_Lat.data), lon=(["rlat","rlon"],model_Lon.data)), coords=dict(time=TimeY)) #, dims=["time","rlat","rlon"])

Matrix_2['%s'%variab_output] = Matrix_2['variab']
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
Matrix_2['rotated_pole'] = model_lsm['rotated_pole']
Matrix_2 = Matrix_2.sel(time=slice("1979-01-01", "2024-01-01"))

if variab_output == 'z830':
	Matrix_2['z830'] = Matrix_2.z830.assign_attrs(units = 'm')
	Matrix_2['z830'] = Matrix_2.z830.assign_attrs(long_name = 'Depth of critical density level 830 kg/m3')
	Matrix_2['z830'] = Matrix_2.z830.assign_attrs(standard_name = 'z830')
	Matrix_2 = Matrix_2.assign_attrs(title = 'Depth of z830')

Matrix_2.to_netcdf(fout)  # save output file


