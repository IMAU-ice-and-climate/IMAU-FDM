import os
import glob
import pandas as pd
import xarray as xr
import numpy as np
import math
import datetime as dt
from datetime import timedelta, datetime

os.chdir('/perm/rusv/IMAU_FDM/scripts/scripts/')

variab = "zs"   #FirnAir 

dt = 0  #1 (0=no detrending, 1=detrending)
atzero = 1  #0 (0=start at zero,1=do not start at zero)  
nyears = 44 

dir1 = "/perm/rusv/IMAU_FDM/code_CESM2/future/DATA/"
dir2 = "/scratch/rusv/data/input/ant27_era_files/averages/"
dir3 = "/ec/res4/scratch/rusv/backup_HPC/output/ant27_exp_MO_2023/"
dir5 = "/scratch/rusv/data/output/NCL/ANT27/ERA5/"

pre_fname1 = "ECMWF_ANT27_1979-2023_CB"
fname1 = dir1 + "IN_ll_ANT27_CB1979-2020_nor_MO.txt"

Grid1 = pd.read_csv(fname1,header=None)

fname_lsm = dir2 + "ANT27_Masks.nc"
model_lsm = xr.open_dataset(fname_lsm)
new_LSM = model_lsm['IceMask']
new_Lat = model_lsm['lat']
new_Lon = model_lsm['lon']

new_Latr = np.linspace(0,239,240)
new_Lonr = np.linspace(0,261,262)

new_n_lat = len(new_LSM)
new_n_lon = len(new_LSM[0])

print(new_LSM)
print(new_n_lat)
print(new_n_lon)

fname_example = dir3 + pre_fname1 + "_1D_1.nc"
forT = xr.open_dataset(fname_example)
nsteps = np.int32(nyears*365.25/10+1)
Time_zs = forT['zs'][0:nsteps]  
 
print(Time_zs)
ntime = len(Time_zs)

tdata = forT.time*10
taxis = tdata
taxis_day = (taxis%365+1).astype(str)
taxis_year = ((taxis//365)+1979).astype(str)

print(ntime+0)

Matrix_1 = np.empty([ntime,new_n_lat,new_n_lon])
Matrix_1[:,:,:] = np.nan

TimeY = np.empty([ntime])

for t in range(0, ntime):
        TimeY[t] = 1979.02738 + (t/36.5250)       
print(TimeY[ntime-1]+0)

for nnn in range(0,18136):
        numb = nnn + 1
        print(numb)
        fname = dir3 + pre_fname1 + "_1D_" + str(numb) + ".nc"
       
        if os.path.exists(fname):
                ec_out = xr.open_dataset(fname)
                ds = np.array(ec_out[variab][:nsteps])
                xx = Grid1.iloc[nnn,5].astype(int)
                yy = Grid1.iloc[nnn,6].astype(int)
                print(xx)
                print(yy)
                if atzero == 0:
                    ds[:] = ds[:] - ds[0]
                if dt == 1:
                    rc0 = (ds[1533]-ds[0])/(1533)   #1533 is amount of timesteps in ANT spin up, check for Greenland
                    for sss in range(0,1533):
                        Matrix_1[sss,xx,yy] = ds[sss]-rc0*sss
                    jump = (Matrix_1[1532,xx,yy]-ds[1532]) 
                    for sss in range(1533,nsteps): 
                        Matrix_1[sss,xx,yy] = Matrix_1[sss,xx,yy] + jump 			
                else:
                    Matrix_1[:,xx,yy] = ds[:]
 
fout = dir5 +"FDMv12A_" + variab + "1979-2023_ERA5_ANT27_not_detrended_check.nc"

Matrix_2 = xr.Dataset(data_vars=dict(zs=(["time","rlat","rlon"], Matrix_1.data), lat=(["rlat","rlon"],new_Lat.data), lon=(["rlat","rlon"],new_Lon.data)), coords=dict(time=TimeY)) #, dims=["time","rlat","rlon"])

print(Matrix_2)
print(Matrix_2.rlat)
print(Matrix_2.rlon)

#assign attributes etc.
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

Matrix_2.to_netcdf(fout)


