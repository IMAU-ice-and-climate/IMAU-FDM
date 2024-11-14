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
nyears = 20 #spinup time or total runtime?

dir_ref = "/ec/res4/scratch/nld4814/FGRN055_era055/"
dir_mask = "/ec/res4/scratch/nld4814/FGRN055_era055/mask/"

dir3 = "/ec/res4/scratch/rusv/backup_HPC/output/ant27_exp_MO_2023/"
dir5 = "/scratch/rusv/data/output/NCL/ANT27/ERA5/"
pre_fname1 = "ECMWF_ANT27_1979-2023_CB"

## load pointlist
file_pointlist = dir_ref + "IN_ll_FGRN055_GrIS_GIC_implicit.txt"
pointlist = pd.read_csv(file_pointlist,header=None)

## load mask
fname_lsm = dir2 + "ANT27_Masks.nc"
model_lsm = xr.open_dataset(fname_lsm)
new_LSM = model_lsm['IceMask']
new_Lat = model_lsm['lat']
new_Lon = model_lsm['lon']