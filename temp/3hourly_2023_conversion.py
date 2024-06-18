import numpy as np
import xarray as xr
import os 

#from Sanne

os.chdir('/ec/res4/scratch/nld4814/RACMO23p2_GR/raw/FGRN055/yearly')

varlist=["evap","ff10m","precip","sndiv","snowfall","snowmelt"]
file_append=".KNMI-2023.FGRN055.BN_RACMO2.3p2_ERA5_3h_1940_FGRN055.3H.nc"

for v in varlist:
	data = xr.open_dataset(v+"_ANT27_forFDM_Year2023.nc")

data_precip = xr.open_dataset('precip_ANT27_forFDM_Year2023.nc')
data_snowfall = xr.open_dataset('snowfall_ANT27_forFDM_Year2023.nc')
data_evap = xr.open_dataset('evap_ANT27_forFDM_Year2023.nc')
data_snowmelt = xr.open_dataset('snowmelt_ANT27_forFDM_Year2023.nc')

data_precip['precip'] = data_precip['precip']/(3600*3)
data_precip.to_netcdf('precip_ANT27_forFDM_Year2023b.nc')
data_snowfall['snowfall'] = data_snowfall['snowfall']/(3600*3)
data_snowfall.to_netcdf('snowfall_ANT27_forFDM_Year2023b.nc')

data_evap['evap'] = data_evap['evap']/(3600*3)
data_evap.to_netcdf('evap_ANT27_forFDM_Year2023b.nc')

data_snowmelt['snowmelt'] = data_snowmelt['snowmelt']/(3600*3)
data_snowmelt.to_netcdf('snowmelt_ANT27_forFDM_Year2023b.nc')
