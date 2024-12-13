#/bin/sh

rsync -avz --partial --append-verify /Volumes/imau01/rapid/RACMO2.3p2/FGRN055_1940/3-hourly/snowmelt.KNMI-2021*.nc hpc-login:/ec/res4/scratch/nld4814/RACMO23p2_GR/raw/FGRN055/decadal/
rsync -avz --partial --append-verify /Volumes/imau01/rapid/RACMO2.3p2/FGRN055_1940/3-hourly/snowfall.KNMI-2021*.nc hpc-login:/ec/res4/scratch/nld4814/RACMO23p2_GR/raw/FGRN055/decadal/
rsync -avz --partial --append-verify /Volumes/imau01/rapid/RACMO2.3p2/FGRN055_1940/3-hourly/tskin.KNMI-2021*.nc hpc-login:/ec/res4/scratch/nld4814/RACMO23p2_GR/raw/FGRN055/decadal/
rsync -avz --partial --append-verify /Volumes/imau01/rapid/RACMO2.3p2/FGRN055_1940/3-hourly/evap.KNMI-2021*.nc hpc-login:/ec/res4/scratch/nld4814/RACMO23p2_GR/raw/FGRN055/decadal/
rsync -avz --partial --append-verify /Volumes/imau01/rapid/RACMO2.3p2/FGRN055_1940/3-hourly/precip.KNMI-2021*.nc hpc-login:/ec/res4/scratch/nld4814/RACMO23p2_GR/raw/FGRN055/decadal/
rsync -avz --partial --append-verify /Volumes/imau01/rapid/RACMO2.3p2/FGRN055_1940/3-hourly/ff10m.KNMI-2021*.nc hpc-login:/ec/res4/scratch/nld4814/RACMO23p2_GR/raw/FGRN055/decadal/
rsync -avz --partial --append-verify /Volumes/imau01/rapid/RACMO2.3p2/FGRN055_1940/3-hourly/sndiv.KNMI-2021*.nc hpc-login:/ec/res4/scratch/nld4814/RACMO23p2_GR/raw/FGRN055/decadal/

rsync -avz --partial --append-verify /Volumes/imau01/rapid/RACMO2.3p2/FGRN055_1940/3-hourly/NC_FDM-3h/snowmelt.KNMI-2023*.nc hpc-login:/ec/res4/scratch/nld4814/RACMO23p2_GR/raw/FGRN055/yearly/
rsync -avz --partial --append-verify /Volumes/imau01/rapid/RACMO2.3p2/FGRN055_1940/3-hourly/NC_FDM-3h/snowfall.KNMI-2023*.nc --append hpc-login:/ec/res4/scratch/nld4814/RACMO23p2_GR/raw/FGRN055/yearly/
rsync -avz --partial --append-verify N/Volumes/imau01/rapid/RACMO2.3p2/FGRN055_1940/3-hourly/NC_FDM-3h/tskin.KNMI-2023*.nc hpc-login:/ec/res4/scratch/nld4814/RACMO23p2_GR/raw/FGRN055/yearly/
rsync -avz --partial --append-verify /Volumes/imau01/rapid/RACMO2.3p2/FGRN055_1940/3-hourly/NC_FDM-3h/evap.KNMI-2023*.nc hpc-login:/ec/res4/scratch/nld4814/RACMO23p2_GR/raw/FGRN055/yearly/
rsync -avz --partial --append-verify /Volumes/imau01/rapid/RACMO2.3p2/FGRN055_1940/3-hourly/NC_FDM-3h/precip.KNMI-2023*.nc hpc-login:/ec/res4/scratch/nld4814/RACMO23p2_GR/raw/FGRN055/yearly/
rsync -avz --partial --append-verify /Volumes/imau01/rapid/RACMO2.3p2/FGRN055_1940/3-hourly/NC_FDM-3h/ff10m.KNMI-2023*.nc hpc-login:/ec/res4/scratch/nld4814/RACMO23p2_GR/raw/FGRN055/yearly/
rsync -avz --partial --append-verify /Volumes/imau01/rapid/RACMO2.3p2/FGRN055_1940/3-hourly/NC_FDM-3h/sndiv.KNMI-2023*.nc hpc-login:/ec/res4/scratch/nld4814/RACMO23p2_GR/raw/FGRN055/yearly/
