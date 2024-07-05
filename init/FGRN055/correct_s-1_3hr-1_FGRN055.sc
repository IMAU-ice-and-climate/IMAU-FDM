#!/bin/ksh -x

raw_dir="$SCRATCH/FM_Data/INPUT/FGRN055_rawfiles/"
ave_dir="$SCRATCH/FM_Data/INPUT/FGRN055_averages/"
files_dir="$SCRATCH/FM_Data/INPUT/FGRN055_files/"

nco_dir="/usr/local/apps/nco/4.3.7/bin/"
ncks=${nco_dir}"ncks"
ncap2=${nco_dir}"ncap2"
ncrename=${nco_dir}"ncrename"

varlist=$1 #"precip snowfall evap snowmelt sndiv"

for varname in $varlist; do

  temp1=${raw_dir}${varname}"_stemp1.nc"
  temp2=${raw_dir}${varname}"_stemp2.nc"
  temp3=${raw_dir}${varname}"_stemp3.nc"

  (( part = 1 ))
  while [ $part -lt 74 ]; do
    if [ $part -eq 24]; then
      (( part = part + 1 ))
    fi

    fname=${files_dir}${varname}"_FGRN055_60-17_p"${part}".nc"
    fname_bu=${fname}"_backup"
    cp ${fname} ${fname_bu}
    
    fname_cor=${fname}"_corrected"
  
    ${ncap2} -s ${varname}2=${varname}*3600*3 ${fname} ${temp1}
    ${ncks} -x -v ${varname} ${temp1} ${temp2}
    ${ncrename} -v ${varname}2,${varname} ${temp2} ${temp3}
    mv ${temp3} ${fname_cor}
    rm ${temp1} ${temp2}

    mv ${fname_cor} ${fname}
  
    (( part = part + 1 ))
  done

#  fname=${ave_dir}${varname}"_FGRN055_60-79_ave.nc"
#  fname_bu=${fname}"_backup"
#  cp ${fname} ${fname_bu}
    
#  fname_cor=${fname}"_corrected"

#  ${ncap2} -s ${varname}2=${varname}*3600*3 ${fname} ${temp1}
#  ${ncks} -x -v ${varname} ${temp1} ${temp2}
#  ${ncrename} -v ${varname}2,${varname} ${temp2} ${temp3}
#  mv ${temp3} ${fname_cor}
#  rm ${temp1} ${temp2}

#  mv ${fname_cor} ${fname}

done

