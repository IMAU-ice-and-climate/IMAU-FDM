#!/bin/bash 

# load the right compiler environment
# prgenvswitchto cray

# Find the LAT, LON, and NOR for the current point ($1)
ccab=$1
exe_id=$2
cpoint=$3

line2=`head -$cpoint $p2input | tail -1`
lat=`echo $line2 | cut -d',' -f1`
lon=`echo $line2 | cut -d',' -f2`
nor=`echo $line2 | cut -d',' -f3`
ImpExp=`echo $line2 | cut -d',' -f4`

# Construct an initialisation file for the current point depending on the domain
#
MSscript=$p2ms/model_settings_${domain}_${cpoint}.txt
if [[ "$domain" == "ANT27" ]]; then

cat << EOS > $MSscript
!-	MODEL SETTINGS FOR THE FIRN DENSIFICATION MODEL
!----------------------------------------------------
42	! nyears; simulation time [yr]
42	! nyears; simulation time during the spin-up [yr]
180	! dtmodelExp; timestep in model with explicit T-scheme [s]
3600	! dtmodelImp; timestep in model with implicit T-scheme [s]
$ImpExp ! ImpExp; implicit or explicit scheme, based on melt or not. (1=Implicit, 2=Explicit)
10800	! dtobs; timestep in input data [s]
31557600  ! dtSnow; duration of the running average used in the snow parameterisation [s]

0.15	! dzmax; vertical model resolution [m]
1.	! initdepth; initial depth of firn profile [m]
0.5	! th; theta (if theta=0.5, it is a Crank Nicolson scheme) 
1 	! startasice; indicates the initial rho-profile (1=linear, 2=ice)
3	! begintT; indicates the inital T-profile (1=winter, 2=summer, 3=linear)
$nor	! numberrepeat; number of times the data series is repeated for the initial rho profile is constructed

864000	! writeinspeed; frequency of writing speed components to file [s]
2592000 ! writeinprof; frequency of writing of firn profiles to file [s]
3000	! proflayers; number of output layer in prof file 
864000  ! writeindetail; frequency of writing of detailed firn profiles to file [s]
500	! detlayers; number of output layers in detail file
0.04	! detthick; thickness of the output layer in detail file [m]

$lat    ! beginLon; indicates the begin longitude gridpoint
$lon    ! beginLat; indicates the begin latitude gridpoint
EOS

fi

if [[ "$domain" == "XPEN055" ]]; then

cat << EOS > $MSscript
!-	MODEL SETTINGS FOR THE FIRN DENSIFICATION MODEL
!----------------------------------------------------
38	! nyears; simulation time [yr]
38	! nyearsSU; simulation time during the spin-up [yr]
180	! dtmodelExp; timestep in model with explicit T-scheme [s]
3600	! dtmodelImp; timestep in model with implicit T-scheme [s]
$ImpExp ! ImpExp; implicit or explicit scheme, based on melt or not. (1=Implicit, 2=Explicit)
10800	! dtobs; timestep in input data [s]

0.15	! dzmax; vertical model resolution [m]
1.	! initdepth; initial depth of firn profile [m]
0.5	! th; theta (if theta=0.5 , it is a Crank Nicolson scheme) 
1 	! startasice; indicates the initial rho-profile (1=linear, 2=ice)
3	! begintT; indicates the inital T-profile (1=winter, 2=summer, 3=linear)
$nor	! numberrepeat; number of times the data series is repeated for the initial rho profile is constructed

864000	! writeinspeed; frequency of writing speed components to file (in seconds)
2592000 ! writeinprof; frequency of writing of firn profiles to file (in seconds)
3000	! proflayers; number of output layer in prof file
864000  ! writeindetail; frequency of writing of detailed firn profiles to file
500	! detlayers; number of output layers in detail file
0.04	! detthick; thickness of the output layer in detail file

$lat    ! beginLon; indicates the begin longitude gridpoint
$lon    ! beginLat; indicates the begin latitude gridpoint
EOS

fi

if [ $domain = "FGRN11" ] ; then

cat << EOS > $MSscript
!-	MODEL SETTINGS FOR THE FIRN DENSIFICATION MODEL
!----------------------------------------------------
65	! nyears; simulation time [yr]
20	! nyears; simulation time during the spin-up [yr]
180	! dtmodelExp; timestep in model with explicit T-scheme [s]
3600	! dtmodelImp; timestep in model with implicit T-scheme [s]
$ImpExp ! ImpExp; implicit or explicit scheme, based on melt or not. (1=Implicit, 2=Explicit)
10800	! dtobs; timestep in input data [s]

0.15	! dzmax; vertical model resolution [m]
1.	! initdepth; initial depth of firn profile [m]
0.5	! th; theta (if theta=0.5 , it is a Crank Nicolson scheme) 
1 	! startasice; indicates the initial rho-profile (1=linear, 2=ice)
3	! begintT; indicates the inital T-profile (1=winter, 2=summer, 3=linear)
$nor	! numberrepeat; number of times the data series is repeated for the initial rho profile is constructed

864000	! writeinspeed; frequency of writing speed components to file (in seconds)
2592000 ! writeinprof; frequency of writing of firn profiles to file (in seconds)
3000	! proflayers; number of output layer in prof file
864000  ! writeindetail; frequency of writing of detailed firn profiles to file
500	! detlayers; number of output layers in detail file
0.04	! detthick; thickness of the output layer in detail file

$lat    ! beginLon; indicates the begin longitude gridpoint
$lon    ! beginLat; indicates the begin latitude gridpoint
EOS

fi

if [[ "$domain" == "FGRN055" ]]; then

cat << EOS > $MSscript
!-	MODEL SETTINGS FOR THE FIRN DENSIFICATION MODEL
!----------------------------------------------------
64	! nyears; simulation time [yr]
20	! nyears; simulation time during the spin-up [yr]
180	! dtmodelExp; timestep in model with explicit T-scheme [s]
900	! dtmodelImp; timestep in model with implicit T-scheme [s]
$ImpExp ! ImpExp; implicit or explicit scheme, based on melt or not. (1=Implicit, 2=Explicit)
10800	! dtobs; timestep in input data [s]
31557600  ! dtSnow; duration of the running average used in the snow parameterisation [s]

0.15	! dzmax; vertical model resolution [m]
1.	! initdepth; initial depth of firn profile [m]
0.5	! th; theta (if theta=0.5 , it is a Crank Nicolson scheme) 
1 	! startasice; indicates the initial rho-profile (1=linear, 2=ice)
3	! begintT; indicates the inital T-profile (1=winter, 2=summer, 3=linear)
$nor	! numberrepeat; number of times the data series is repeated for the initial rho profile is constructed

86400	! writeinspeed; frequency of writing speed components to file (in seconds)
2592000 ! writeinprof; frequency of writing of firn profiles to file (in seconds)
3000	! proflayers; number of output layer in prof file
86400  ! writeindetail; frequency of writing of detailed firn profiles to file
500	! detlayers; number of output layers in detail file
0.04	! detthick; thickness of the output layer in detail file

$lat    ! beginLon; indicates the begin longitude gridpoint
$lon    ! beginLat; indicates the begin latitude gridpoint
EOS

fi

if [ $domain = "FGRN11_ssp585" ] ; then

cat << EOS > $MSscript
!-	MODEL SETTINGS FOR THE FIRN DENSIFICATION MODEL
!----------------------------------------------------
150	! nyears; simulation time [yr]
20	! nyears; simulation time during the spin-up [yr]
180	! dtmodelExp; timestep in model with explicit T-scheme [s]
900	! dtmodelImp; timestep in model with implicit T-scheme [s]
$ImpExp ! ImpExp; implicit or explicit scheme, based on melt or not. (1=Implicit, 2=Explicit)
10800	! dtobs; timestep in input data [s]
31557600  ! dtSnow; duration of the running average used in the snow parameterisation [s]

0.15	! dzmax; vertical model resolution [m]
1.	! initdepth; initial depth of firn profile [m]
0.5	! th; theta (if theta=0.5 , it is a Crank Nicolson scheme) 
1 	! startasice; indicates the initial rho-profile (1=linear, 2=ice)
3	! begintT; indicates the inital T-profile (1=winter, 2=summer, 3=linear)
$nor	! numberrepeat; number of times the data series is repeated for the initial rho profile is constructed

864000	! writeinspeed; frequency of writing speed components to file (in seconds)
2592000 ! writeinprof; frequency of writing of firn profiles to file (in seconds)
3000	! proflayers; number of output layer in prof file
864000  ! writeindetail; frequency of writing of detailed firn profiles to file
500	! detlayers; number of output layers in detail file
0.04	! detthick; thickness of the output layer in detail file

$lat    ! beginLon; indicates the begin longitude gridpoint
$lon    ! beginLat; indicates the begin latitude gridpoint
EOS

fi
if [[ "$domain" == "PAT055" ]]; then

cat << EOS > $MSscript
!-	MODEL SETTINGS FOR THE FIRN DENSIFICATION MODEL
!----------------------------------------------------
34	! nyears; simulation time [yr]
34	! nyears; simulation time during the spin-up [yr]
180	! dtmodelExp; timestep in model with explicit T-scheme [s]
3600	! dtmodelImp; timestep in model with implicit T-scheme [s]
$ImpExp ! ImpExp; implicit or explicit scheme, based on melt or not. (1=Implicit, 2=Explicit)
21600	! dtobs; timestep in input data [s]

0.15	! dzmax; vertical model resolution [m]
20.	! initdepth; initial depth of firn profile [m]
0.5	! th; theta (if theta=0.5 , it is a Crank Nicolson scheme) 
1 	! startasice; indicates the initial rho-profile (1=linear, 2=ice)
3	! begintT; indicates the inital T-profile (1=winter, 2=summer, 3=linear)
$nor	! numberrepeat; number of times the data series is repeated for the initial rho profile is constructed

864000	! writeinspeed; frequency of writing speed components to file (in seconds)
2592000 ! writeinprof; frequency of writing of firn profiles to file (in seconds)
3000	! proflayers; number of output layer in prof file
864000  ! writeindetail; frequency of writing of detailed firn profiles to file
500	! detlayers; number of output layers in detail file
0.04	! detthick; thickness of the output layer in detail file

$lat    ! beginLon; indicates the begin longitude gridpoint
$lon    ! beginLat; indicates the begin latitude gridpoint
EOS

fi

if [[ $domain == "XDML055" ]]; then

cat << EOS > $MSscript
!-	MODEL SETTINGS FOR THE FIRN DENSIFICATION MODEL
!----------------------------------------------------
36	! nyears; simulation time [yr]
36	! nyears; simulation time during the spin-up [yr]
180	! dtmodelExp; timestep in model with explicit T-scheme [s]
3600	! dtmodelImp; timestep in model with implicit T-scheme [s]
$ImpExp ! ImpExp; implicit or explicit scheme, based on melt or not. (1=Implicit, 2=Explicit)
10800	! dtobs; timestep in input data [s]

0.15	! dzmax; vertical model resolution [m]
20.	! initdepth; initial depth of firn profile [m]
0.5	! th; theta (if theta=0.5 , it is a Crank Nicolson scheme) 
1 	! startasice; indicates the initial rho-profile (1=linear, 2=ice)
3	! begintT; indicates the inital T-profile (1=winter, 2=summer, 3=linear)
$nor	! numberrepeat; number of times the data series is repeated for the initial rho profile is constructed

864000	! writeinspeed; frequency of writing speed components to file (in seconds)
2592000 ! writeinprof; frequency of writing of firn profiles to file (in seconds)
3000	! proflayers; number of output layer in prof file
864000  ! writeindetail; frequency of writing of detailed firn profiles to file
500	! detlayers; number of output layers in detail file
0.04	! detthick; thickness of the output layer in detail file

$lat    ! beginLon; indicates the begin longitude gridpoint
$lon    ! beginLat; indicates the begin latitude gridpoint
EOS

fi

if [[ $domain == "ASE055" ]]; then

cat << EOS > $MSscript
!-	MODEL SETTINGS FOR THE FIRN DENSIFICATION MODEL
!----------------------------------------------------
37	! nyears; simulation time [yr]
37	! nyears; simulation time during the spin-up [yr]
180	! dtmodelExp; timestep in model with explicit T-scheme [s]
3600	! dtmodelImp; timestep in model with implicit T-scheme [s]
$ImpExp ! ImpExp; implicit or explicit scheme, based on melt or not. (1=Implicit, 2=Explicit)
10800	! dtobs; timestep in input data [s]

0.15	! dzmax; vertical model resolution [m]
20.	! initdepth; initial depth of firn profile [m]
0.5	! th; theta (if theta=0.5 , it is a Crank Nicolson scheme) 
1 	! startasice; indicates the initial rho-profile (1=linear, 2=ice)
3	! begintT; indicates the inital T-profile (1=winter, 2=summer, 3=linear)
$nor	! numberrepeat; number of times the data series is repeated for the initial rho profile is constructed

864000	! writeinspeed; frequency of writing speed components to file (in seconds)
2592000 ! writeinprof; frequency of writing of firn profiles to file (in seconds)
3000	! proflayers; number of output layer in prof file
864000  ! writeindetail; frequency of writing of detailed firn profiles to file
500	! detlayers; number of output layers in detail file
0.04	! detthick; thickness of the output layer in detail file

$lat    ! beginLon; indicates the begin longitude gridpoint
$lon    ! beginLat; indicates the begin latitude gridpoint
EOS

fi

if [[ $domain = "DMIS055" ]]; then

cat << EOS > $MSscript
!-	MODEL SETTINGS FOR THE FIRN DENSIFICATION MODEL
!----------------------------------------------------
38	! nyears; simulation time [yr]
38	! nyearsSU; simulation time during the spin-up [yr]
180	! dtmodelExp; timestep in model with explicit T-scheme [s]
3600	! dtmodelImp; timestep in model with implicit T-scheme [s]
$ImpExp ! ImpExp; implicit or explicit scheme, based on melt or not. (1=Implicit, 2=Explicit)
10800	! dtobs; timestep in input data [s]

0.15	! dzmax; vertical model resolution [m]
20.	! initdepth; initial depth of firn profile [m]
0.5	! th; theta (if theta=0.5 , it is a Crank Nicolson scheme) 
1 	! startasice; indicates the initial rho-profile (1=linear, 2=ice)
3	! begintT; indicates the inital T-profile (1=winter, 2=summer, 3=linear)
$nor	! numberrepeat; number of times the data series is repeated for the initial rho profile is constructed

864000	! writeinspeed; frequency of writing speed components to file (in seconds)
2592000 ! writeinprof; frequency of writing of firn profiles to file (in seconds)
3000	! proflayers; number of output layer in prof file
864000  ! writeindetail; frequency of writing of detailed firn profiles to file
500	! detlayers; number of output layers in detail file
0.04	! detthick; thickness of the output layer in detail file

$lat    ! beginLon; indicates the begin longitude gridpoint
$lon    ! beginLat; indicates the begin latitude gridpoint
EOS

fi


# Run the model
log_fname=${p2logs}/log_IMAU-FDM_${ccab}_${cpoint}.out
echo "$(date +%c) ${EC_FARM_ID}: We launch the model for ${cpoint} with:"
echo "$exe_id $usern $cpoint $domain $filename_part1 $ini_filename &> ${log_fname}"
$exe_id $usern $cpoint $domain $filename_part1 $ini_filename &> ${log_fname}

echo "$(date +%c) ${EC_FARM_ID}: Model run complete, report back..."
# report back that we are ready
echo "$(date +%c) ${EC_FARM_ID}: IMAU_FDM has completed gridpoint $cpoint" > $readydir/$cpoint

echo "$(date +%c) ${EC_FARM_ID}: Ready: Exit"
exit 0

