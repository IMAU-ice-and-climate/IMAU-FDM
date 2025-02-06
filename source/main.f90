!------------------------------------------------------
!- MAIN PROGRAM OF THE FIRN DENSIFICATION MODEL
!--------------
!- Cleaned up to become IMAU-FDM v1.2
!- Cleaned up to become IMAU-FDM v1.0 (SL: 12-2014)
!- Adapted for ECMWF use (01-2012)
!- Made by: Stefan Ligtenberg (02-2010) 
!- based on the firn model of Michiel Helsen (2004-2007)
!
!- Current TO DOs
!   - initialise from restart
!   - fix age tracking through Year
!   - why are thread trackeres from distribute going to the restart folder?
!------------------------------------------------------

program main

    use openNetCDF, only: Load_Ave_Forcing, Load_TimeSeries_Forcing, Restart_From_Spinup, Restart_From_Run
    use output, only: Save_out_1D, Save_out_2D, Save_out_2Ddetail, Save_out_spinup, Save_out_run
    use initialise_variables, only: Define_Constants, Get_All_Command_Line_Arg, Get_Forcing_Dims, Get_Model_Settings, Calc_Output_Freq, Init_TimeStep_Var, Init_Prof_Var, &
    Init_Output_Var, Alloc_Forcing_Var
    use initialise_model, only: Init_Density_Prof, Init_Temp_Prof, Interpol_Forcing, Find_Grid, Index_Ave_Forcing
    use time_loop, only: Time_Loop_SpinUp, Time_Loop_Main
    
    implicit none
    
    integer :: ind_z_surf, ind_lon, ind_lat, Nt_forcing, Nlat, Nlon, Nt_model_interpol, Nt_model_tot, Nt_model_spinup, dtSnow
    integer :: ImpExp, dtmodel, dtmodelImp, dtmodelExp, dtobs
    integer :: writeinprof, writeinspeed, writeindetail
    integer :: numOutputProf, numOutputSpeed, numOutputDetail, outputProf, outputSpeed, outputDetail
    integer :: proflayers, detlayers, BeginT, startasice, IceShelf
    integer :: nyears, nyearsSU
    integer :: prev_nt
    integer, parameter :: ind_z_max = 20000
    
    character*255 :: username, point_numb, domain, fname_p1, project_name, restart_type
    
    double precision :: dzmax, initdepth, th, rho0_init, rhoi, R, pi, Ec, Eg, g, Lh, detthick
    double precision :: tsav, acav, ffav, lon_current, lat_current
    
    double precision, dimension(ind_z_max) :: Rho, M, T, Depth, Mlwc, DZ, DenRho, Refreeze, Year
    double precision, dimension(:), allocatable :: SnowMelt, PreTot, PreSol, PreLiq
    double precision, dimension(:), allocatable :: Sublim, SnowDrif, TempSurf, FF10m
    double precision, dimension(:), allocatable :: TempFM, PsolFM, PliqFM, SublFM, MeltFM, DrifFM, Rho0FM
    double precision, dimension(:,:), allocatable :: AveTsurf, AveAcc, AveWind, AveMelt
    double precision, dimension(:,:), allocatable :: ISM, LSM, Latitude, Longitude

    double precision, dimension(:,:), allocatable :: out_1D
    double precision, dimension(:,:), allocatable :: out_2D_dens, out_2D_temp, out_2D_lwc, out_2D_depth, out_2D_dRho, out_2D_year
    double precision, dimension(:,:), allocatable :: out_2D_det_dens, out_2D_det_temp, out_2D_det_lwc, out_2D_det_refreeze
    
    print *, " "
    print *, "------------------------------------"
    print *, "----- FIRN DENSIFICATION MODEL -----"
    print *, "------------------------------------"
    print *, " "

    ! TKTK RESTART IN PROGRESS
    
    call Get_All_Command_Line_Arg(username, point_numb, domain, fname_p1, project_name, restart_type)
    
    ! Read in the model settings, input settings and constants
    call Get_Model_Settings(dtSnow, nyears, nyearsSU, dtmodelImp, dtmodelExp, ImpExp, dtobs, ind_z_surf, startasice, &
        beginT, writeinprof, writeinspeed, writeindetail, proflayers, detlayers, detthick, dzmax, initdepth, th, &
        lon_current, lat_current, point_numb, username, domain, project_name)
    
    ! Read in resolution of the forcing data
    call Get_Forcing_Dims(Nlon, Nlat, Nt_forcing, domain, username)

    ! Load in physical and mathematical constants
    call Define_Constants(rhoi, R, pi, Ec, Eg, g, Lh)

    ! Determine model time step and amount of model time steps
    call Init_TimeStep_Var(dtobs, dtmodel, dtmodelImp, dtmodelExp, Nt_forcing, Nt_model_interpol, Nt_model_tot, Nt_model_spinup, ImpExp, nyearsSU)
    
    call Alloc_Forcing_Var(SnowMelt, PreTot, PreSol, PreLiq, Sublim, TempSurf, SnowDrif, FF10m, AveTsurf, LSM, ISM, &
        Latitude, Longitude, AveAcc, AveWind, AveMelt, TempFM, PsolFM, PliqFM, SublFM, MeltFM, DrifFM, Rho0FM, &
        Nt_forcing, Nt_model_tot, Nlon, Nlat)

    ! Get variables from the NetCDF files
    call Load_Ave_Forcing(AveTsurf, AveAcc, AveWind, AveMelt, LSM, Nlat, Nlon, Latitude, &
        Longitude, ISM, username, domain)
    
    print *, "Read all averaged values"
    print *, " "
    
    ! Find corresponding indices of the grid point for a given latitude and longitude
    call Find_Grid(ind_lon, ind_lat, lon_current, lat_current, Latitude, Longitude, LSM, Nlon, Nlat)

    ! Read averages for the current grid point		
    call Index_Ave_Forcing(AveTsurf, AveAcc, AveWind, AveMelt, ISM, tsav, acav, ffav, IceShelf, Nlon, Nlat, ind_lon, ind_lat)

    print *, "------ Point number: ", trim(point_numb), "------"
    print *, " Run for Lon: ", lon_current, " and Lat: ", lat_current
    print *, " Grid indices lon: ", ind_lon, ", lat: ", ind_lat
    print *, " Grounded (0) or Floating (1) ice: ", IceShelf
    print *, " Implicit (1) or Explicit (2) scheme: ", ImpExp
    print *, "------------------------------------"
    print *, " "
    
    call Load_TimeSeries_Forcing(SnowMelt, PreTot, PreSol, PreLiq, Sublim, SnowDrif, TempSurf, FF10m, Nt_forcing, ind_lon, ind_lat, username, &
        domain, dtobs)

    print *, "Got all variables from the NetCDF files"
    print *, " "

    call Init_Prof_Var(ind_z_surf, Rho, M, T, Depth, Mlwc, DZ, DenRho, Refreeze, Year, dzmax)

    ! Get variables needed for outputting data
    call Calc_Output_Freq(dtmodel, nyears, writeinprof, writeinspeed, writeindetail, numOutputProf, &
    numOutputSpeed, numOutputDetail, outputProf, outputSpeed, outputDetail)
    
    call Init_Output_Var(out_1D, out_2D_dens, out_2D_temp, out_2D_lwc, out_2D_depth, out_2D_dRho, out_2D_year, &
        out_2D_det_dens, out_2D_det_temp, out_2D_det_lwc, out_2D_det_refreeze, outputSpeed, outputProf, outputDetail, &
        proflayers, detlayers)
    
    ! Interpolate the RACMO forcing data to firn model time step
    call Interpol_Forcing(TempSurf, PreSol, PreLiq, Sublim, SnowMelt, SnowDrif, FF10m, TempFM, PsolFM, PliqFM, SublFM, &
        MeltFM, DrifFM, Rho0FM, Nt_forcing, Nt_model_interpol, Nt_model_tot, dtSnow, dtmodel, domain)
    

    if ( restart_type == "spinup" ) then
        call Restart_From_Spinup(ind_z_max, ind_z_surf, Rho, M, T, Depth, Mlwc, DZ, DenRho, Refreeze, username, point_numb, fname_p1, project_name)
        prev_nt=1
    else if ( restart_type == "none" ) then ! do spinup
        ! Construct an initial firn layer (T-, rho-, dz-, and M-profile)
        rho0_init = Rho0FM(1)
        call Init_Density_Prof(ind_z_max, ind_z_surf, dzmax, rho0_init, rhoi, R, Ec, Eg, g, acav, tsav, DZ, Rho, M, domain)

        call Init_Temp_Prof(ind_z_max, ind_z_surf, beginT, tsav, pi, T, Rho, Depth, rhoi)
    
        ! Spin up the model to a 'steady state'
        call Time_Loop_SpinUp(Nt_model_tot, Nt_model_spinup, ind_z_max, ind_z_surf, dtmodel, R, Ec, Eg, g, Lh, rhoi, acav, &
            th, dzmax, M, T, DZ, Rho, DenRho, Depth, Mlwc, Refreeze, Year, TempFM, PSolFM, PLiqFM, SublFM, MeltFM, DrifFM, Rho0FM, &
            IceShelf, ImpExp, nyears, nyearsSU, domain, project_name)

        ! Write intitial profile to NetCDF-file and prepare output arrays
        call Save_out_spinup(ind_z_max, ind_z_surf, Rho, M, T, Depth, Mlwc, Year, point_numb, fname_p1, username, project_name)
        prev_nt=1
    else if ( restart_type == "run" ) then ! start from previous run, so no spinup
        call Restart_From_Run(prev_nt, ind_z_max, ind_z_surf, Rho, M, T, Depth, Mlwc, DZ, Year, DenRho, Refreeze, username, &
                                point_numb, fname_p1, project_name)
    else
        print *, "Restart type not recognized: ", restart_type
    endif

    ! Call subprogram for spin-up and the time-integration			
    call Time_Loop_Main(dtmodel, ImpExp, Nt_model_tot, nyears, ind_z_max, ind_z_surf, numOutputSpeed, numOutputProf, numOutputDetail, &
        outputSpeed, outputProf, outputDetail, th, R, Ec, Eg, g, Lh, dzmax, rhoi, proflayers, detlayers, detthick, acav, IceShelf, &
        TempFM, PsolFM, PliqFM, SublFM, MeltFM, DrifFM, Rho0FM, Rho, M, T, Depth, Mlwc, DZ, DenRho, Refreeze, Year, &
        fname_p1, username, domain, out_1D, out_2D_dens, out_2D_temp, out_2D_lwc, out_2D_depth, out_2D_dRho, out_2D_year, &
        out_2D_det_dens, out_2D_det_temp, out_2D_det_lwc, out_2D_det_refreeze, prev_nt, restart_type)

    ! Write output to netcdf files
    call Save_out_1D(outputSpeed, point_numb, fname_p1, username, out_1D, project_name)
    call Save_out_2D(outputProf, proflayers, out_2D_dens, out_2D_temp, out_2D_lwc, out_2D_depth, out_2D_dRho, out_2D_year, &
        point_numb, fname_p1, username, project_name)
    call Save_out_2Ddetail(outputDetail, detlayers, detthick, out_2D_det_dens, out_2D_det_temp, out_2D_det_lwc, &
        out_2D_det_refreeze, point_numb, fname_p1, username, project_name)
    call Save_out_run(Nt_model_tot, ind_z_max, ind_z_surf, Rho, M, T, Depth, Mlwc, Year, &
        DenRho, Refreeze, point_numb, fname_p1, username, project_name)
    
    print *, "Written output data to files"
    
end program main
