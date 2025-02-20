!------------------------------------------------------
!- MAIN PROGRAM OF THE FIRN DENSIFICATION MODEL
!--------------
!- Cleaned up to become IMAU-FDM v1.0 (SL: 12-2014)
!- Adapted for ECMWF use (01-2012)
!- Made by: Stefan Ligtenberg (02-2010) 
!- based on the firn model of Michiel Helsen (2004-2007)
!------------------------------------------------------

program main

    use openNetCDF, only: Load_Ave_Forcing, Load_TimeSeries_Forcing
    use output, only: Write_Initial_Output, Save_out_1D, Save_out_2D, Save_out_2Ddetail, Save_out_restart
    use initialise_variables, only: Define_Constants, Get_All_Command_Line_Arg, Get_Forcing_Dims, Get_Model_Settings, Define_Paths, Calc_Output_Freq, &
        Init_TimeStep_Var, Init_Prof_Var, Init_Output_Var, Alloc_Forcing_Var
    use initialise_model, only: Init_Density_Prof, Init_Temp_Prof, Interpol_Forcing, Find_Grid, Index_Ave_Forcing
    use time_loop, only: Time_Loop_SpinUp, Time_Loop_Main
    
    implicit none
    
    integer :: ind_z_surf, ind_lon, ind_lat, Nt_forcing, Nlat, Nlon, Nt_model_interpol, Nt_model_tot, Nt_model_spinup, Nlat_timeseries, Nlon_timeseries, dtSnow
    integer :: ImpExp, dtmodel, dtmodelImp, dtmodelExp, dtobs
    integer :: writeinprof, writeinspeed, writeindetail
    integer :: numOutputProf, numOutputSpeed, numOutputDetail, outputProf, outputSpeed, outputDetail
    integer :: proflayers, detlayers, BeginT, startasice, IceShelf
    integer :: nyears, nyearsSU
    integer, parameter :: ind_z_max = 20000
    
    character*255 :: username, point_numb, domain, prefix_output, ini_fname, project_name
    character*255 :: path_settings, path_forcing_dims, path_forcing_mask, path_forcing_averages, path_forcing_timeseries
    character*255 :: path_in_restart, path_out_restart, path_out_ini, path_out_1d, path_out_2d, path_out_2ddet
    character*255 :: fname_settings, fname_mask, suffix_forcing_averages, prefix_forcing_timeseries, fname_out_restart, fname_out_ini, fname_out_1d, fname_out_2d, fname_out_2ddet
    character*255 :: fname_forcing_dims, fname_shelfmask
    
    double precision :: DZ_max, initdepth, th, rho0_init, rhoi, R, pi, Ec, Ec2, Eg, g, Lh, kg, detthick
    double precision :: tsav, acav, ffav, lon_current, lat_current
    
    double precision, dimension(ind_z_max) :: Rho, M, T, rgrain2, Depth, Mlwc, DZ, DenRho, Refreeze, Year
    
    double precision, dimension(:), allocatable :: SnowMelt, PreTot, PreSol, PreLiq
    double precision, dimension(:), allocatable :: Sublim, SnowDrif, TempSurf, FF10m
    double precision, dimension(:), allocatable :: TempFM, PsolFM, PliqFM, SublFM, MeltFM, DrifFM, Rho0FM
    double precision, dimension(:,:), allocatable :: AveTsurf, AveAcc, AveWind, AveSubl, AveSnowDrif, AveMelt
    double precision, dimension(:,:), allocatable :: ISM, LSM, Latitude, Longitude
    double precision, dimension(:,:), allocatable :: out_1D
    double precision, dimension(:,:), allocatable :: out_2D_dens, out_2D_temp, out_2D_lwc, out_2D_depth, out_2D_dRho, out_2D_year, out_2D_rgrain
    double precision, dimension(:,:), allocatable :: out_2D_det_dens, out_2D_det_temp, out_2D_det_lwc, out_2D_det_refreeze, out_2D_det_rgrain

    print *, " "
    print *, "------------------------------------"
    print *, "----- FIRN DENSIFICATION MODEL -----"
    print *, "------------------------------------"
    print *, " "
    
    call Get_All_Command_Line_Arg(username, point_numb, domain, prefix_output, ini_fname, project_name)

    call Define_Paths(username, prefix_output, point_numb, path_settings, path_forcing_dims, path_forcing_mask, path_forcing_averages, path_forcing_timeseries, &
        path_in_restart, path_out_restart, path_out_ini, path_out_1d, path_out_2d, path_out_2ddet, fname_settings, fname_forcing_dims, fname_mask, &
        suffix_forcing_averages, prefix_forcing_timeseries, fname_out_restart, fname_out_ini, fname_out_1d, fname_out_2d, fname_out_2ddet, &
        Nlat_timeseries, Nlon_timeseries, project_name, domain)
    
    ! Read in the model settings, input settings and constants
    call Get_Model_Settings(dtSnow, nyears, nyearsSU, dtmodelImp, dtmodelExp, ImpExp, dtobs, ind_z_surf, startasice, &
        beginT, writeinprof, writeinspeed, writeindetail, proflayers, detlayers, detthick, DZ_max, initdepth, th, &
        lon_current, lat_current, path_settings, fname_settings, project_name)
    
    ! Read in resolution of the forcing data
    call Get_Forcing_Dims(Nlon, Nlat, Nt_forcing, path_forcing_dims, fname_forcing_dims)

    ! Load in physical and mathematical constants
    call Define_Constants(rhoi, R, pi, Ec, Ec2, Eg, g, Lh, kg)

    ! Determine model time step and amount of model time steps
    call Init_TimeStep_Var(dtobs, dtmodel, dtmodelImp, dtmodelExp, Nt_forcing, Nt_model_interpol, Nt_model_tot, Nt_model_spinup, ImpExp, nyearsSU)
    
    call Alloc_Forcing_Var(SnowMelt, PreTot, PreSol, PreLiq, Sublim, TempSurf, SnowDrif, FF10m, LSM, ISM, &
        Latitude, Longitude, AveAcc, AveWind, AveMelt, AveTsurf, AveSubl, AveSnowDrif, TempFM, PsolFM, PliqFM, SublFM, MeltFM, DrifFM, Rho0FM, &
        Nt_forcing, Nt_model_tot, Nlon, Nlat)

    ! Get variables from the NetCDF files
    call Load_Ave_Forcing(AveTsurf, AveAcc, AveWind, AveSubl, AveSnowDrif, AveMelt, LSM, Nlat, Nlon, Latitude, Longitude, &
        ISM, domain, path_forcing_mask, fname_mask, fname_shelfmask, path_forcing_averages, suffix_forcing_averages)
    
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
    
    call Load_TimeSeries_Forcing(SnowMelt, PreTot, PreSol,PreLiq, Sublim, SnowDrif, TempSurf, FF10m, Nt_forcing, &
        ind_lon, ind_lat, path_forcing_timeseries, prefix_forcing_timeseries, dtobs, Nlon, Nlat_timeseries, Nlon_timeseries)

    print *, "Got all variables from the NetCDF files"
    print *, " "

    call Init_Prof_Var(ind_z_surf, Rho, M, T, rgrain2, Depth, Mlwc, DZ, DenRho, Refreeze, Year, DZ_max)

    ! Get variables needed for outputting data
    call Calc_Output_Freq(dtmodel, nyears, writeinprof, writeinspeed, writeindetail, numOutputProf, &
    numOutputSpeed, numOutputDetail, outputProf, outputSpeed, outputDetail)
    
    call Init_Output_Var(out_1D, out_2D_dens, out_2D_temp, out_2D_lwc, out_2D_depth, out_2D_dRho, out_2D_year, out_2D_rgrain, &
        out_2D_det_dens, out_2D_det_temp, out_2D_det_lwc, out_2D_det_refreeze, out_2D_det_rgrain, outputSpeed, outputProf, outputDetail, &
        proflayers, detlayers)
    
    ! Interpolate the RACMO forcing data to firn model time step
    call Interpol_Forcing(TempSurf, PreSol, PreLiq, Sublim, SnowMelt, SnowDrif, FF10m, TempFM, PsolFM, PliqFM, SublFM, &
        MeltFM, DrifFM, Rho0FM, Nt_forcing, Nt_model_interpol, Nt_model_tot, dtSnow, dtmodel, domain)
    
    ! Construct an initial firn layer
    rho0_init = Rho0FM(1)
    call Init_Density_Prof(ind_z_max, ind_z_surf, DZ_max, rho0_init, rhoi, R, Ec, Ec2, Eg, g, kg, acav, tsav, DZ, Rho, M, rgrain2)

    call Init_Temp_Prof(ind_z_max, ind_z_surf, beginT, tsav, pi, T, Rho, Depth, rhoi)

    ! Spin up the model to a 'steady state'
    call Time_Loop_SpinUp(Nt_model_tot, Nt_model_spinup, ind_z_max, ind_z_surf, dtmodel, R, Ec, Ec2, Eg, g, Lh, kg, rhoi, acav, th, DZ_max, M, T, rgrain2, DZ, Rho, &
        DenRho, Depth, Mlwc, Refreeze, Year, TempFM, PSolFM, PLiqFM, SublFM, MeltFM, DrifFM, Rho0FM, IceShelf, &
        ImpExp, nyears, nyearsSU, domain)

    ! Write intitial profile to NetCDF-file and prepare output arrays
    call Write_Initial_Output(ind_z_max, ind_z_surf, Rho, M, T, Depth, Mlwc, Year, path_out_ini, fname_out_ini)
    
    ! Call subprogram for spin-up and the time-integration
    call Time_Loop_Main(dtmodel, ImpExp, Nt_model_tot, nyears, ind_z_max, ind_z_surf, numOutputSpeed, numOutputProf, numOutputDetail, &
        outputSpeed, outputProf, outputDetail, th, R, Ec, Ec2, Eg, g, Lh, kg, DZ_max, rhoi, proflayers, detlayers, detthick, acav, IceShelf, &
        TempFM, PsolFM, PliqFM, SublFM, MeltFM, DrifFM, Rho0FM, Rho, M, T, rgrain2, Depth, Mlwc, DZ, DenRho, Refreeze, Year, username, domain, &
        out_1D, out_2D_dens, out_2D_temp, out_2D_lwc, out_2D_depth, out_2D_dRho, out_2D_year, out_2D_rgrain, out_2D_det_dens, out_2D_det_temp, out_2D_det_lwc, &
        out_2D_det_refreeze, out_2D_det_rgrain)
    
    ! Write output to netcdf files
    call Save_out_1D(outputSpeed, path_out_1d, fname_out_1d, out_1D)
    call Save_out_2D(outputProf, proflayers, out_2D_dens, out_2D_temp, out_2D_lwc, out_2D_depth, out_2D_dRho, &
        out_2D_year, out_2D_rgrain, path_out_2d, fname_out_2d)
    call Save_out_2Ddetail(outputDetail, detlayers, detthick, out_2D_det_dens, out_2D_det_temp, out_2D_det_lwc, &
        out_2D_det_refreeze, out_2D_det_rgrain, path_out_2ddet, fname_out_2ddet)
    call Save_out_restart(Nt_model_tot, ind_z_max, ind_z_surf, Rho, DenRho, M, T, Depth, Mlwc, Year, Refreeze, DZ, path_out_restart, fname_out_restart)
    
    print *, "Written output data to files"
    
end program main
