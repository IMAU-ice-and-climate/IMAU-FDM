!------------------------------------------------------
!- MAIN PROGRAM OF THE FIRN DENSIFICATION MODEL
!--------------
!- 2026 updates: 
!      - Submodule of distributor (distributor.f90)
!- Cleaned up to become IMAU-FDM v1.2
!- Cleaned up to become IMAU-FDM v1.0 (SL: 12-2014)
!- Adapted for ECMWF use (01-2012)
!- Made by: Stefan Ligtenberg (02-2010) 
!- based on the firn model of Michiel Helsen (2004-2007)
!
!- Current TO DOs
!   - fix age tracking through Year
!------------------------------------------------------

module model_main

    use model_settings
    use openNetCDF, only: Load_Mask, Load_Ave_Forcing, Load_TimeSeries_Forcing, Restart_From_Spinup, Restart_From_Run
    use output, only: Save_out_1D, Save_out_2D, Save_out_2Ddetail, Save_out_spinup, Save_out_run
    use initialise_variables, only: Calc_Output_Freq, Init_TimeStep_Var, Init_Prof_Var, Init_Output_Var, Alloc_Forcing_Var
    use initialise_model, only: Init_Density_Prof, Init_Temp_Prof, Interpol_Forcing, Index_Ave_Forcing
    use time_loop, only: Time_Loop_SpinUp, Time_Loop_Main
    
    implicit none

contains

subroutine Run_Model(cur_lat, cur_lon, ind_lat, ind_lon)

    double precision, intent(in) :: cur_lat, cur_lon
    integer,          intent(in) :: ind_lat, ind_lon

    integer :: Nt_forcing, Nt_model_interpol, Nt_model_tot, Nt_model_spinup
    integer :: writeinprof, writeinspeed, writeindetail
    integer :: numOutputProf, numOutputSpeed, numOutputDetail, outputProf, outputSpeed, outputDetail
    integer :: IceShelf
    integer :: dtmodel
    integer :: prev_nt

    double precision :: rho0_init, tsav, acav, ffav

    double precision, dimension(:), allocatable :: Rho, M, T, Depth, Mlwc, DZ, DenRho, Refreeze, Year
    double precision, dimension(:), allocatable :: SnowMelt, PreTot, PreSol, PreLiq
    double precision, dimension(:), allocatable :: Sublim, SnowDrif, TempSurf, FF10m
    double precision, dimension(:), allocatable :: TempFM, PsolFM, PliqFM, SublFM, MeltFM, DrifFM, Rho0FM
    double precision, dimension(:,:), allocatable :: AveTsurf, AveAcc, AveWind, AveMelt
    double precision, dimension(:,:), allocatable :: ISM, LSM, Latitude, Longitude

    double precision, dimension(:,:), allocatable :: out_1D
    double precision, dimension(:,:), allocatable :: out_2D_dens, out_2D_temp, out_2D_lwc, out_2D_depth, out_2D_dRho, out_2D_year
    double precision, dimension(:,:), allocatable :: out_2D_det_dens, out_2D_det_temp, out_2D_det_lwc, out_2D_det_refreeze

    write(log_unit, *) " "
    write(log_unit, *) "------------------------------------"
    write(log_unit, *) "----- FIRN DENSIFICATION MODEL -----"
    write(log_unit, *) "------------------------------------"
    write(log_unit, *) " "
    flush(log_unit)

    ! Define model settings and physics, runs first because decides whether or not to run as an example point
    call Load_Model_Settings()

    allocate(Rho(ind_z_max), M(ind_z_max), T(ind_z_max), Depth(ind_z_max), Mlwc(ind_z_max), &
             DZ(ind_z_max), DenRho(ind_z_max), Refreeze(ind_z_max), Year(ind_z_max))

    ! Build per-point output and restart filenames
    call Define_Filenames()

    ! Defines constants used throughout the model
    call Load_Constants()
    
    ! Determine model time step and amount of model time steps
    Nt_forcing = config%forcing_dimensions%Nt_forcing

    call Init_TimeStep_Var(dtmodel, Nt_forcing, Nt_model_interpol, Nt_model_tot, Nt_model_spinup)
    
    call Alloc_Forcing_Var(SnowMelt, PreTot, PreSol, PreLiq, Sublim, TempSurf, SnowDrif, FF10m, AveTsurf, LSM, ISM, &
        Latitude, Longitude, AveAcc, AveWind, AveMelt, TempFM, PsolFM, PliqFM, SublFM, MeltFM, DrifFM, Rho0FM, &
        Nt_forcing, Nt_model_tot)

    ! Get variables from the NetCDF files

    call Load_Mask(LSM, Latitude, Longitude, ISM)

    call Load_Ave_Forcing(AveTsurf, AveAcc, AveWind, AveMelt)
    
    write(log_unit, *) "Read all averaged values"
    write(log_unit, *) " "
    
    ! Read averages for the current grid point
    call Index_Ave_Forcing(AveTsurf, AveAcc, AveWind, AveMelt, ISM, tsav, acav, ffav, IceShelf, ind_lon, ind_lat)

    write(log_unit, *) "------ Point number: ", trim(point_numb), "------"
    write(log_unit, *) " Run for Lon: ", cur_lon, " and Lat: ", cur_lat
    write(log_unit, *) " Grid indices lon: ", ind_lon, ", lat: ", ind_lat
    write(log_unit, *) " Grounded (0) or Floating (1) ice: ", IceShelf
    write(log_unit, *) " Implicit (1) or Explicit (2) scheme: ", config%model_choices%ImpExp
    write(log_unit, *) "------------------------------------"
    write(log_unit, *) " "
    
    call Load_TimeSeries_Forcing(SnowMelt, PreTot, PreSol,PreLiq, Sublim, SnowDrif, TempSurf, FF10m, &
    ind_lon, ind_lat)

    write(log_unit, *) "Got all variables from the NetCDF files"
    write(log_unit, *) " "
    flush(log_unit)

    call Init_Prof_Var(Rho, M, T, Depth, Mlwc, DZ, DenRho, Refreeze, Year)

    ! Get variables needed for outputting data
    call Calc_Output_Freq(dtmodel, Nt_forcing, numOutputProf, numOutputSpeed, numOutputDetail, outputProf, outputSpeed, &
        outputDetail)

    call Init_Output_Var(out_1D, out_2D_dens, out_2D_temp, out_2D_lwc, out_2D_depth, out_2D_dRho, out_2D_year, &
        out_2D_det_dens, out_2D_det_temp, out_2D_det_lwc, out_2D_det_refreeze, outputSpeed, outputProf, outputDetail)
    
    ! Interpolate the RACMO forcing data to firn model time step
    call Interpol_Forcing(TempSurf, PreSol, PreLiq, Sublim, SnowMelt, SnowDrif, FF10m, TempFM, PsolFM, PliqFM, SublFM, &
        MeltFM, DrifFM, Rho0FM, Nt_forcing, Nt_model_interpol, Nt_model_tot, dtmodel)

    

    if ( restart_type == "spinup" ) then
        call Restart_From_Spinup(Rho, M, T, Depth, Mlwc, DZ, DenRho, Refreeze)
        prev_nt=1
    else if ( restart_type == "none" ) then ! do spinup
        ! Construct an initial firn layer (T-, rho-, dz-, and M-profile)
        rho0_init = Rho0FM(1)
        call Init_Density_Prof(rho0_init, acav, tsav, DZ, Rho, M)

        call Init_Temp_Prof(tsav, T, Rho, Depth)

        ! Spin up the model to a 'steady state'
        call Time_Loop_SpinUp(Nt_model_tot, Nt_model_spinup, dtmodel, acav, ffav, M, T, DZ, Rho, DenRho, Depth, &
            Mlwc, Refreeze, Year, TempFM, PSolFM, PLiqFM, SublFM, MeltFM, DrifFM, Rho0FM, &
            IceShelf)

        ! Write intitial profile to NetCDF-file and prepare output arrays
        call Save_out_spinup(Rho, M, T, Depth, Mlwc, Year)
        prev_nt=1

    else if ( restart_type == "run" ) then ! start from previous run, so no spinup

        call Restart_From_Run(prev_nt, Rho, M, T, Depth, Mlwc, DZ, Year, DenRho, Refreeze)
    else
        write(log_unit, *) "Restart type not recognized: ", restart_type

    endif

    flush(log_unit)
    ! Call subprogram for spin-up and the time-integration
    call Time_Loop_Main(dtmodel, Nt_model_tot, numOutputSpeed, numOutputProf, numOutputDetail, &
        outputSpeed, outputProf, outputDetail, acav, ffav, IceShelf, TempFM, PsolFM, PliqFM, SublFM, MeltFM, DrifFM, Rho0FM, &
        Rho, M, T, Depth, Mlwc, DZ, DenRho, Refreeze, Year, out_1D, out_2D_dens, out_2D_temp, out_2D_lwc, out_2D_depth, out_2D_dRho, &
        out_2D_year, out_2D_det_dens, out_2D_det_temp, out_2D_det_lwc, out_2D_det_refreeze, prev_nt)

    ! Write output to netcdf files
    call Save_out_1D(outputSpeed, out_1D)
    call Save_out_2D(out_2D_dens, out_2D_temp, out_2D_lwc, out_2D_depth, out_2D_dRho, out_2D_year)
    call Save_out_2Ddetail(out_2D_det_dens, out_2D_det_temp, out_2D_det_lwc, out_2D_det_refreeze)
    call Save_out_run(Nt_model_tot, Rho, M, T, Depth, Mlwc, Year, DenRho, Refreeze)
    
    write(log_unit, *) "Written output data to files"
    flush(log_unit)

end subroutine Run_Model
    

! ******************************************************* 


end module model_main
