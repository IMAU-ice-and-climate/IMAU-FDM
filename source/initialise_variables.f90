module initialise_variables
    !*** Subroutines for allocating memory and assigning initial values to variables ***!

     use model_settings

    implicit none

    private

    public :: Get_All_Command_Line_Arg, Get_Model_Settings_and_Forcing_Dimensions, Init_TimeStep_Var, Calc_Output_Freq, &
        Init_Prof_Var, Init_Output_Var, Alloc_Forcing_Var
    
contains


! *******************************************************


subroutine Get_All_Command_Line_Arg(username, point_numb, domain_name, prefix_output, project_name, restart_type)
    !*** Get all command line arguments ***!

    ! declare arguments
    character*255, intent(out) :: username, point_numb, domain_name, prefix_output, project_name, restart_type

    ! 1: ECMWF username 
    ! 2: Simulation number, corresponding to the line number in the IN-file.
    ! 3: Domain name, similar to the RACMO forcing (e.g. ANT27) - set to global in define_paths_and_constants
    ! 4: General part of output filename
    ! 5: Project name, for defining file paths
    ! 6: Restart type, where none= do spinup; spinup = restart from spinup

    call get_command_argument(1, username)
    call get_command_argument(2, point_numb)
    call get_command_argument(3, domain_name)
    call get_command_argument(4, prefix_output)
    call get_command_argument(5, project_name)
    call get_command_argument(6, restart_type)
    
end subroutine Get_All_Command_Line_Arg



!*** TKTKTK: move this to model_settings and make all variables global ***!

subroutine Get_Model_Settings_and_Forcing_Dimensions(dtSnow, nyears, nyearsSU, dtmodelImp, dtmodelExp, ImpExp, dtobs, &
    ind_z_surf, startasice, beginT, writeinprof, writeinspeed, writeindetail, proflayers, detlayers, detthick, dzmax, &
    initdepth, th, lon_current, lat_current, Nlon, Nlat, Nlon_timeseries, Nt_forcing)
    !*** Load model settings from file ***!
    
    ! declare arguments
    integer, intent(out) :: writeinprof, writeinspeed, writeindetail, dtSnow, nyears, nyearsSU, dtmodelImp, &
        dtmodelExp, ImpExp, dtobs, ind_z_surf, startasice, beginT, proflayers, detlayers, Nlon, Nlat, Nlon_timeseries, &
        Nt_forcing
    double precision, intent(out) :: dzmax, initdepth, th, lon_current, lat_current, detthick

    ! declare local variables
    integer :: NoR

    print *, "Path to input settings file:"
    print *, trim(path_settings)//trim(fname_settings)
    print *, " "

    open(unit=12, file=trim(path_settings)//trim(fname_settings))

    ! read model settings and forcing dimensions
    read (12,*)
    read (12,*)
    read (12,*) nyears              ! simulation time [yr]
    read (12,*) nyearsSU            ! simulation time during the spin up period [yr]
    read (12,*) dtmodelExp          ! time step in model with explicit T-scheme [s] 
    read (12,*) dtmodelImp          ! time step in model with implicit T-scheme [s]
    read (12,*) ImpExp              ! Impicit or Explicit scheme (1=Implicit/fast, 2= Explicit/slow)
    read (12,*) dtobs               ! time step in input data [s]
    read (12,*) dtSnow              ! duration of running average of variables used in snow parameterisation [s]
    read (12,*)
    read (12,*) dzmax               ! vertical model resolution [m]
    read (12,*) initdepth           ! initial depth of firn profile [m]
    read (12,*) th                  ! theta (if theta=0.5 , it is a Crank Nicolson scheme) 
    read (12,*) startasice          ! indicates the initial rho-profile (1=linear, 2=ice)
    read (12,*) beginT              ! indicates the inital T-profile (0=winter, 1=summer, 2=linear)
    read (12,*) NoR                 ! number of times the data series is repeated for the initial rho profile is constructed
    read (12,*) 
    read (12,*) writeinspeed        ! frequency of writing speed components to file
    read (12,*) writeinprof         ! frequency of writing of firn profiles to file
    read (12,*) proflayers          ! number of output layers in the prof file
    read (12,*) writeindetail       ! frequency of writing to detailed firn profiles to file
    read (12,*) detlayers           ! number of output layers in the detailed prof file
    read (12,*) detthick            ! thickness of output layer in the detailed prof file
    read (12,*)
    read (12,*) lon_current         ! Lon; indicates the longitude gridpoint
    read (12,*) lat_current         ! Lat; indicates the latitude gridpoint
    read (12,*)
    read (12,*) Nlon                 ! Number of longitude points in forcing
    read (12,*) Nlat                 ! Number of latitude points in forcing
    read (12,*) Nlon_timeseries           ! num of longitude bands (set during input pre-processing)
    read (12,*) Nt_forcing           ! Number of timesteps in forcing


    ind_z_surf = nint(initdepth/dzmax)  ! initial amount of vertical layers

    close(12)
    
    print *, "Loaded model settings"
    print *, " "

end subroutine Get_Model_Settings_and_Forcing_Dimensions


! *******************************************************


subroutine Init_TimeStep_Var(dtobs, dtmodel, dtmodelImp, dtmodelExp, Nt_forcing, Nt_model_interpol, Nt_model_tot, Nt_model_spinup, ImpExp, nyearsSU)
    !*** Initialise time step variables ***!

    ! declare arguments
    integer, intent(in) :: dtobs, dtmodelImp, dtmodelExp, nyearsSU, Nt_forcing, ImpExp
    integer, intent(out) :: Nt_model_interpol, Nt_model_tot, Nt_model_spinup
    integer, intent(inout) :: dtmodel

    if (ImpExp == 2) then
        dtmodel = dtmodelExp
    else
        dtmodel = dtmodelImp
    endif

    ! Number of IMAU-FDM time steps per forcing time step
    Nt_model_interpol = dtobs/dtmodel
    ! Total number of IMAU-FDM time steps
    Nt_model_tot = Nt_forcing*Nt_model_interpol
    ! Total number of time steps per spin-up period
    Nt_model_spinup = INT( (REAL(nyearsSU)*seconds_per_year)/(REAL(dtmodel)) )
    if (Nt_model_spinup > Nt_model_tot) Nt_model_spinup = Nt_model_tot

    print *, "dtmodel: ", dtmodel
    print *, 'Nt_model_interpol: ', Nt_model_interpol
    print *, 'Nt_model_tot: ', Nt_model_tot
    print *, 'Nt_model_spinup: ', Nt_model_spinup
    print *, ' '

end subroutine Init_TimeStep_Var


subroutine Calc_Output_Freq(dtmodel, nyears, writeinprof, writeinspeed, writeindetail, dtobs, Nt_forcing, numOutputProf, &
    numOutputSpeed, numOutputDetail, outputProf, outputSpeed, outputDetail)
    !*** Calculate the output frequency ***!

    ! declare arguments
    integer, intent(in) :: dtmodel, nyears, writeinprof, writeinspeed, writeindetail, dtobs, Nt_forcing
    integer, intent(out) :: numOutputProf, numOutputSpeed, numOutputDetail
    integer, intent(out) :: outputProf, outputSpeed, outputDetail
    integer, parameter :: double_kind = selected_real_kind( p=15, r=200 )

    integer(double_kind) :: outputProf2, outputProf2_num

    numOutputProf = writeinprof / dtmodel
    numOutputSpeed = writeinspeed / dtmodel
    numOutputDetail = writeindetail / dtmodel
    
    outputProf = INT(REAL(Nt_forcing * dtobs) / REAL(writeinprof))
    outputSpeed = INT(REAL(Nt_forcing * dtobs) / REAL(writeinspeed))
    outputDetail = INT(REAL(Nt_forcing * dtobs) / REAL(writeindetail))
    
    outputProf2_num = nyears * seconds_per_year
    outputProf2 = outputProf2_num / writeinprof

    print *, " "
    print *, "Output variables"
    print *, "Prof, Speed, Detail"
    print *, "writein...", writeinprof, writeinspeed, writeindetail
    print *, "numOutput...", numOutputProf, numOutputSpeed, numOutputDetail
    print *, "output...", outputProf, outputSpeed, outputDetail
    print *, " "

end subroutine Calc_Output_Freq


! *******************************************************


subroutine Init_Prof_Var(ind_z_surf, Rho, M, T, Depth, Mlwc, DZ, DenRho, Refreeze, Year, DZ_max)
    !*** Initialise firn profile variables with zeros ***!

    ! declare arguments
    integer, intent(in) :: ind_z_surf
    double precision, intent(in) :: DZ_max
    double precision, dimension(:), intent(out) :: Rho, M, T, Depth, Mlwc, DZ, DenRho, Refreeze, Year

    ! declare local variables
    integer :: ind_z

    Rho(:) = 0.
    M(:) = 0.
    T(:) = 0.
    Depth(:) = 0.
    Mlwc(:) = 0.
    DZ(:) = 0.
    DenRho(:) = 0.
    Refreeze(:) = 0.
    Year(:) = 0.
    
    DZ(1:ind_z_surf) = DZ_max

    Year(1:ind_z_surf) = -999.
    Depth(ind_z_surf) = 0.5*DZ_max
    do ind_z = (ind_z_surf-1), 1, -1
        Depth(ind_z) = Depth(ind_z+1) + 0.5*DZ(ind_z+1) + 0.5*DZ(ind_z)
    enddo

end subroutine Init_Prof_Var


! *******************************************************


subroutine Init_Output_Var(out_1D, out_2D_dens, out_2D_temp, out_2D_lwc, out_2D_depth, out_2D_dRho, out_2D_year, &
        out_2D_det_dens, out_2D_det_temp, out_2D_det_lwc, out_2D_det_refreeze, outputSpeed, outputProf, outputDetail, &
        proflayers, detlayers)
    !*** Initialise output variables with NaNs ***!

    ! declare arguments
    integer, intent(in) :: outputSpeed, outputProf, outputDetail, proflayers, detlayers
    double precision, dimension(:,:), allocatable, intent(out) :: out_1D, out_2D_dens, out_2D_temp, out_2D_lwc, out_2D_depth, out_2D_dRho, &
        out_2D_year, out_2D_det_dens, out_2D_det_temp, out_2D_det_lwc, out_2D_det_refreeze

    ! allocate memory to variables
    allocate(out_1D((outputSpeed), 18))
    allocate(out_2D_dens((outputProf), proflayers))
    allocate(out_2D_temp((outputProf), proflayers))
    allocate(out_2D_lwc((outputProf), proflayers))
    allocate(out_2D_depth((outputProf), proflayers))
    allocate(out_2D_dRho((outputProf), proflayers))
    allocate(out_2D_year((outputProf), proflayers))
    allocate(out_2D_det_dens((outputDetail), detlayers))
    allocate(out_2D_det_temp((outputDetail), detlayers))
    allocate(out_2D_det_lwc((outputDetail), detlayers))
    allocate(out_2D_det_refreeze((outputDetail), detlayers))
    
    ! Set missing value
    out_1D(:,:) = NaN_value
    out_2D_dens(:,:) = NaN_value
    out_2D_temp(:,:) = NaN_value
    out_2D_lwc(:,:) = NaN_value
    out_2D_depth(:,:) = NaN_value
    out_2D_dRho(:,:) = NaN_value
    out_2D_year(:,:) = NaN_value
    out_2D_det_dens(:,:) = NaN_value
    out_2D_det_temp(:,:) = NaN_value
    out_2D_det_lwc(:,:) = NaN_value
    out_2D_det_refreeze(:,:) = NaN_value

end subroutine Init_Output_Var


! *******************************************************

subroutine Alloc_Forcing_Var(SnowMelt, PreTot, PreSol, PreLiq, Sublim, TempSurf, SnowDrif, FF10m, AveTsurf, LSM, ISM, &
    Latitude, Longitude, AveAcc, AveWind, AveMelt, TempFM, PsolFM, PliqFM, SublFM, MeltFM, DrifFM, Rho0FM, &
    Nt_forcing, Nt_model_tot, Nlon, Nlat)
    !*** Allocate memory to forcing variables ***!

    ! declare arguments
    integer, intent(in) :: Nt_forcing, Nt_model_tot, Nlon, Nlat
    double precision, dimension(:), allocatable, intent(out) :: SnowMelt, PreTot, PreSol, PreLiq, Sublim, TempSurf, SnowDrif, FF10m
    double precision, dimension(:), allocatable, intent(out) :: MeltFM, PsolFM, PliqFM, SublFM, TempFM, DrifFM, Rho0FM
    double precision, dimension(:,:), allocatable, intent(out) :: AveTsurf, LSM, ISM, Latitude, Longitude, AveAcc, AveWind, AveMelt

    ! Adjust matrices to correct dimensions of the netCDF files

    ! Forcing time series variables
    allocate(SnowMelt(Nt_forcing))
    allocate(PreTot(Nt_forcing))
    allocate(PreSol(Nt_forcing))
    allocate(PreLiq(Nt_forcing))
    allocate(Sublim(Nt_forcing))
    allocate(TempSurf(Nt_forcing))
    allocate(SnowDrif(Nt_forcing))
    allocate(FF10m(Nt_forcing))
    
    ! Averaged variables
    allocate(AveTsurf(Nlon,Nlat))
    allocate(LSM(Nlon,Nlat))
    allocate(ISM(Nlon,Nlat))
    allocate(Latitude(Nlon,Nlat))
    allocate(Longitude(Nlon,Nlat))
    allocate(AveAcc(Nlon,Nlat))
    allocate(AveWind(Nlon,Nlat))
    allocate(AveMelt(Nlon,Nlat))

    ! Interpolated forcing time series variables
    allocate(MeltFM(Nt_model_tot))
    allocate(PsolFM(Nt_model_tot))
    allocate(PliqFM(Nt_model_tot))
    allocate(SublFM(Nt_model_tot))
    allocate(TempFM(Nt_model_tot))
    allocate(DrifFM(Nt_model_tot))
    allocate(Rho0FM(Nt_model_tot))
    
end subroutine Alloc_Forcing_Var

    
end module initialise_variables
