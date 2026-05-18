module initialise_variables
    !*** Subroutines for allocating memory and assigning initial values to variables ***!

    use model_settings

    implicit none

    private

    public :: Get_Model_Settings_and_Forcing_Dimensions, Init_TimeStep_Var, Calc_Output_Freq, &
        Init_Prof_Var, Init_Output_Var, Alloc_Forcing_Var
    
contains


!*** TKTKTK: move this to model_settings and make all variables global ***!


subroutine Get_Model_Settings_and_Forcing_Dimensions(dtobs, ind_z_surf, lon_current, lat_current, Nlon, Nlat, &
    Nlon_timeseries, Nt_forcing)

    !*** Load model settings from parsed config ***!

    ! declare arguments
    integer, intent(out) :: dtobs, ind_z_surf, Nlon, Nlat, Nlon_timeseries, Nt_forcing
    double precision, intent(out) :: lon_current, lat_current

    print *, "Loaded model settings"
    print *, " "

    if (.not. allocated(config)) then
        error stop 'Model settings must be loaded before initializing variables'
    end if

    dtobs = 10800               ! time step in input data [s]
    lon_current = -43.5083      ! indicates the longitude gridpoint  ##TODO: These vary, where to set this?
    lat_current = 60.1478       ! indicates the latitude gridpoint  ##TODO: These vary, where to set this?
    Nlon = 438                  ! Number of longitude points in forcing
    Nlat = 566                  ! Number of latitude points in forcing
    Nlon_timeseries = 6         ! num of longitude bands (set during input pre-processing)
    Nt_forcing = 246424         ! Number of timesteps in forcing
     
    ind_z_surf = nint(config%model_physics%initdepth/config%general_settings%dzmax)  ! initial amount of vertical layers

end subroutine Get_Model_Settings_and_Forcing_Dimensions


! *******************************************************


subroutine Init_TimeStep_Var(dtobs, dtmodel, Nt_forcing, Nt_model_interpol, Nt_model_tot, Nt_model_spinup)
    !*** Initialise time step variables ***!

    ! declare arguments
    integer, intent(in) :: dtobs, Nt_forcing
    integer, intent(out) :: Nt_model_interpol, Nt_model_tot, Nt_model_spinup
    integer, intent(inout) :: dtmodel

    if (config%general_settings%ImpExp == 2) then
        dtmodel = config%general_settings%dtmodelExp
    else
        dtmodel = config%general_settings%dtmodelImp
    endif

    ! Number of IMAU-FDM time steps per forcing time step
    Nt_model_interpol = dtobs/dtmodel
    ! Total number of IMAU-FDM time steps
    Nt_model_tot = Nt_forcing*Nt_model_interpol
    ! Total number of time steps per spin-up period
    Nt_model_spinup = INT( (REAL(config%general_settings%nyears_spinup)*seconds_per_year)/(REAL(dtmodel)) )
    if (Nt_model_spinup > Nt_model_tot) Nt_model_spinup = Nt_model_tot

    print *, "dtmodel: ", dtmodel
    print *, 'Nt_model_interpol: ', Nt_model_interpol
    print *, 'Nt_model_tot: ', Nt_model_tot
    print *, 'Nt_model_spinup: ', Nt_model_spinup
    print *, ' '

end subroutine Init_TimeStep_Var


subroutine Calc_Output_Freq(dtmodel, dtobs, Nt_forcing, numOutputProf, &
    numOutputSpeed, numOutputDetail, outputProf, outputSpeed, outputDetail)
    !*** Calculate the output frequency ***!

    ! declare arguments
    integer, intent(in) :: dtmodel, dtobs, Nt_forcing
    integer, intent(out) :: numOutputProf, numOutputSpeed, numOutputDetail
    integer, intent(out) :: outputProf, outputSpeed, outputDetail

    ! Declare 64-bit integers locally for calculations
    integer(kind=8) :: Nt_forcing_64, dtobs_64, writeinspeed_64, writeinprof_64, writeindetail_64
    integer(kind=8) :: tempOutputProf, tempOutputSpeed, tempOutputDetail

    ! Convert inputs to 64-bit 
    ! if integers are larger than 2147483647 (32-bit integer), 64-bit integers are needed which can go up to 9223372036854775807
    ! Nt_forcing * dtobs > 2147483647, so 64-bit conversion
    Nt_forcing_64 = int(Nt_forcing, kind=8)
    dtobs_64 = int(dtobs, kind=8)
    writeinspeed_64 = int(config%output_dimensions%writeinspeed, kind=8)
    writeinprof_64 = int(config%output_dimensions%writeinprof, kind=8)
    writeindetail_64 = int(config%output_dimensions%writeindetail, kind=8)

    ! Calculate at what timestep output is produced
    numOutputProf = config%output_dimensions%writeinprof / dtmodel
    numOutputSpeed = config%output_dimensions%writeinspeed / dtmodel
    numOutputDetail = config%output_dimensions%writeindetail / dtmodel

    ! Calculate outputs using 64-bit integers
    tempOutputProf = (Nt_forcing_64 * dtobs_64) / writeinprof_64
    tempOutputSpeed = (Nt_forcing_64 * dtobs_64) / writeinspeed_64
    tempOutputDetail = (Nt_forcing_64 * dtobs_64) / writeindetail_64

    ! Assign back to 32-bit outputs (assuming values fit)
    outputProf = int(tempOutputProf)
    outputSpeed = int(tempOutputSpeed)
    outputDetail = int(tempOutputDetail)

    print *, " "
    print *, "Output variables"
    print *, "Prof, Speed, Detail"
    print *, "writein...", config%output_dimensions%writeinprof, config%output_dimensions%writeinspeed, config%output_dimensions%writeindetail
    print *, "numOutput...", numOutputProf, numOutputSpeed, numOutputDetail
    print *, "output...", outputProf, outputSpeed, outputDetail
    print *, " "

end subroutine Calc_Output_Freq


! *******************************************************


subroutine Init_Prof_Var(ind_z_surf, Rho, M, T, Depth, Mlwc, DZ, DenRho, Refreeze, Year)
    !*** Initialise firn profile variables with zeros ***!

    ! declare arguments
    integer, intent(in) :: ind_z_surf
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
    
    DZ(1:ind_z_surf) = config%general_settings%dzmax

    Year(1:ind_z_surf) = -999.
    Depth(ind_z_surf) = 0.5*config%general_settings%dzmax
    do ind_z = (ind_z_surf-1), 1, -1
        Depth(ind_z) = Depth(ind_z+1) + 0.5*DZ(ind_z+1) + 0.5*DZ(ind_z)
    enddo

end subroutine Init_Prof_Var


! *******************************************************


subroutine Init_Output_Var(out_1D, out_2D_dens, out_2D_temp, out_2D_lwc, out_2D_depth, out_2D_dRho, out_2D_year, &
    out_2D_det_dens, out_2D_det_temp, out_2D_det_lwc, out_2D_det_refreeze, outputSpeed, outputProf, outputDetail)
    !*** Initialise output variables with NaNs ***!

    ! declare arguments
    integer, intent(in) :: outputSpeed, outputProf, outputDetail
    double precision, dimension(:,:), allocatable, intent(out) :: out_1D, out_2D_dens, out_2D_temp, out_2D_lwc, out_2D_depth, out_2D_dRho, &
        out_2D_year, out_2D_det_dens, out_2D_det_temp, out_2D_det_lwc, out_2D_det_refreeze

    ! allocate memory to variables
    allocate(out_1D((outputSpeed), 18))
    allocate(out_2D_dens((outputProf), config%output_dimensions%proflayers))
    allocate(out_2D_temp((outputProf), config%output_dimensions%proflayers))
    allocate(out_2D_lwc((outputProf), config%output_dimensions%proflayers))
    allocate(out_2D_depth((outputProf), config%output_dimensions%proflayers))
    allocate(out_2D_dRho((outputProf), config%output_dimensions%proflayers))
    allocate(out_2D_year((outputProf), config%output_dimensions%proflayers))
    allocate(out_2D_det_dens((outputDetail), config%output_dimensions%detlayers))
    allocate(out_2D_det_temp((outputDetail), config%output_dimensions%detlayers))
    allocate(out_2D_det_lwc((outputDetail), config%output_dimensions%detlayers))
    allocate(out_2D_det_refreeze((outputDetail), config%output_dimensions%detlayers))
    
    ! Set missing value
    out_1D(:,:) = const%NaN_value
    out_2D_dens(:,:) = const%NaN_value
    out_2D_temp(:,:) = const%NaN_value
    out_2D_lwc(:,:) = const%NaN_value
    out_2D_depth(:,:) = const%NaN_value
    out_2D_dRho(:,:) = const%NaN_value
    out_2D_year(:,:) = const%NaN_value
    out_2D_det_dens(:,:) = const%NaN_value
    out_2D_det_temp(:,:) = const%NaN_value
    out_2D_det_lwc(:,:) = const%NaN_value
    out_2D_det_refreeze(:,:) = const%NaN_value

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
