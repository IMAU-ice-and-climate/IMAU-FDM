module initialise_variables
    !*** Subroutines for allocating memory and assigning initial values to variables ***!

    use, intrinsic :: iso_fortran_env, only: stderr => error_unit
    use tomlf, only: toml_table, toml_parse, get_value

    use model_settings

    implicit none

    private

    public :: Init_TimeStep_Var, Calc_Output_Freq, &
        Init_Prof_Var, Init_Output_Var, Alloc_Forcing_Var
    
contains

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

    write(log_unit, *) "dtmodel: ", dtmodel
    write(log_unit, *) 'Nt_model_interpol: ', Nt_model_interpol
    write(log_unit, *) 'Nt_model_tot: ', Nt_model_tot
    write(log_unit, *) 'Nt_model_spinup: ', Nt_model_spinup
    write(log_unit, *) ' '

end subroutine Init_TimeStep_Var

! *******************************************************

subroutine Calc_Output_Freq(dtmodel, writeinprof, writeinspeed, writeindetail, dtobs, Nt_forcing, numOutputProf, &
    numOutputSpeed, numOutputDetail, outputProf, outputSpeed, outputDetail)
    !*** Calculate the output frequency ***!

    ! declare arguments
    integer, intent(in) :: dtmodel, writeinprof, writeinspeed, writeindetail, dtobs, Nt_forcing
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
    writeinspeed_64 = int(writeinspeed, kind=8)
    writeinprof_64 = int(writeinprof, kind=8)
    writeindetail_64 = int(writeindetail, kind=8)

    ! Calculate at what timestep output is produced
    numOutputProf = writeinprof / dtmodel
    numOutputSpeed = writeinspeed / dtmodel
    numOutputDetail = writeindetail / dtmodel

    ! Calculate outputs using 64-bit integers
    tempOutputProf = (Nt_forcing_64 * dtobs_64) / writeinprof_64
    tempOutputSpeed = (Nt_forcing_64 * dtobs_64) / writeinspeed_64
    tempOutputDetail = (Nt_forcing_64 * dtobs_64) / writeindetail_64

    ! Assign back to 32-bit outputs (assuming values fit)
    outputProf = int(tempOutputProf)
    outputSpeed = int(tempOutputSpeed)
    outputDetail = int(tempOutputDetail)

    write(log_unit, *) " "
    write(log_unit, *) "Output variables"
    write(log_unit, *) "Prof, Speed, Detail"
    write(log_unit, *) "writein...", writeinprof, writeinspeed, writeindetail
    write(log_unit, *) "numOutput...", numOutputProf, numOutputSpeed, numOutputDetail
    write(log_unit, *) "output...", outputProf, outputSpeed, outputDetail
    write(log_unit, *) " "

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
