module initialise_variables
    !*** Subroutines for allocating memory and assigning initial values to variables ***!

    implicit none
    private

    public :: Define_Constants, Get_All_Command_Line_Arg, Get_Model_Settings_and_Forcing_Dimensions, Init_TimeStep_Var, Calc_Output_Freq, Init_Prof_Var, &
        Init_Output_Var, Alloc_Forcing_Var
    
contains


! *******************************************************


subroutine Define_Constants(rhoi, R, pi, Ec, Eg, g, Lh)
    !*** Define physical constants ***!
    
    double precision, intent(out) :: rhoi, R, pi, Ec, Eg, g, Lh
    
    pi = asin(1.)*2         ! pi = 3.1415
    R = 8.3145              ! gas constant [J mole-1 K-1]
    g = 9.81                ! gravitational acceleration [m s-2]
    rhoi = 917.             ! density of ice [kg m-3]
    Ec = 60000.             ! activation energy for creep [J mole-1]
    Eg = 42400.             ! activation energy for grain growth [J mole-1]
    Lh = 333500.            ! latent heat of fusion [J kg-1]

end subroutine Define_Constants


! *******************************************************


subroutine Get_All_Command_Line_Arg(username, point_numb, domain, fname_p1, project_name, restart_type)
    !*** Get all command line arguments ***!

    ! declare arguments
    character*255, intent(out) :: username, point_numb, domain, fname_p1, project_name, restart_type

    ! 1: ECMWF username (e.g. nmg)
    ! 2: Simulation number, corresponding to the line number in the IN-file.
    ! 3: Domain name, similar to the RACMO forcing (e.g. ANT27)
    ! 4: General part of output filename
    ! 5: Optional, name of the initialization file
    ! 6: Path to output - **TKTK IN PROGRESS**
    call get_command_argument(1, username)
    call get_command_argument(2, point_numb)
    call get_command_argument(3, domain)
    call get_command_argument(4, fname_p1)
    call get_command_argument(5, project_name)
    call get_command_argument(6, restart_type)
    
    
end subroutine Get_All_Command_Line_Arg


! *******************************************************

! Decomissioned on 2025-02-06, input_settings values moved to start_model_ccab file and read in in "Get_Model_Settings"

! subroutine Get_Forcing_Dims(Nlon, Nlat, Nt_forcing, domain, username)

!     !*** Load dimensions of forcing data from file ***!
    
!     ! declare arguments
!     integer, intent(out) :: Nlon, Nlat, Nt_forcing
!     character*255 :: pad
!     character*255, intent(in) :: domain, username

!     !pad = ""//trim(cwd)//"/../reference/"//trim(domain)//"/input_settings_"//trim(domain)//".txt"

!     pad = "/perm/"//trim(username)//"/code/IMAU-FDM/reference/"//trim(domain)//"/input_settings_"//trim(domain)//".txt"

!     print *, "Path to model dimensions file:"
!     print *, trim(pad)
    
!     open(unit=12,file=trim(pad))

!     read(12,*)
!     read(12,*)
!     read(12,*) Nlon
!     read(12,*) Nlat
!     read(12,*) Nt_forcing

!     print *, "Read input settings"
    
!     close(12)

! end subroutine Get_Forcing_Dims


! *******************************************************


subroutine Get_Model_Settings_and_Forcing_Dimensions(dtSnow, nyears, nyearsSU, dtmodelImp, dtmodelExp, ImpExp, dtobs, ind_z_surf, startasice, &
    beginT, writeinprof, writeinspeed, writeindetail, proflayers, detlayers, detthick, dzmax, initdepth, th, &
    lon_current, lat_current, Nlon, Nlat, Nt_forcing, point_numb, username, domain, project_name)
    !*** Load model settings from file ***!
    
    ! declare arguments
    integer, intent(out) :: writeinprof, writeinspeed, writeindetail, dtSnow, nyears, nyearsSU, dtmodelImp, &
        dtmodelExp, ImpExp, dtobs, ind_z_surf, startasice, beginT, proflayers, detlayers, Nlon, Nlat, Nt_forcing
    double precision, intent(out) :: dzmax, initdepth, th, lon_current, lat_current, detthick
    character*255, intent(in) :: point_numb, username, domain, project_name

    ! declare local variables
    character*255 :: pad
    integer :: NoR

    pad = "/ec/res4/scratch/"//trim(username)//"/"//trim(project_name)//"/ms_files/"
    
    print *, "Path to input settings file:"
    print *, trim(pad)//"model_settings_"//trim(domain)//"_"//trim(point_numb)//".txt"
    print *, " "

    open(unit=12, file=trim(pad)//"model_settings_"//trim(domain)//"_"//trim(point_numb)//".txt")

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
    if (domain .NE. "none") read (12,*) lon_current         ! Lon; indicates the longitude gridpoint
    if (domain .NE. "none") read (12,*) lat_current         ! Lat; indicates the latitude gridpoint
    read (12,*)
    if (domain .NE. "none") read (12,*) Nlon                 ! Number of longitude points in forcing
    if (domain .NE. "none") read (12,*) Nlat                 ! Number of latitude points in forcing
    read (12,*) Nt_forcing           ! Number of timesteps in forcing


    ind_z_surf = nint(initdepth/dzmax)  ! initial amount of vertical layers

    close(11)
    
    print *, "Loaded model settings"
    print *, " "

end subroutine Get_Model_Settings_and_Forcing_Dimensions


! *******************************************************


subroutine Init_TimeStep_Var(dtobs, dtmodel, dtmodelImp, dtmodelExp, Nt_forcing, Nt_model_interpol, Nt_model_tot, Nt_model_spinup, ImpExp, nyearsSU)
    !*** Initialise time step variables ***!

    ! declare arguments
    integer, intent(in) :: dtmodelImp, dtmodelExp, nyearsSU, ImpExp, Nt_forcing, dtobs
    integer, intent(out) :: Nt_model_interpol, Nt_model_tot, Nt_model_spinup
    integer, intent(inout) :: dtmodel

    if (ImpExp == 2) then
        dtmodel = dtmodelExp   !You could use this for dry locations in Antarctica, but should not be called dtmodelExp because explicit was meant for wet locations
    else
        dtmodel = dtmodelImp
    endif

    ! Number of IMAU-FDM time steps per forcing time step
    Nt_model_interpol = dtobs/dtmodel
    ! Total number of IMAU-FDM time steps
    Nt_model_tot = Nt_forcing*Nt_model_interpol
    ! Total number of time steps per spin-up period
    ! Nt_model_spinup = INT( REAL(Nt_model_tot)*(REAL(nyearsSU)*365.25*24.*3600.)/REAL((REAL(Nt_forcing)*REAL(dtobs))) )
    Nt_model_spinup = INT ((REAL(nyearsSU)*365.25*24.*3600.))/REAL(dtmodel)
    
    if (Nt_model_spinup > Nt_model_tot) then
        Nt_model_spinup = Nt_model_tot
    elseif (Nt_model_spinup < 0.) then
        Nt_model_spinup = Nt_model_tot
    endif

    print *, "dtmodel: ", dtmodel
    print *, 'Nt_model_interpol: ', Nt_model_interpol
    print *, 'Nt_model_tot: ', Nt_model_tot
    print *, 'Nt_model_spinup: ', Nt_model_spinup
    print *, ' '

end subroutine Init_TimeStep_Var


subroutine Calc_Output_Freq(dtmodel, nyears, writeinprof, writeinspeed, writeindetail, numOutputProf, &
    numOutputSpeed, numOutputDetail, outputProf, outputSpeed, outputDetail)
    !*** Calculate the output frequency ***!

    ! declare arguments
    integer, intent(in) :: dtmodel, nyears, writeinprof, writeinspeed, writeindetail
    integer, intent(out) :: numOutputProf, numOutputSpeed, numOutputDetail
    integer, intent(out) :: outputProf, outputSpeed, outputDetail
    integer, parameter :: double_kind = selected_real_kind( p=15, r=200 ) !a bit over the top
    integer(double_kind) :: spy ! prevents overflow when calculating output writing intervals
    
    spy = 3600.*24.*365.
    
    numOutputProf = writeinprof / dtmodel
    numOutputSpeed = writeinspeed / dtmodel
    numOutputDetail = writeindetail / dtmodel

    outputProf = INT((REAL(nyears * spy) / REAL(writeinprof)))
    outputSpeed = INT((REAL(nyears * spy) / REAL(writeinspeed)))
    outputDetail = INT((REAL(nyears * spy) / REAL(writeindetail)))

    print *, " "
    print *, "Output variables"
    print *, "Prof, Speed, Detail"
    print *, "writein...", writeinprof, writeinspeed, writeindetail
    print *, "numOutput...", numOutputProf, numOutputSpeed, numOutputDetail
    print *, "output...", outputProf, outputSpeed, outputDetail
    print *, " "

end subroutine Calc_Output_Freq


! *******************************************************


subroutine Init_Prof_Var(ind_z_surf, Rho, M, T, Depth, Mlwc, DZ, DenRho, Refreeze, Year, dzmax)
    !*** Initialise firn profile variables with zeros ***!

    ! declare arguments
    integer, intent(in) :: ind_z_surf
    double precision, intent(in) :: dzmax
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
    
    DZ(1:ind_z_surf) = dzmax
    Year(1:ind_z_surf) = -999.
    Depth(ind_z_surf) = 0.5*dzmax
    do ind_z = (ind_z_surf-1), 1, -1
        Depth(ind_z) = Depth(ind_z+1) + 0.5*DZ(ind_z+1) + 0.5*DZ(ind_z)
    enddo

end subroutine Init_Prof_Var


! *******************************************************


subroutine Init_Output_Var(out_1D, out_2D_dens, out_2D_temp, out_2D_lwc, out_2D_depth, &
    out_2D_dRho, out_2D_year, out_2D_det_dens, out_2D_det_temp, out_2D_det_lwc, &
    out_2D_det_refreeze, outputSpeed, outputProf, outputDetail, proflayers, detlayers)
    !*** Initialise output variables with NaNs ***!

    ! declare arguments
    integer, intent(in) :: outputSpeed, outputProf, outputDetail, proflayers, detlayers
    double precision, dimension(:,:), allocatable, intent(out) :: out_1D, out_2D_dens, out_2D_temp, out_2D_lwc, &
        out_2D_depth, out_2D_dRho, out_2D_year, out_2D_det_dens, out_2D_det_temp, out_2D_det_lwc, out_2D_det_refreeze

    ! declare local variables
    double precision :: NaN_value

    ! This is the missing value for doubles as used in the NCL scripts
    NaN_value = 9.96921e+36

    ! allocate memory to variables
    allocate(out_1D((outputSpeed+50), 18))
    allocate(out_2D_dens((outputProf+50), proflayers))
    allocate(out_2D_temp((outputProf+50), proflayers))
    allocate(out_2D_lwc((outputProf+50), proflayers))
    allocate(out_2D_depth((outputProf+50), proflayers))
    allocate(out_2D_dRho((outputProf+50), proflayers))
    allocate(out_2D_year((outputProf+50), proflayers))
    allocate(out_2D_det_dens((outputDetail+50), detlayers))
    allocate(out_2D_det_temp((outputDetail+50), detlayers))
    allocate(out_2D_det_lwc((outputDetail+50), detlayers))
    allocate(out_2D_det_refreeze((outputDetail+50), detlayers))
    
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
    Nt_forcing, Nt_model_tot, Nlon, Nlat, domain)
    !*** Allocate memory to forcing variables ***!

    ! declare arguments
    integer, intent(in) :: Nt_forcing, Nt_model_tot, Nlon, Nlat
    double precision, dimension(:), allocatable, intent(out) :: SnowMelt, PreTot, PreSol, PreLiq, Sublim, TempSurf, SnowDrif, FF10m
    double precision, dimension(:), allocatable, intent(out) :: MeltFM, PsolFM, PliqFM, SublFM, TempFM, DrifFM, Rho0FM
    double precision, dimension(:,:), allocatable, intent(out) :: AveTsurf, LSM, ISM, Latitude, Longitude, AveAcc, AveWind, AveMelt
    character*255 :: domain
    
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
    
    if (domain .NE. "none") then

    ! Averaged variables
    allocate(AveTsurf(Nlon,Nlat))
    allocate(LSM(Nlon,Nlat))
    allocate(ISM(Nlon,Nlat))
    allocate(Latitude(Nlon,Nlat))
    allocate(Longitude(Nlon,Nlat))
    allocate(AveAcc(Nlon,Nlat))
    allocate(AveWind(Nlon,Nlat))
    allocate(AveMelt(Nlon,Nlat))

    endif

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
