module model_settings
    !*** Subroutines for setting paths and constants for global use throughout the model ***!

    use, intrinsic :: iso_fortran_env, only: stderr => error_unit
    use tomlf
    implicit none

    private
    

    public :: Define_Job, Define_Filenames, Define_Constants, Define_Settings

    ! Module-level variables that can be accessed by other modules
    
    ! read in from distributor.f90
    public :: point_numb
    public :: settings_dir
    public :: log_unit

    ! read in from toml files
    public :: domain, forcing, restart_type
    public :: project_name
    public :: model_version
    public :: code_dir, data_dir
    public :: rhoi, rho_ocean, Tmelt, NaN_value, R, pi, Ec, Eg, g, Lh, seconds_per_year, ts_minimum, det2d_minimum
    public :: do_MO_fit

    ! currently read in from toml, will be set using netcdf
    public :: start_ts_year, end_ts_year, start_ave_year, end_ave_year, model_first_timestep, model_last_timestep
    
    ! created in model_settings
    
    public :: reference_dir
    public :: input_dir, input_averages_dir, input_timeseries_dir
    public :: project_dir, restart_dir, output_dir
    public :: fname_mask
    public :: prefix_forcing_timeseries, suffix_forcing_timeseries, suffix_forcing_averages
    public :: fname_restart_from_spinup, prefix_fname_run, suffix_fname_run, fname_restart_from_previous_run
    public :: fname_out_1d, fname_out_2d, fname_out_2ddet
    public :: iceshelf_var

    ! Declare the module variables
    ! note
    !   Variables read directly from TOML via get_value → character(len=:), allocatable
    !   Everything else (paths built by concatenation, command-line args, hardcoded strings) → character(len=512)

    ! read in from distributor.f90
    character(len=512) :: point_numb, settings_dir
    integer :: log_unit = 6   ! 6 = stdout; set to a file unit by distributor worker before each Run_Model call
    
    ! read in from toml files
    character(len=:), allocatable :: domain, forcing, restart_type
    character(len=:), allocatable :: project_name
    character(len=:), allocatable :: model_version
    character(len=:), allocatable :: code_dir, data_dir
    logical :: do_MO_fit 
    double precision :: rhoi, rho_ocean, Tmelt, NaN_value, R, pi, Ec, Eg, g, Lh, seconds_per_year, ts_minimum, det2d_minimum, days_per_year

    ! currently read in from toml, will be set using netcdf
    character(len=:), allocatable :: start_ts_year, end_ts_year, start_ave_year, end_ave_year, model_first_timestep, model_last_timestep

    ! created in model_settings
    character(len=512) :: reference_dir
    character(len=512) :: input_dir, input_averages_dir, input_timeseries_dir
    character(len=512) :: project_dir, output_dir, restart_dir
    character(len=512) :: fname_mask
    character(len=512) :: prefix_forcing_timeseries, suffix_forcing_timeseries, suffix_forcing_averages
    character(len=512) :: fname_restart_from_spinup, prefix_fname_run, suffix_fname_run, fname_restart_from_previous_run
    character(len=512) :: fname_out_1d, fname_out_2d, fname_out_2ddet
    character(len=512) :: iceshelf_var
    
contains

subroutine Define_Settings(nyearsSU, nyears, dtobs, Nlon, Nlat, Nlon_timeseries, Nt_forcing, &
    initdepth, startasice, begintT, &
    ImpExp, dtmodelImp, dtmodelExp, dtSnow, dzmax, th, &
    writeinspeed, writeinprof, writeindetail, proflayers, detlayers, detthick, &
    ind_z_surf)

    !*** Replaced and extended initialise_variables/Get_Model_Settings_and_Forcing_Dimensions ***!
    !*** Define model settings and physics from model.toml ***!

    ! Output arguments
    integer,          intent(out) :: nyearsSU, nyears, dtobs, Nlon, Nlat, Nlon_timeseries, Nt_forcing
    double precision, intent(out) :: initdepth
    integer,          intent(out) :: startasice, begintT
    integer,          intent(out) :: ImpExp, dtmodelImp, dtmodelExp, dtSnow
    double precision, intent(out) :: dzmax, th
    integer,          intent(out) :: writeinspeed, writeinprof, writeindetail, proflayers, detlayers
    double precision, intent(out) :: detthick
    integer,          intent(out) :: ind_z_surf

    ! Local variables
    character(len=512)            :: settings_path_model
    type(toml_table), allocatable :: table
    type(toml_table), pointer     :: child

    settings_path_model = trim(settings_dir)//"model.toml"
    call Load_TOML(settings_path_model, table)

    call get_value(table, 'fitting', child, requested=.false.)
    if (associated(child)) then
        call get_value(child, 'DO_MO_FIT', do_MO_fit)
    end if

    call get_value(table, 'forcing', child, requested=.false.)
    if (associated(child)) then
        call get_value(child, 'nyears_spinup', nyearsSU)
        call get_value(child, 'nyears_forcing', nyears)
        call get_value(child, 'dtobs', dtobs)
        call get_value(child, 'Nlon', Nlon)
        call get_value(child, 'Nlat', Nlat)
        call get_value(child, 'Nlon_timeseries', Nlon_timeseries)
        call get_value(child, 'Nt_forcing', Nt_forcing)
    end if

    call get_value(table, 'initialization', child, requested=.false.)
    if (associated(child)) then
        call get_value(child, 'startasice', startasice)
        call get_value(child, 'begintT', begintT)
        call get_value(child, 'initdepth', initdepth)
    end if

    call get_value(table, 'model', child, requested=.false.)
    if (associated(child)) then
        call get_value(child, 'ImpExp', ImpExp)
        call get_value(child, 'dtmodelImp', dtmodelImp)
        call get_value(child, 'dtmodelExp', dtmodelExp)
        call get_value(child, 'dtSnow', dtSnow)
        call get_value(child, 'dzmax', dzmax)
        call get_value(child, 'th', th)
    end if

    call get_value(table, 'output', child, requested=.false.)
    if (associated(child)) then
        call get_value(child, 'writeinspeed', writeinspeed)
        call get_value(child, 'writeinprof', writeinprof)
        call get_value(child, 'writeindetail', writeindetail)
        call get_value(child, 'proflayers', proflayers)
        call get_value(child, 'detlayers', detlayers)
        call get_value(child, 'detthick', detthick)
    end if

    call get_value(table, 'thresholds', child, requested=.false.)
    if (associated(child)) then
        call get_value(child, 'ts_minimum', ts_minimum)
        call get_value(child, 'det2d_minimum', det2d_minimum)
    end if

    ind_z_surf = nint(initdepth/dzmax)  ! initial amount of vertical layers

    write(log_unit, *) "Loaded model settings"
    write(log_unit, *) " "

    deallocate(table)

end subroutine Define_Settings

subroutine Define_Job()
    !*** Read run.toml and build all directory paths ***!

    character(len=512)            :: settings_path_run
    type(toml_table), allocatable :: table
    type(toml_table), pointer     :: child

    settings_path_run = trim(settings_dir)//"run.toml"
    call Load_TOML(settings_path_run, table)

    call get_value(table, 'job', child, requested=.false.)
    if (associated(child)) then
        call get_value(child, 'project_name', project_name)
        call get_value(child, 'forcing', forcing)
        call get_value(child, 'domain', domain)
        call get_value(child, 'restart_type', restart_type)

    end if

    call get_value(table, 'metadata', child, requested=.false.)
    if (associated(child)) then
        call get_value(child, 'model_version', model_version)
        call get_value(child, 'start_ts_year', start_ts_year)
        call get_value(child, 'end_ts_year', end_ts_year)
        call get_value(child, 'start_ave_year', start_ave_year)
        call get_value(child, 'end_ave_year', end_ave_year)
        call get_value(child, 'model_first_timestep', model_first_timestep)
        call get_value(child, 'model_last_timestep', model_last_timestep)

    end if

    call get_value(table, 'directories', child, requested=.false.)
    if (associated(child)) then
        call get_value(child, 'code_dir', code_dir)
        call get_value(child, 'data_dir', data_dir)

    end if

    project_dir   = trim(data_dir)//trim(project_name)//"/"
    output_dir    = trim(project_dir)//"output/"
    restart_dir   = trim(project_dir)//"restart/"
    reference_dir = trim(code_dir)//"reference/"//trim(domain)//"/"
    input_dir     = trim(data_dir)//trim(domain)//"_"//trim(forcing)//"/input/"

    if (trim(project_name) == "example") then
        input_dir   = trim(data_dir)//"input/"
        restart_dir = trim(code_dir)//"example/restart/"
    end if

    input_averages_dir   = trim(input_dir)//"averages/"
    input_timeseries_dir = trim(input_dir)//"timeseries/"

    write(log_unit, *) "Path to IMAU-FDM: ", code_dir
    write(log_unit, *) "Path to project directory: ", output_dir
    write(log_unit, *) "Path to input directory: ", input_dir
    write(log_unit, *) "Path to restart directory: ", restart_dir

    deallocate(table)

end subroutine Define_Job


! *******************************************************


subroutine Define_Filenames()
    !*** Build all output, restart, and input filenames ***!

    character(len=512) :: domain_forcing

    domain_forcing = trim(domain)//"_"//trim(forcing)

    fname_mask                = trim(domain)//"_Masks.nc"
    prefix_forcing_timeseries = "_"//trim(domain_forcing)//"_"//trim(start_ts_year)//"-"//trim(end_ts_year)//"_p"
    suffix_forcing_timeseries = ".nc"
    suffix_forcing_averages   = "_"//trim(domain_forcing)//"-"//trim(start_ts_year)//"_"// &
                                trim(start_ave_year)//"-"//trim(end_ave_year)//"_ave.nc"

    fname_out_1d    = trim(domain_forcing)//"_1D_"//trim(point_numb)//".nc"
    fname_out_2d    = trim(domain_forcing)//"_2D_"//trim(point_numb)//".nc"
    fname_out_2ddet = trim(domain_forcing)//"_2Ddetail_"//trim(point_numb)//".nc"

    fname_restart_from_spinup       = trim(domain_forcing)//"_restart_from_spinup_"//trim(point_numb)//".nc"
    prefix_fname_run                = trim(domain_forcing)//"_initialize_from_"
    suffix_fname_run                = "_run_"//trim(point_numb)//".nc"
    fname_restart_from_previous_run = trim(prefix_fname_run)//trim(end_ts_year)//trim(suffix_fname_run)

    if (trim(project_name) == "example") then
        prefix_forcing_timeseries = "_"//trim(domain_forcing)//"_"//trim(start_ts_year)//"-"//trim(end_ts_year)
        suffix_forcing_timeseries = "_point_"//trim(point_numb)//".nc"
    end if

    if (domain == "ANT27") then
        iceshelf_var = "IceShelve"
    else
        write(log_unit, *) "No iceshelf_var defined for domain, ", domain
        write(log_unit, *) "so not loading ice shelf mask."
    end if


end subroutine Define_Filenames



! ******************************************************* 



subroutine Define_Constants()

    character(len=512)            :: settings_path_constants
    type(toml_table), allocatable :: table
    type(toml_table), pointer     :: child

    settings_path_constants = trim(settings_dir)//"constants.toml"
    call Load_TOML(settings_path_constants, table)

    ! find section.
    call get_value(table, 'constants', child, requested=.false.)
    
    if (associated(child)) then
        call get_value(child, 'R', R)                           ! gas constant [J mole-1 K-1]
        call get_value(child, 'g', g)                           ! gravitational acceleration [m s-2]
        call get_value(child, 'rhoi', rhoi)                     ! density of ice [kg m-3]
        call get_value(child, 'rho_ocean', rho_ocean)           ! density of ocean water [kg m-3]
        call get_value(child, 'Tmelt', Tmelt)                   ! melting point of ice [K]
        call get_value(child, 'Ec', Ec)                         ! activation energy for creep [J mole-1]
        call get_value(child, 'Eg', Eg)                         ! activation energy for grain growth [J mole-1]
        call get_value(child, 'Lh', Lh)                         ! latent heat of fusion [J kg-1]
        call get_value(child, 'days_per_year', days_per_year)   ! days per year [days]
        call get_value(child, 'NaN_value', NaN_value)           ! missing value for doubles as used in the NCL scripts

        write(log_unit, *) " "
        write(log_unit, *) "days_per_year: ", days_per_year

    end if

    deallocate(table)

     pi = asin(1.)*2                             ! pi = 3.1415
    seconds_per_year = 3600.*24.*days_per_year   ! seconds per year [s]

end subroutine Define_Constants


! *******************************************************


subroutine Load_TOML(file_path, table)
    !*** Open, parse, and validate a TOML file; stops on any error ***!
    character(len=*),                intent(in)  :: file_path
    type(toml_table), allocatable,   intent(out) :: table

    integer :: fu, rc
    logical :: file_exists

    inquire(file=file_path, exist=file_exists)
    if (.not. file_exists) then
        write(stderr, '("Error: TOML file ", a, " not found")') trim(file_path)
        stop
    end if

    open(action='read', file=file_path, iostat=rc, newunit=fu)
    if (rc /= 0) then
        write(stderr, '("Error: Reading TOML file ", a, " failed")') trim(file_path)
        stop
    end if

    call toml_parse(table, fu)
    close(fu)

    if (.not. allocated(table)) then
        write(stderr, '("Error: Parsing TOML file ", a, " failed")') trim(file_path)
        stop
    end if

end subroutine Load_TOML


end module model_settings
