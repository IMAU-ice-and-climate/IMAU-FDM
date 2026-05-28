module model_settings
    !*** Subroutines for setting paths and constants for global use throughout the model ***!

    use, intrinsic :: iso_fortran_env, only: stderr => error_unit
    use tomlf
    implicit none

    private ! none?
    
    type, public :: Constants
        double precision :: R  ! gas constant [J mole-1 K-1]
        double precision :: g  ! gravitational acceleration [m s-2]
        double precision :: rhoi  ! density of ice [kg m-3]
        double precision :: rho_ocean  ! density of ocean water [kg m-3]
        double precision :: Tmelt  ! melting point of ice [K]
        double precision :: Ec  ! activation energy for creep [J mole-1]
        double precision :: Eg  ! activation energy for grain growth [J mole-1]
        double precision :: Lh  ! latent heat of fusion [J kg-1]
        double precision :: days_per_year  ! days per year TKTKTK: remove once timestamps incorprated
        double precision :: NaN_value  ! missing value for doubles as used in the NCL scripts
    end type Constants

    type, public :: forcing_dimensions_t

        integer :: nyears_spinup ! number of years spinup is performed for
        integer :: nyears_forcing ! number of years of forcing
        integer :: dtobs !length of forcing timestep
        integer :: Nlon ! size of longitude grid
        integer :: Nlat ! size of latitude grid
        integer :: Nlon_timeseries ! number of longitude bands in each parts file
        integer :: Nt_forcing !number of timesteps in forcing

    end type forcing_dimensions_t

    type, public :: initialization_t

        integer :: startasice
        integer :: beginT
        double precision :: initdepth

    end type initialization_t

    type, public :: model_choices_t

        integer :: ImpExp
        integer :: dtmodelImp
        integer :: dtmodelExp
        integer :: dtSnow
        double precision :: dzmax
        double precision :: th
        double precision :: ts_minimum
        double precision :: det2d_minimum

    end type model_choices

    type, public :: output_dimensions_t
        integer :: writeinspeed
        integer :: writeinprof
        integer :: proflayers
        integer :: writeindetail
        integer :: detlayers
        double precision :: detthick
    end type output_dimensions_t

    type, public :: run_settings_t

        character(len=:), allocatable :: restart_type
        character(len=:), allocatable :: domain
        character(len=:), allocatable :: forcing

    end type general_settings_t

    type, public :: Settings !TKTKTK explain and what about Constants?
        type(forcing_t)  :: forcing_dimensions
        type(initialization_t)  :: initialization
        type(model_choices_t)  :: model_choices
        type(output_dimensions_t)  :: output_dimensions
        type(run_settings_t):: run_settings
    end type Settings

    public :: Read_Job, Read_Constants, Load_Constants, Read_Settings, Load_Model_Settings, Define_Filenames

    ! Module-level variables that can be accessed by other modules
    
    ! set in distributor.f90
    public :: point_numb
    public :: log_unit

    ! subfolder of code_dir; defined in launch_job file
    public :: settings_dir

    ! set in run.tmol
    public :: project_name
    public :: model_version
    public :: code_dir, data_dir

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

    public :: pi, seconds_per_year
    public :: config, const

    ! Declare the module variables
    ! note
    !   Variables read directly from TOML via get_value → character(len=:), allocatable
    !   Everything else (paths built by concatenation, command-line args, hardcoded strings) → character(len=512)

    ! read in from distributor.f90
    character(len=512) :: point_numb, settings_dir
    integer :: log_unit = 6   ! 6 = stdout; set to a file unit by distributor worker before each Run_Model call
    
    ! read in from toml files
    character(len=:), allocatable :: project_name
    character(len=:), allocatable :: model_version
    character(len=:), allocatable :: code_dir, data_dir

    ! currently read in from toml, will be set using netcdf
    character(len=:), allocatable :: start_ts_year, end_ts_year, start_ave_year, end_ave_year, model_first_timestep, model_last_timestep

    ! created in here in model_settings
    character(len=512) :: reference_dir
    character(len=512) :: input_dir, input_averages_dir, input_timeseries_dir
    character(len=512) :: project_dir, output_dir, restart_dir
    character(len=512) :: fname_mask
    character(len=512) :: prefix_forcing_timeseries, suffix_forcing_timeseries, suffix_forcing_averages
    character(len=512) :: fname_restart_from_spinup, prefix_fname_run, suffix_fname_run, fname_restart_from_previous_run
    character(len=512) :: fname_out_1d, fname_out_2d, fname_out_2ddet
    character(len=512) :: iceshelf_var
    double precision :: pi, seconds_per_year

    ! structures
    type(Settings), allocatable :: config
    type(Constants), allocatable :: const

contains

!TKTKTK: why have this be separate from subroutine Load_Model_Settings()
subroutine Read_Settings(table, settings_out)
    !*** Define model settings and physics from model.toml ***!

    type(toml_table), allocatable, intent(inout) :: table
    type(toml_table), pointer :: child
    type(Settings), intent(out), allocatable :: settings_out

    allocate(settings_out)

    ! Model physics
    nullify(child)
    call get_value(table, 'model_physics', child, requested=.false.)
    if (associated(child)) then
        call get_value(child, 'do_mo_fit', settings_out%model_physics%do_MO_fit)
    end if

    ! Forcing settings --> TKTKTK: to be read in using netcdf
    nullify(child)
    call get_value(table, 'forcing_dimensions', child, requested=.false.)
    if (associated(child)) then
        call get_value(child, 'nyears_spinup', settings_out%forcing_dimensions%nyears_spinup)
        call get_value(child, 'nyears_forcing', settings_out%forcing_dimensions%nyears_forcing)
        call get_value(child, 'dtobs', settings_out%forcing_dimensions%dtobs)
        call get_value(child, 'Nlon', settings_out%forcing_dimensions%Nlon)
        call get_value(child, 'Nlat', settings_out%forcing_dimensions%Nlat)
        call get_value(child, 'Nlon_timeseries', settings_out%forcing_dimensions%Nlon_timeseries)
        call get_value(child, 'Nt_forcing', settings_out%forcing_dimensions%Nt_forcing)
    else
        write(std_err, *) "Cannot find section >forcing_dimensions< in model.toml file."
        stop
    end if

    ! Initialization
    nullify(child)
    call get_value(table, 'initialization', child, requested=.false.)
    if (associated(child)) then
        call get_value(child, 'startasice', settings_out%initialization%startasice)
        call get_value(child, 'beginT', settings_out%initialization%beginT)
        call get_value(child, 'initdepth', settings_out%initialization%initdepth)
    else
        write(std_err, *) "Cannot find section >initialization< in model.toml file."
        stop
    end if

    ! Model timesteps and dimensions
    nullify(child)
    call get_value(table, 'model_choices', child, requested=.false.)
    if (associated(child)) then
        call get_value(child, 'ImpExp', settings_out%model_choices%ImpExp)
        call get_value(child, 'dtmodelImp', settings_out%model_choices%dtmodelImp)
        call get_value(child, 'dtmodelExp', settings_out%model_choices%dtmodelExp)
        call get_value(child, 'dtSnow', settings_out%model_choices%dtSnow)
        call get_value(child, 'dzmax', settings_out%model_choices%dzmax)
        call get_value(child, 'th', settings_out%model_choices%th)
        call get_value(child, 'ts_minimum',    settings_out%model_choices%ts_minimum)
        call get_value(child, 'det2d_minimum', settings_out%model_choices%det2d_minimum)
    else
        write(std_err, *) "Cannot find section >model_choices< in model.toml file."
        stop
    end if

    ! Output dimensions
    nullify(child)
    call get_value(table, 'output_dimensions', child, requested=.false.)
    if (associated(child)) then
        call get_value(child, 'writeinspeed', settings_out%output_dimensions%writeinspeed)
        call get_value(child, 'writeinprof', settings_out%output_dimensions%writeinprof)
        call get_value(child, 'proflayers', settings_out%output_dimensions%proflayers)
        call get_value(child, 'writeindetail', settings_out%output_dimensions%writeindetail)
        call get_value(child, 'detlayers', settings_out%output_dimensions%detlayers)
        call get_value(child, 'detthick', settings_out%output_dimensions%detthick)
    else
        write(std_err, *) "Cannot find section >output_dimensions< in model.toml file."
        stop
    end if

    ind_z_surf = nint(initdepth/dzmax)  ! initial amount of vertical layers

    write(log_unit, *) "Loaded model settings"
    write(log_unit, *) " "

end subroutine Read_Settings


subroutine Load_Model_Settings()
    !*** Read runtime model settings from TOML and expose commonly used values ***!

    integer :: fu, rc
    logical :: file_exists
    type(toml_error), allocatable :: parse_error
    type(toml_table), allocatable :: table

    settings_path_model = trim(settings_dir)//"model.toml"
    call Load_TOML(settings_path_model, table)

    if (allocated(config)) then
        deallocate(config)
    end if

    call read_settings(table, config) 

    deallocate(table)

end subroutine Load_Model_Settings

subroutine Read_Job()
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
    else
        write(std_err, *) "Cannot find section >job< in run.toml file."
        stop
    end if

    if (allocated(config)) then
        if (allocated(config%general_settings%domain)) domain = config%run_settings%domain
        if (allocated(config%general_settings%forcing)) forcing = config%run_settings%forcing
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
    else
        write(std_err, *) "Cannot find section >metadata< in run.toml file."
        stop
    end if

    call get_value(table, 'directories', child, requested=.false.)

    if (associated(child)) then
        call get_value(child, 'code_dir', code_dir)
        call get_value(child, 'data_dir', data_dir)
    else
        write(std_err, *) "Cannot find section >directories< in run.toml file."
        stop
    end if

    code_dir      = trim(adjustl(code_dir))
    data_dir      = trim(adjustl(data_dir))
    project_dir   = trim(data_dir)//trim(project_name)//"/"
    output_dir    = trim(project_dir)//"output/"
    restart_dir   = trim(project_dir)//"restart/"
    reference_dir = trim(code_dir)//"reference/"//trim(config%run_settings%domain)//"/"
    input_dir     = trim(data_dir)//trim(config%run_settings%domain)//"_"//trim(config%run_settings%forcing)//"/input/"

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

end subroutine Read_Job

subroutine Read_Constants(table, const)

    type(toml_table), allocatable, intent(inout) :: table
    type(toml_table), pointer :: child
    type(Constants), intent(out), allocatable :: const

    call get_value(table, "constants", child, requested=.false.)
    if (.not. associated(child)) then
        write(stderr, '("Cannot find section >constants< in constants.toml file.")')
        stop
    end if

    allocate(const)

    call get_value(child, "R", const%R)
    call get_value(child, "g", const%g)
    call get_value(child, "rhoi", const%rhoi)
    call get_value(child, "rho_ocean", const%rho_ocean)
    call get_value(child, "Tmelt", const%Tmelt)
    call get_value(child, "Ec", const%Ec)
    call get_value(child, "Eg", const%Eg)
    call get_value(child, "Lh", const%Lh)
    call get_value(child, "days_per_year", const%days_per_year)
    call get_value(child, "NaN_value", const%NaN_value)

    const%pi = asin(1.)*2 ! pi = 3.1415
    const%seconds_per_year = 3600.*24.*const%days_per_year   ! seconds per year [s]

end subroutine Read_Constants

! *******************************************************

subroutine Load_Constants()
    !*** Read runtime model settings from TOML and expose commonly used values ***!

    integer :: fu, rc
    logical :: file_exists
    type(toml_error), allocatable :: parse_error
    type(toml_table), allocatable :: table

    settings_path_constants = trim(settings_dir)//"constants.toml"
    call Load_TOML(settings_path_constants, table)

    if (allocated(const)) then
        deallocate(const)
    end if

    call read_settings(table, const) 

    deallocate(table)

end subroutine Load_Model_Settings

! *******************************************************


subroutine Define_Filenames()
    !*** Build all output, restart, and input filenames ***!

    character(len=512) :: domain_forcing

    !TKTKTK: is this the correct way of using the struct variables?
    domain_forcing = trim(config%run_settings%domain)//"_"//trim(config%run_settings%forcing)

    fname_mask                = trim(config%run_settings%domain)//"_Masks.nc"
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

    if (config%run_settings%domain == "ANT27") then
        iceshelf_var = "IceShelve"
    else
        write(log_unit, *) "No iceshelf_var defined for domain, ", config%run_settings%domain
        write(log_unit, *) "so not loading ice shelf mask."
    end if


end subroutine Define_Filenames


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
