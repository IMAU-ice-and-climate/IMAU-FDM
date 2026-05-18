module model_settings
    !*** Subroutines for setting paths and constants for global use throughout the model ***!

    use, intrinsic :: iso_fortran_env, only: stderr => error_unit
    use tomlf
    implicit none

    private
    

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
        double precision :: seconds_per_year  ! days per year TKTKTK: remove once timestamps incorprated
        double precision :: NaN_value  ! missing value for doubles as used in the NCL scripts
        double precision :: pi ! pi = 3.1415
    end type Constants


    type, public :: minimum_values_t
        double precision :: ts_minimum
        double precision :: det2d_minimum
    end type minimum_values_t

    type, public :: general_settings_t
        character(len=:), allocatable :: run_location
        character(len=:), allocatable :: restart_type
        character(len=:), allocatable :: domain
        character(len=:), allocatable :: forcing
        integer :: save_output
        integer :: nyears_spinup
        integer :: ImpExp
        integer :: dtmodelExp
        integer :: dtmodelImp
        integer :: dtSnow
        double precision :: dzmax
    end type general_settings_t

    type, public :: model_physics_t
        logical :: do_MO_fit
        double precision :: initdepth
        double precision :: th
        integer :: startasice
        integer :: beginT
        integer :: NoR
    end type model_physics_t

    type, public :: output_dimensions_t
        integer :: writeinspeed
        integer :: writeinprof
        integer :: proflayers
        integer :: writeindetail
        integer :: detlayers
        double precision :: detthick
    end type output_dimensions_t

    type, public :: Settings
        type(minimum_values_t)  :: minimum_values
        type(general_settings_t):: general_settings
        type(model_physics_t):: model_physics
        type(output_dimensions_t) :: output_dimensions
    end type Settings


    public :: Define_Paths, Define_Constants, read_settings, Get_All_Command_Line_Arg, read_constants, Load_Model_Settings

    ! Module-level variables that can be accessed by other modules
    
    public :: point_numb
    public :: domain, paths_file, model_settings_file, username, prefix_output, project_name, forcing, data_dir
    public :: path_settings, path_forcing_dims, path_forcing_mask, path_forcing_averages
    public :: path_forcing_timeseries, path_restart, path_out_1d, path_out_2d, path_out_2ddet
    public :: fname_settings, fname_forcing_dims, fname_mask, suffix_forcing_averages
    public :: prefix_forcing_timeseries, suffix_forcing_timeseries, fname_restart_from_spinup 
    public :: prefix_fname_ini, suffix_fname_ini, fname_restart_from_previous_run, fname_out_1d
    public :: fname_out_2d, fname_out_2ddet, iceshelf_var
    public :: seconds_per_year
    public :: config, const
    public :: model_first_timestep, model_last_timestep

    ! Declare the module variables
    ! note
    !   Variables read directly from TOML via get_value → character(len=:), allocatable
    !   Everything else (paths built by concatenation, command-line args, hardcoded strings) → character(len=512)

    character(len=512) :: domain, paths_file, point_numb, username, prefix_output, project_name, forcing, data_dir
    character(len=512) :: path_settings, path_forcing_dims, path_forcing_mask, path_forcing_averages
    character(len=512) :: path_forcing_timeseries, path_out_1d, path_out_2d, path_out_2ddet
    character(len=512) :: fname_settings, fname_forcing_dims, fname_mask, suffix_forcing_averages
    character(len=512) :: prefix_forcing_timeseries, suffix_forcing_timeseries, fname_restart_from_spinup
    character(len=512) :: prefix_fname_ini, suffix_fname_ini, fname_restart_from_previous_run, fname_out_1d
    character(len=512) :: fname_out_2d, fname_out_2ddet, iceshelf_var
    character(len=512) :: model_first_timestep, model_last_timestep
    double precision :: seconds_per_year
    character(len=:), allocatable :: path_restart, restart_type, model_settings_file
    type(Settings), allocatable :: config
    type(Constants), allocatable :: const
contains


subroutine read_settings(table, settings_out)
    !*** Define model settings and physics ***!

    ! Local variables toml tables

    type(toml_table), allocatable, intent(inout) :: table
    type(toml_table), pointer :: child
    type(Settings), intent(out), allocatable :: settings_out

    allocate(settings_out)

    ! Minimum_values
    nullify(child)
    call get_value(table, 'minimum_values', child, requested=.false.)
    if (associated(child)) then
        call get_value(child, 'ts_minimum',    settings_out%minimum_values%ts_minimum)
        call get_value(child, 'det2d_minimum', settings_out%minimum_values%det2d_minimum)
    end if

    ! General_settings
    nullify(child)
    call get_value(table, 'general_settings', child, requested=.false.)
    if (associated(child)) then
        call get_value(child, 'run_location', settings_out%general_settings%run_location)
        call get_value(child, 'restart_type', settings_out%general_settings%restart_type)
        call get_value(child, 'save_output', settings_out%general_settings%save_output)
        call get_value(child, 'nyears_spinup', settings_out%general_settings%nyears_spinup)
        call get_value(child, 'domain', settings_out%general_settings%domain)
        call get_value(child, 'forcing', settings_out%general_settings%forcing)
        call get_value(child, 'ImpExp', settings_out%general_settings%ImpExp)
        call get_value(child, 'dtmodelExp', settings_out%general_settings%dtmodelExp)
        call get_value(child, 'dtmodelImp', settings_out%general_settings%dtmodelImp)
        call get_value(child, 'dtSnow', settings_out%general_settings%dtSnow)
        call get_value(child, 'dzmax', settings_out%general_settings%dzmax)
    end if

    ! Model physics
    nullify(child)
    call get_value(table, 'model_physics', child, requested=.false.)
    if (associated(child)) then
        call get_value(child, 'DO_MO_FIT', settings_out%model_physics%do_MO_fit)
        call get_value(child, 'initdepth', settings_out%model_physics%initdepth)
        call get_value(child, 'th', settings_out%model_physics%th)
        call get_value(child, 'startasice', settings_out%model_physics%startasice)
        call get_value(child, 'beginT', settings_out%model_physics%beginT)
        call get_value(child, 'NoR', settings_out%model_physics%NoR)
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
    end if
end subroutine read_settings


subroutine Load_Model_Settings()
    !*** Read runtime model settings from TOML and expose commonly used values ***!

    integer :: fu, rc
    logical :: file_exists
    type(toml_error), allocatable :: parse_error
    type(toml_table), allocatable :: table

    inquire(file=model_settings_file, exist=file_exists)

    if (.not. file_exists) then
        write(stderr, '("Error: TOML file ", a, " not found")') model_settings_file
        stop
    end if

    open(action='read', file=model_settings_file, iostat=rc, newunit=fu)

    if (rc /= 0) then
        write(stderr, '("Error: Reading TOML file ", a, " failed")') model_settings_file
        stop
    end if

    call toml_parse(table, fu, parse_error)
    close(fu)

    if (.not. allocated(table)) then
        if (allocated(parse_error)) then
            write(stderr, '("Error: TOML parsing failed for ", a)') model_settings_file
            write(stderr, '(a)') parse_error%message
        else
        write(stderr, '("Error: Parsing failed")')
        end if
        stop
    end if

    if (allocated(config)) then
        deallocate(config)
    end if

    call read_settings(table, config)

    deallocate(table)

end subroutine Load_Model_Settings


subroutine read_constants(table, const)

    type(toml_table), allocatable, intent(inout) :: table
    type(toml_table), pointer :: child
    type(Constants), intent(out), allocatable :: const

    call get_value(table, "constants", child, requested=.false.)
    if (.not. associated(child)) then
        write(stderr, '("Cannot find section constants in toml file.")')
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


end subroutine

    
! *******************************************************


subroutine Get_All_Command_Line_Arg()
    !*** Get all command line arguments ***!

    ! 1: Paths toml file
    ! 2: Model settings toml file
    ! 3: Simulation number, corresponding to the line number in the IN-file.
    ! 6: Restart type, where none= do spinup; spinup = restart from spinup

    print *, "Getting command line arguments..."

    username="ndl4814"
    domain="FGRN055"
    prefix_output="FGRN055_era055"
    project_name="test-new-distributor"
    restart_type="none"
    
    call get_command_argument(1, paths_file)
    call get_command_argument(2, model_settings_file)
    call get_command_argument(3, point_numb)
    call get_command_argument(4, prefix_output)
    call get_command_argument(5, project_name)

    if (trim(project_name) == "example") then

        print *, "Running example point 1 from sample data"

    else
        print *, "Running user defined point"
        ! Reads in current point number, config files, prefix, and project_name for path setting
    
    end if
    

   
end subroutine Get_All_Command_Line_Arg


! *******************************************************


subroutine Define_Paths()
    !*** Definition of all paths used to read in or write out files ***!

    ! Local variables toml tables

    integer                       :: fu, rc, i
    logical                       :: file_exists
    type(toml_table), allocatable :: table
    type(toml_table), pointer     :: child

    ! define local variables
    character*255 :: start_ts_year, end_ts_year, start_ave_year, end_ave_year
    character(len=:), allocatable :: code_dir, input_dir, output_dir
    character(len=:), allocatable :: project_name_toml, username_toml, data_dir_toml, forcing_toml, domain_toml

    inquire (file=paths_file, exist=file_exists)
    
    if (.not. file_exists) then
        write (stderr, '("Error: TOML file ", a, " not found")') paths_file
        stop
    end if

    open (action='read', file=paths_file, iostat=rc, newunit=fu)

    if (rc /= 0) then
        write (stderr, '("Error: Reading TOML file ", a, " failed")') paths_file
        stop
    end if

    call toml_parse(table, fu)
    close (fu)

    if (.not. allocated(table)) then
        write (stderr, '("Error: Parsing failed")')
        stop
    end if

    ! parameters ## to be moved to config files

    start_ts_year = "1939"
    end_ts_year = "2023"
    start_ave_year = "1940"
    end_ave_year = "1970"

    model_first_timestep = "1939-09-01T00:00:00"
    model_last_timestep = "2023-12-31T21:00:00"
    ! paths

    ! find section.
    
    call get_value(table, 'directories', child, requested=.false.)
    
    if (.not. associated(child)) then
        write(stderr, '("Cannot find section constants in toml file.")')
        stop
    end if

    call get_value(child, 'project_name', project_name_toml)
    project_name = trim(adjustl(project_name_toml))

    print *, " "
    print *, "Project name: ", project_name

 
    call get_value(child, 'username', username_toml)
    call get_value(child, 'code_dir', code_dir)
    call get_value(child, 'data_dir', data_dir_toml)
    call get_value(child, 'forcing', forcing_toml)
    call get_value(child, 'domain', domain_toml)

    username = trim(adjustl(username_toml))
    data_dir = trim(adjustl(data_dir_toml))
    forcing = trim(adjustl(forcing_toml))
    domain = trim(adjustl(domain_toml))

    print *, "Username: ", username
    print *, " "

    if (allocated(config)) then
        if (allocated(config%general_settings%domain)) domain = config%general_settings%domain
        if (allocated(config%general_settings%forcing)) forcing = config%general_settings%forcing
    end if

    ! integer :: i

    data_dir = trim(adjustl(data_dir))

    output_dir = trim(data_dir)//"output/"

    if (trim(project_name) == "example") then 

        input_dir = trim(data_dir)//"input/"
        path_restart =  trim(code_dir)//"example/restart/"
        prefix_forcing_timeseries = "_"//trim(prefix_output)//"_"//trim(start_ts_year)//"-"//trim(end_ts_year)
        suffix_forcing_timeseries = "_point_"//trim(point_numb)//".nc"

    else
        input_dir = trim(data_dir)//trim(domain)//"_"//trim(forcing)//"/input/"
        path_restart = "/ec/res4/scratch/"//trim(username)//"/restart/"//trim(project_name)//"/"
        prefix_forcing_timeseries = "_"//trim(prefix_output)//"_"//trim(start_ts_year)//"-"//trim(end_ts_year)//"_p"
        suffix_forcing_timeseries = ".nc"

    end if

    print *, "Path to IMAU-FDM: ", code_dir
    print *, "Path to output directory: ", output_dir
    print *, "Path to input directory: ", input_dir
    print *, "Path to restart directory: ", path_restart

    path_settings = trim(data_dir)//"ms_files/"
    path_forcing_dims =  trim(code_dir)//"reference/"//trim(domain)//"/"
    path_forcing_mask = trim(code_dir)//"reference/"//trim(domain)//"/"
    path_forcing_averages = trim(input_dir)//"averages/"
    path_forcing_timeseries = trim(input_dir)//"timeseries/"
    path_out_1d = trim(output_dir)//"output/"
    path_out_2d = trim(output_dir)//"output/"
    path_out_2ddet = trim(output_dir)//"output/"

    ! define filenames
    fname_settings = "model_settings_"//trim(domain)//"_"//trim(point_numb)//".txt"
    fname_forcing_dims = "input_settings_"//trim(domain)//".txt"
    fname_mask = trim(domain)//"_Masks.nc"
    suffix_forcing_averages = "_"//trim(prefix_output)//"-"//trim(start_ts_year)//"_"//trim(start_ave_year)//"-"//trim(end_ave_year)//"_ave.nc"
    fname_restart_from_spinup = trim(prefix_output)//"_restart_from_spinup_"//trim(point_numb)//".nc"
    prefix_fname_ini = trim(prefix_output)//"_initialize_from_"
    suffix_fname_ini = "_run_"//trim(point_numb)//".nc"
    fname_restart_from_previous_run = trim(prefix_fname_ini)//trim(end_ts_year)//trim(suffix_fname_ini)
    fname_out_1d = trim(prefix_output)//"_1D_"//trim(point_numb)//".nc"
    fname_out_2d = trim(prefix_output)//"_2D_"//trim(point_numb)//".nc"
    fname_out_2ddet = trim(prefix_output)//"_2Ddetail_"//trim(point_numb)//".nc"

    if (domain == "ANT27") then
        iceshelf_var = "IceShelve"
    else
        print *, "No iceshelf_var defined for domain, ", domain
        print *, "so not loading ice shelf mask."
    end if
    
    deallocate(table)

end subroutine Define_Paths



! ******************************************************* 



subroutine Define_Constants()
    !*** Define physical constants ***!
    
    ! pi = asin(1.)*2         ! pi = 3.1415
    ! R = 8.3145              ! gas constant [J mole-1 K-1]
    ! g = 9.81                ! gravitational acceleration [m s-2]
    ! rhoi = 917.             ! density of ice [kg m-3]
    ! rho_ocean = 1027.       ! density of ocean water [kg m-3]
    ! Tmelt = 273.15          ! melting point of ice [K]
    ! Ec = 60000.             ! activation energy for creep [J mole-1]
    ! Eg = 42400.             ! activation energy for grain growth [J mole-1]
    ! Lh = 333500.            ! latent heat of fusion [J kg-1]
    ! days_per_year = 365.25    ! days per year [days]
    ! seconds_per_year = 3600.*24.*days_per_year   ! seconds per year [s]
    ! NaN_value = 9.96921e+36 ! missing value for doubles as used in the NCL scripts

    ! kg = 1.3E-7             ! rate constant for grain growth [m2 s-1]
    ! Ec2 = 49000.            ! activation energy grain boundary diffusion [J mole-1]
    ! rgrain_refreeze = 1.E-3    ! grain size refrozen ice [m]

end subroutine Define_Constants



! ******************************************************* 





end module model_settings
