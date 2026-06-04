module model_settings
    !*** Subroutines for setting paths and constants for global use throughout the model ***!

    use, intrinsic :: iso_fortran_env, only: stderr => error_unit
    use tomlf
    use netcdf, only: nf90_open, nf90_close, nf90_inq_dimid, nf90_inquire_dimension, &
        nf90_inq_varid, nf90_inquire, nf90_inquire_variable, nf90_get_var, &
        nf90_get_att, nf90_global, &
        nf90_noerr, nf90_strerror, nf90_max_name
    implicit none

    private ! none?

    integer, parameter, public :: ind_z_max = 20000
    
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
        real(kind=4) :: NaN_value  ! float32 fill/missing value, matching the nf90_real output vars (single so the missing_value attr type matches the variable; was double -> netCDF "cannot safely cast" warning)
        double precision :: pi
        double precision :: seconds_per_year
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

    ! Single source for run/forcing metadata. The date fields are read from the
    ! forcing netCDF files (timeseries date_bnds, averages global attrs) rather
    ! than from TOML; model_version still comes from run.toml [metadata].
    type, public :: metadata_t
        character(5)  :: start_ts_year
        character(5)  :: end_ts_year
        character(5)  :: start_ave_year
        character(5)  :: end_ave_year
        character(30) :: model_first_timestep
        character(30) :: model_last_timestep
        character(30) :: model_version
        integer :: num_time_steps
        integer :: num_full_years
        integer :: seconds_per_timestep
    end type metadata_t

    type, public :: initialization_t

        integer :: startasice
        integer :: beginT
        double precision :: initdepth

    end type initialization_t

    ! Numerical / discretization settings (scheme, timesteps, resolution, thresholds).
    ! Physics-parameterisation choices live in model_physics_t.
    type, public :: model_choices_t

        integer :: ImpExp
        integer :: dtmodelImp
        integer :: dtmodelExp
        integer :: dtSnow
        double precision :: dzmax
        double precision :: th
        double precision :: ts_minimum
        double precision :: det2d_minimum

    end type model_choices_t

    ! Physics-parameterisation choices (densification refit, LWC scheme).
    ! Numerical/discretization settings live in model_choices_t.
    type, public :: model_physics_t

        logical :: do_MO_fit
        character(len=:), allocatable :: LWC_avail

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
        type(forcing_dimensions_t) :: forcing_dimensions
        type(initialization_t)     :: initialization
        type(model_choices_t)      :: model_choices   ! numerical/discretization settings: scheme (ImpExp), model timesteps, dzmax, theta, ts_minimum, det2d_minimum
        type(model_physics_t)      :: model_physics   ! physics-parameterisation choices: do_MO_fit, LWC_avail
        type(output_dimensions_t)  :: output_dimensions

    end type Settings

    public :: Read_Job, Read_Constants, Load_Constants, Read_Settings, Load_Model_Settings, Define_Filenames

    ! Module-level variables that can be accessed by other modules
    
    ! set in distributor.f90
    public :: point_numb
    public :: log_unit

    ! subfolder of code_dir; defined in launch_job file
    public :: settings_dir

    ! set in run.tmol
    public :: project_name, domain, forcing, restart_type
    public :: code_dir, forcing_root, output_root

    ! single metadata struct: years/timesteps read from netCDF, model_version from run.toml
    public :: metadata

    ! created in model_settings
    public :: reference_dir
    public :: forcing_dir, forcing_averages_dir, forcing_timeseries_dir
    public :: project_dir, restart_in_dir, restart_out_dir, output_dir
    public :: fname_mask
    public :: prefix_forcing_timeseries, suffix_forcing_timeseries, suffix_forcing_averages
    public :: fname_restart_from_spinup, prefix_fname_run, suffix_fname_run, fname_restart_from_previous_run
    public :: fname_out_1d, fname_out_2d, fname_out_2ddet
    public :: iceshelf_var
    public :: ind_z_surf
    public :: config, const

    ! netCDF helpers (moved here from openNetCDF so dimensions can be read while building settings)
    public :: Handle_Error, read_forcing_metadata, read_averages_metadata, Set_Forcing_Dimensions

    ! Declare the module variables
    ! note
    !   Variables read directly from TOML via get_value → character(len=:), allocatable
    !   Everything else (paths built by concatenation, command-line args, hardcoded strings) → character(len=512)

    ! read in from distributor.f90
    character(len=512) :: point_numb, settings_dir
    integer :: log_unit = 6   ! 6 = stdout; set to a file unit by distributor worker before each Run_Model call
    
    ! read in from toml files
    character(len=:), allocatable :: project_name, domain, forcing, restart_type
    character(len=:), allocatable :: code_dir, forcing_root, output_root

    ! created in here in model_settings
    character(len=512) :: reference_dir
    character(len=512) :: forcing_dir, forcing_averages_dir, forcing_timeseries_dir
    character(len=512) :: project_dir, output_dir, restart_in_dir, restart_out_dir
    character(len=512) :: fname_mask
    character(len=512) :: prefix_forcing_timeseries, suffix_forcing_timeseries, suffix_forcing_averages
    character(len=512) :: fname_restart_from_spinup, prefix_fname_run, suffix_fname_run, fname_restart_from_previous_run
    character(len=512) :: fname_out_1d, fname_out_2d, fname_out_2ddet
    character(len=512) :: iceshelf_var
    integer :: ind_z_surf

    ! structures
    type(Settings), allocatable :: config
    type(Constants), allocatable :: const
    type(metadata_t) :: metadata

contains

subroutine Handle_Error(stat, msg)
    !*** error stop for netCDF ***!

    integer, intent(in) :: stat
    character(len=*), intent(in) :: msg

    if (stat /= nf90_noerr) then
        write(log_unit, *) 'netCDF error (', msg, '): ', nf90_strerror(stat)
        stop
    endif

end subroutine Handle_Error


subroutine read_forcing_metadata(netcdf_file)
    !*** Read time metadata + forcing dimensions from a timeseries netCDF file.   ***!
    !*** Fills metadata% start_ts_year, end_ts_year, model_first/last_timestep    ***!
    !***   and config%forcing_dimensions% Nft_forcing, Nlon_timeseries, dtobs,     ***!
    !***   nyears_forcing.                                                         ***!

    character(len=*), intent(in) :: netcdf_file

    integer :: ncid, bnds_varid, dim_id, idim, dtg_varid
    integer :: counts(2), bound_dimids(2)
    integer, allocatable :: bounds(:, :), dtg(:)
    integer :: n_lon_ts, dtobs
    integer :: start_date, end_date, start_year, end_year
    character(len=nf90_max_name) :: dim_name
    integer :: diff_days, diff_hours, num_years

    call Handle_Error(nf90_open(trim(netcdf_file), 0, ncid), "open timeseries metadata file")

    ! --- Time range: first/last verifying date from date_bnds (yyyymmdd) ---
    call Handle_Error(nf90_inq_varid(ncid, "date_bnds", bnds_varid), "inq date_bnds")
    call Handle_Error(nf90_inq_varid(ncid, "dtg", dtg_varid), "Reading var id of dtg")

    call Handle_Error(nf90_inquire_variable(ncid, bnds_varid, dimids=bound_dimids), "date_bnds dimids")
    do idim = 1, 2
        call Handle_Error(nf90_inquire_dimension(ncid, bound_dimids(idim), len=counts(idim)), "date_bnds counts")
    end do
    allocate(bounds(counts(1), counts(2)))
    allocate(dtg(counts(2)))

    call Handle_Error(nf90_get_var(ncid, bnds_varid, bounds), "get date_bnds")
    call Handle_Error(nf90_get_var(ncid, dtg_varid, dtg), "Reading starting dates")

    if (.not. dtg(2)/10000 == dtg(1)/10000) then
        print*, "Error getting seconds per year: timesteps of more than one month are unsupported."
        stop
    end if

    ! dtobs can be calculated in fortran 
    diff_days = dtg(2)/100 - dtg(1)/100
    diff_hours = mod(dtg(2), 100) - mod(dtg(1), 100)
    metadata%seconds_per_timestep = 3600*diff_hours + 24*3600*diff_days
 
    metadata%num_time_steps = counts(2)
    start_date = bounds(1, 1)
    end_date   = bounds(1, counts(2))
    start_year = start_date / 10000
    end_year   = end_date / 10000

    num_years = bounds(2, counts(2)) / 10000 - bounds(1, 1) / 10000 - 1
    if (mod(bounds(1, 1), 10000) == 101) then
        num_years = num_years + 1
    end if

    metadata%num_full_years = num_years

     write(metadata%num_full_years, "(i2)") num_years
     write(metadata%num_full_years, "(i2)") num_years

    write(metadata%start_ts_year, "(i4)") start_year
    write(metadata%end_ts_year, "(i4)") end_year
    write(metadata%model_first_timestep, "(i4,'-',i2.2,'-',i2.2,'T00:00:00')") &
        start_year, mod(start_date/100, 100), mod(start_date, 100)
    write(metadata%model_last_timestep, "(i4,'-',i2.2,'-',i2.2,'T21:00:00')") &
        end_year, mod(end_date/100, 100), mod(end_date, 100)

    ! Nt_forcing = number of timesteps = length of the time (non-bnds) dimension
    config%forcing_dimensions%Nt_forcing = counts(2)
    ! nyears_forcing = inclusive span of forcing years
    config%forcing_dimensions%nyears_forcing = (end_year - start_year) + 1

    ! --- Nlon_timeseries: rlon dimension of this parts file ---
    call Handle_Error(nf90_inq_dimid(ncid, "rlon", dim_id), "inq rlon (timeseries)")
    call Handle_Error(nf90_inquire_dimension(ncid, dim_id, dim_name, n_lon_ts), "get rlon len")
    config%forcing_dimensions%Nlon_timeseries = n_lon_ts

    ! --- dtobs: forcing timestep [s], stamped as a global attr by preprocessing ---
    call Handle_Error(nf90_get_att(ncid, nf90_global, "dtobs", dtobs), "get dtobs attr")
    config%forcing_dimensions%dtobs = dtobs

    call Handle_Error(nf90_close(ncid), "close timeseries metadata file")

end subroutine read_forcing_metadata


subroutine read_averages_metadata(netcdf_file)
    !*** Read grid dimensions (rlon->Nlon, rlat->Nlat) and the spin-up averaging  ***!
    !*** period (start_ave_year/end_ave_year global attrs) from an averages file.  ***!

    character(len=*), intent(in) :: netcdf_file

    integer :: ncid, dim_id
    integer :: n_lon, n_lat, start_ave, end_ave
    character(len=nf90_max_name) :: dim_name

    call Handle_Error(nf90_open(trim(netcdf_file), 0, ncid), "open averages metadata file")

    call Handle_Error(nf90_inq_dimid(ncid, "rlat", dim_id), "inq rlat")
    call Handle_Error(nf90_inquire_dimension(ncid, dim_id, dim_name, n_lat), "get rlat len")
    call Handle_Error(nf90_inq_dimid(ncid, "rlon", dim_id), "inq rlon (averages)")
    call Handle_Error(nf90_inquire_dimension(ncid, dim_id, dim_name, n_lon), "get rlon len")
    config%forcing_dimensions%Nlon = n_lon
    config%forcing_dimensions%Nlat = n_lat

    ! spin-up averaging period: global attrs stamped by preprocessing
    call Handle_Error(nf90_get_att(ncid, nf90_global, "start_ave_year", start_ave), "get start_ave_year attr")
    call Handle_Error(nf90_get_att(ncid, nf90_global, "end_ave_year", end_ave), "get end_ave_year attr")
    write(metadata%start_ave_year, "(i4)") start_ave
    write(metadata%end_ave_year, "(i4)") end_ave

    ! Spin-up length is definitional: the climatology is averaged over the INCLUSIVE
    ! window [start_ave_year, end_ave_year] (see pre-process-RACMO/make_fdm_averages.py,
    ! which loops range(start, end+1)), so nyears_spinup = end_ave_year - start_ave_year + 1.
    ! e.g. 1940..1969 -> 30 years. No longer read from model.toml.
    config%forcing_dimensions%nyears_spinup = end_ave - start_ave + 1

    call Handle_Error(nf90_close(ncid), "close averages metadata file")

end subroutine read_averages_metadata


subroutine Set_Forcing_Dimensions()
    !*** Populate config%forcing_dimensions + metadata from the forcing netCDF files. ***!
    !*** Replaces the hardcoded forcing dimensions/years in model.toml + run.toml.     ***!
    !*** Must be called after Define_Filenames() (needs the forcing filename parts).    ***!

    character(len=512) :: ts_file, avg_file

    ! Representative timeseries file to read metadata from. For real runs this is
    ! tskin part 1, which exists for every domain (unlike higher part numbers --
    ! e.g. ANT27 has only 17 lon-bands); the example run has a single point file
    ! with no part number.
    if (trim(project_name) == "example") then
        ts_file = trim(forcing_timeseries_dir)//"tskin"//trim(prefix_forcing_timeseries)//trim(suffix_forcing_timeseries)
    else
        ts_file = trim(forcing_timeseries_dir)//"tskin"//trim(prefix_forcing_timeseries)//"1"//trim(suffix_forcing_timeseries)
    end if
    avg_file = trim(forcing_averages_dir)//"tskin"//trim(suffix_forcing_averages)

    call read_forcing_metadata(trim(ts_file))
    call read_averages_metadata(trim(avg_file))

    ! Restart-from-previous-run file is named with the last forcing year, which is
    ! only known once the timeseries metadata has been read.
    fname_restart_from_previous_run = trim(prefix_fname_run)//trim(metadata%end_ts_year)//trim(suffix_fname_run)

    write(log_unit, *) "Forcing metadata read from:"
    write(log_unit, *) "  ", trim(ts_file)
    write(log_unit, *) "  ", trim(avg_file)
    write(log_unit, *) "  Nlon, Nlat       = ", config%forcing_dimensions%Nlon, config%forcing_dimensions%Nlat
    write(log_unit, *) "  Nt_forcing       = ", config%forcing_dimensions%Nt_forcing
    write(log_unit, *) "  Nlon_timeseries  = ", config%forcing_dimensions%Nlon_timeseries
    write(log_unit, *) "  dtobs            = ", config%forcing_dimensions%dtobs
    write(log_unit, *) "  nyears_forcing   = ", config%forcing_dimensions%nyears_forcing
    write(log_unit, *) "  ts years         = ", trim(metadata%start_ts_year), " - ", trim(metadata%end_ts_year)
    write(log_unit, *) "  ave years        = ", trim(metadata%start_ave_year), " - ", trim(metadata%end_ave_year)
    write(log_unit, *) " "

end subroutine Set_Forcing_Dimensions


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
        call get_value(child, 'LWC_avail', settings_out%model_physics%LWC_avail)

        ! stops if an option gets passed that doesn't match the options in water_physics
        if (settings_out%model_physics%LWC_avail /= "Coleou1998_corr" .and. &
            settings_out%model_physics%LWC_avail /= "Coleou1998_1p2") then
            write(stderr, *) "Invalid LWC_avail: '"//trim(settings_out%model_physics%LWC_avail)// &
                "'. Set a valid LWC_avail option (Coleou1998_corr | Coleou1998_1p2) in model.toml."
            stop
        endif
    
    end if

    ! All forcing dimensions + nyears_spinup are read from the forcing netCDF files in
    ! Set_Forcing_Dimensions() (nyears_spinup = end_ave_year - start_ave_year), so
    ! nothing forcing-related is read from model.toml anymore.

    ! Initialization
    nullify(child)
    call get_value(table, 'initialization', child, requested=.false.)
    if (associated(child)) then
        call get_value(child, 'startasice', settings_out%initialization%startasice)
        call get_value(child, 'beginT', settings_out%initialization%beginT)
        call get_value(child, 'initdepth', settings_out%initialization%initdepth)
    else
        write(stderr, *) "Cannot find section >initialization< in model.toml file."
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
        write(stderr, *) "Cannot find section >model_choices< in model.toml file."
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
        write(stderr, *) "Cannot find section >output_dimensions< in model.toml file."
        stop
    end if

    ind_z_surf = nint(settings_out%initialization%initdepth/settings_out%model_choices%dzmax)  ! initial amount of vertical layers

    write(log_unit, *) "Loaded model settings"
    write(log_unit, *) " "

end subroutine Read_Settings


subroutine Load_Model_Settings()
    !*** Read runtime model settings from TOML and expose commonly used values ***!

    character(len=512)            :: settings_path_model
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
    ! TKTKTK: check on domain, forcing, restart_type

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
        write(stderr, *) "Cannot find section >job< in run.toml file."
        stop
    end if

    call get_value(table, 'metadata', child, requested=.false.)
    if (associated(child)) then
        block
            character(len=:), allocatable :: version_str
            call get_value(child, 'model_version', version_str)
            metadata%model_version = version_str
        end block
        ! start/end_ts_year + model_first/last_timestep come from the timeseries netCDF,
        ! start/end_ave_year from the averages netCDF (both in Set_Forcing_Dimensions).
    else
        write(stderr, *) "Cannot find section >metadata< in run.toml file."
        stop
    end if

    call get_value(table, 'directories', child, requested=.false.)

    if (associated(child)) then
        call get_value(child, 'code_dir', code_dir)
        call get_value(child, 'forcing_root', forcing_root)    ! base for forcing input data
        call get_value(child, 'output_root', output_root)      ! base for model output
    else
        write(stderr, *) "Cannot find section >directories< in run.toml file."
        stop
    end if

    code_dir        = trim(adjustl(code_dir))
    forcing_root  = trim(adjustl(forcing_root))
    output_root = trim(adjustl(output_root))
    project_dir   = trim(output_root)//trim(project_name)//"/"
    output_dir    = trim(project_dir)//"output/"
    ! restart_out_dir: where THIS run WRITES its restart files (always its own project).
    ! restart_in_dir : where restart files are READ at startup (may point at another
    ! project via restart_project below).
    restart_out_dir = trim(project_dir)//"restart/"
    restart_in_dir  = trim(project_dir)//"restart/"
    reference_dir = trim(code_dir)//"reference/"//trim(domain)//"/"
    forcing_dir     = trim(forcing_root)//trim(domain)//"_"//trim(forcing)//"/input/"

    if (trim(project_name) == "example") then
        forcing_dir       = trim(forcing_root)//"input/"
        restart_out_dir = trim(code_dir)//"example/restart/"
        restart_in_dir  = trim(code_dir)//"example/restart/"
    end if

    ! Optional: READ initial restart files from a different project. Writes still go to
    ! this run's own restart_out_dir, so the source project is never modified.
    block
        character(len=:), allocatable :: restart_project
        nullify(child)
        call get_value(table, 'job', child, requested=.false.)
        if (associated(child)) then
            call get_value(child, 'restart_project', restart_project, default="")
            if (allocated(restart_project) .and. len_trim(restart_project) > 0) then
                restart_in_dir = trim(output_root)//trim(restart_project)//"/restart/"
                write(log_unit, *) "Reading restart files from project: ", trim(restart_project)
            end if
        end if
    end block

    forcing_averages_dir   = trim(forcing_dir)//"averages/"
    forcing_timeseries_dir = trim(forcing_dir)//"timeseries/"

    write(log_unit, *) "Path to IMAU-FDM: ", code_dir
    write(log_unit, *) "Path to project directory: ", output_dir
    write(log_unit, *) "Path to input directory: ", forcing_dir
    write(log_unit, *) "Restart read dir:  ", restart_in_dir
    write(log_unit, *) "Restart write dir: ", restart_out_dir

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
    !*** Load constants from constnats.toml ***!
    ! TKTKTK: is this needed?

    character(len=512)            :: settings_path_constants
    integer :: fu, rc
    logical :: file_exists
    type(toml_error), allocatable :: parse_error
    type(toml_table), allocatable :: table

    settings_path_constants = trim(settings_dir)//"constants.toml"
    call Load_TOML(settings_path_constants, table)

    if (allocated(const)) then
        deallocate(const)
    end if

    call Read_Constants(table, const) 

    deallocate(table)

end subroutine Load_Constants

! *******************************************************


subroutine Define_Filenames()
    !*** Build all output, restart, and input filenames.                          ***!
    !*** Forcing filenames are date-free; the forcing's dates live in the netCDF   ***!
    !*** metadata (read in Set_Forcing_Dimensions), not in the filename.           ***!

    character(len=512) :: domain_forcing

    domain_forcing = trim(domain)//"_"//trim(forcing)

    fname_mask                = trim(domain)//"_Masks.nc"
    prefix_forcing_timeseries = "_"//trim(domain_forcing)//"_p"
    suffix_forcing_timeseries = ".nc"
    suffix_forcing_averages   = "_"//trim(domain_forcing)//"_ave.nc"

    fname_out_1d    = trim(domain_forcing)//"_1D_"//trim(point_numb)//".nc"
    fname_out_2d    = trim(domain_forcing)//"_2D_"//trim(point_numb)//".nc"
    fname_out_2ddet = trim(domain_forcing)//"_2Ddetail_"//trim(point_numb)//".nc"

    fname_restart_from_spinup = trim(domain_forcing)//"_restart_from_spinup_"//trim(point_numb)//".nc"
    prefix_fname_run          = trim(domain_forcing)//"_initialize_from_"
    suffix_fname_run          = "_run_"//trim(point_numb)//".nc"
    ! fname_restart_from_previous_run is built in Set_Forcing_Dimensions: it needs
    ! the last forcing year (metadata%end_ts_year), read from the timeseries netCDF.

    if (trim(project_name) == "example") then
        prefix_forcing_timeseries = "_"//trim(domain_forcing)
        suffix_forcing_timeseries = "_point_"//trim(point_numb)//".nc"
    end if

    if (trim(domain) == "ANT27") then
        iceshelf_var = "IceShelve"
    else
        write(log_unit, *) "No iceshelf_var defined for domain, ", trim(domain)
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
