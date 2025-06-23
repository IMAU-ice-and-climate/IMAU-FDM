module model_settings
    !*** Subroutines for setting paths and constants for global use throughout the model ***!

    implicit none

    private

    public :: Define_Paths, Define_Constants, Define_Settings

    ! Module-level variables that can be accessed by other modules
    public :: domain
    public :: path_settings, path_forcing_dims, path_forcing_mask, path_forcing_averages
    public :: path_forcing_timeseries, path_restart, path_out_1d, path_out_2d, path_out_2ddet
    public :: fname_settings, fname_forcing_dims, fname_mask, suffix_forcing_averages
    public :: prefix_forcing_timeseries, fname_restart_from_spinup, prefix_fname_ini
    public :: suffix_fname_ini, fname_restart_from_previous_run, fname_out_1d
    public :: fname_out_2d, fname_out_2ddet, iceshelf_var
    public :: rhoi, rho_ocean, Tmelt, NaN_value, R, pi, Ec, Eg, g, Lh, seconds_per_year, ts_minimum, det2d_minimum

    ! Declare the module variables
    character(len=255) :: domain
    character(len=255) :: path_settings, path_forcing_dims, path_forcing_mask, path_forcing_averages
    character(len=255) :: path_forcing_timeseries, path_restart, path_out_1d, path_out_2d, path_out_2ddet
    character(len=255) :: fname_settings, fname_forcing_dims, fname_mask, suffix_forcing_averages
    character(len=255) :: prefix_forcing_timeseries, fname_restart_from_spinup, prefix_fname_ini
    character(len=255) :: suffix_fname_ini, fname_restart_from_previous_run, fname_out_1d
    character(len=255) :: fname_out_2d, fname_out_2ddet, iceshelf_var
    double precision :: rhoi, rho_ocean, Tmelt, NaN_value, R, pi, Ec, Eg, g, Lh, seconds_per_year, ts_minimum, det2d_minimum

contains

subroutine Define_Paths(username, prefix_output, point_numb, project_name, domain_name)
    !*** Definition of all paths used to read in or write out files ***!

    ! declare arguments
    character*255, intent(in) :: username, prefix_output, point_numb, project_name, domain_name

    ! define local variables
    character*255 :: forcing, code_dir, output_dir, input_dir, start_ts_year, end_ts_year, start_ave_year, end_ave_year

    domain = domain_name
    forcing = "era055"
    start_ts_year = "1939"
    end_ts_year = "2023"
    start_ave_year = "1940"
    end_ave_year = "1970"

    code_dir = "/perm/"//trim(username)//"/code/IMAU-FDM/"
    output_dir = "/ec/res4/scratch/"//trim(username)//"/"//trim(project_name)//"/"
    input_dir = "/ec/res4/scratch/"//trim(username)//"/"//trim(domain)//"_"//trim(forcing)//"/input/"
    path_restart = "/ec/res4/scratch/"//trim(username)//"/restart/"//trim(project_name)//"/"

    print *, "Path to IMAU-FDM: ", code_dir
    print *, "Path to output directory: ", output_dir
    print *, "Path to input directory: ", input_dir
    print *, "Path to restart directory: ", path_restart

    path_settings = trim(output_dir)//"ms_files/"
    path_forcing_dims =  trim(code_dir)//"reference/FGRN055/"
    path_forcing_mask = trim(code_dir)//"reference/FGRN055/"
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
    prefix_forcing_timeseries = "_"//trim(prefix_output)//"_"//trim(start_ts_year)//"-"//trim(end_ts_year)//"_p"
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

end subroutine Define_Paths



! ******************************************************* 



subroutine Define_Constants()
    !*** Define physical constants ***!
    
    pi = asin(1.)*2         ! pi = 3.1415
    R = 8.3145              ! gas constant [J mole-1 K-1]
    g = 9.81                ! gravitational acceleration [m s-2]
    rhoi = 917.             ! density of ice [kg m-3]
    rho_ocean = 1027.       ! density of ocean water [kg m-3]
    Tmelt = 273.15          ! melting point of ice [K]
    Ec = 60000.             ! activation energy for creep [J mole-1]
    Eg = 42400.             ! activation energy for grain growth [J mole-1]
    Lh = 333500.            ! latent heat of fusion [J kg-1]

    seconds_per_year = 3600.*24.*365.   ! seconds per year [s]
    NaN_value = 9.96921e+36 ! missing value for doubles as used in the NCL scripts
    
!    ts_minimum = 1.e-04     ! minimum magnitude for timeseries value, set when ts are loaded
!    det2d_minimum = 1.e-05 ! minimum magnitude for refreezing sum in 2ddetail output

    ! kg = 1.3E-7             ! rate constant for grain growth [m2 s-1]
    ! Ec2 = 49000.            ! activation energy grain boundary diffusion [J mole-1]
    ! rgrain_refreeze = 1.E-3    ! grain size refrozen ice [m]

end subroutine Define_Constants



! ******************************************************* 



subroutine Define_Settings
    !*** Under construction ***!
    !*** Define model settings and physics ***!

    ! Model settings
    ts_minimum = 1.e-04     ! minimum magnitude for timeseries value, set when ts are loaded
    det2d_minimum = 1.e-05 ! minimum magnitude for refreezing sum in 2ddetail output

    ! Model physics

end subroutine Define_Settings

end module model_settings