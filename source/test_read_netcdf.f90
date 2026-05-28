program toml_main
    use netcdf, only: nf90_open, nf90_close, nf90_noerr, nf90_nowrite, nf90_max_name, nf90_inq_dimid, nf90_inquire_dimension, &
        nf90_inquire, nf90_inquire_attribute, nf90_inquire_variable, nf90_inq_varid, nf90_get_var
    use openNetCDF, only: read_avg_dimensions, handle_error, read_time_bounds, TimeBounds

    implicit none

    integer :: n_lon, n_lat
    integer :: status, ncid, ndim
    integer :: nvar, nattr, unlim_dim, format_num
    integer :: idim, dim_size
    integer :: iattr
    integer :: ivar
    integer :: bnds_varid

    integer, allocatable :: bounds(:, :)
    integer :: counts(2)
    integer :: bound_dimids(2)
    integer :: start_ts_year, end_ts_year, start_date, end_date

    character(len=nf90_max_name) :: dim_name, var_name
    character(len=50) :: model_first_timestep, model_last_timestep

    type(TimeBounds) :: time_bounds

    call read_avg_dimensions("data/evap_FGRN055_era055-1939_1940-1970_ave.nc", n_lon, n_lat)
    time_bounds = read_time_bounds("data/evap_FGRN055_era055_1939-2023_p24-007.nc")
    print*, time_bounds
    ! print *, n_lat, n_lon

    ! call handle_error(nf90_open("data/evap_FGRN055_era055_1939-2023_p24-007.nc", 0, ncid), "open_netcdf_file")
    ! call handle_error(nf90_inquire(ncid, ndim, nvar, nattr, unlim_dim, format_num), "inquire_net_file")

    ! ! print *, ndim, nvar, nattr, unlim_dim, format_num

    ! call Handle_Error(nf90_inq_varid(ncid, "date_bnds", bnds_varid), "Reading var id day_bnds")
    ! print *, bnds_varid

    ! call Handle_Error(nf90_inquire_variable(ncid, bnds_varid, dimids=bound_dimids), "Getting dim_ids for bounds")
    ! do idim = 1, 2
    !     call Handle_Error(nf90_inquire_dimension(ncid, bound_dimids(idim), len=counts(idim)), "Reading time bound counts")
    ! end do
    ! allocate(bounds(counts(1), counts(2)))
    ! print*, counts
    ! ! call Handle_Error
    ! call Handle_Error(nf90_get_var(ncid, bnds_varid, bounds), "Reading time bounds")

    ! start_date = bounds(1, 1)
    ! end_date = bounds(1, counts(2))
    ! start_ts_year = start_date / 10000
    ! end_ts_year = bounds(1, counts(2)) / 10000
    ! print*, bounds(1, 1), bounds(1, counts(2))
    ! print*, end_ts_year, start_ts_year
    ! write(model_first_timestep, "(i4, '-', i2.2, '-', i2.2, 'T00:00:00')") start_ts_year, mod(start_date/100, 100), mod(start_date, 100)
    ! write(model_last_timestep, "(i4, '-', i2.2, '-', i2.2, 'T21:00:00')") end_ts_year, mod(end_date/100, 100), mod(end_date, 100)
    ! print*, model_first_timestep, model_last_timestep
    ! do idim = 1, ndim
        ! call Handle_Error(nf90_inquire_dimension(ncid, idim, dim_name, dim_size), "inquire_dimension")
        ! print *, idim, dim_name, dim_size
    ! end do
! 
    ! do ivar = 1, nvar
        ! call Handle_Error(nf90_inquire_variable(ncid, ivar, var_name), "Inquire variable")
        ! print *, var_name
    ! end do
    ! do iattr = 1, nattr
        ! call Handle_Error(nf90_inquire_attribute(ncid, ))
end program toml_main