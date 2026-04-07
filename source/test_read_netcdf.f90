program toml_main
    use netcdf, only: nf90_open, nf90_close, nf90_noerr, nf90_nowrite, nf90_max_name, nf90_inq_dimid, nf90_inquire_dimension
    use openNetCDF, only: read_avg_dimensions


    implicit none

    ! integer :: ncid, dim_id, dim_len_lat, dim_len_lon, status
    ! character(len=nf90_max_name) :: dim_name
    integer :: n_lon, n_lat
    call read_avg_dimensions("./data/ff10m_FGRN055_era055-1939_1940-1970_ave.nc", n_lon, n_lat)
    print *, n_lat, n_lon
end program toml_main