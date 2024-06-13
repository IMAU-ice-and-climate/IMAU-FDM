module openNetCDF
    !*** subroutines for loading forcing data

    use netcdf, only : nf90_open, nf90_inq_varid, nf90_inq_dimid, nf90_inquire_dimension, & 
        nf90_close, nf90_get_var, nf90_noerr, nf90_strerror

    implicit none
    private

    public :: Load_Ave_Forcing, Load_TimeSeries_Forcing, Init_From_File, Handle_Error
    
contains


! *******************************************************


subroutine Load_Ave_Forcing(AveTsurf, AveAcc, AveWind, AveSubl, AveSnowDrif, AveMelt, LSM, Nlat, Nlon, Latitude, Longitude, &
    ISM, domain, path_forcing_mask, fname_mask, fname_shelfmask, path_forcing_averages, suffix_forcing_averages)
    !*** Load average forcing data from netCDF files ***!
    
    ! declare arguments
    integer, intent(in) :: Nlat, Nlon
    double precision, dimension(Nlon,Nlat), intent(out) :: LSM, ISM, Latitude, Longitude
    double precision, dimension(Nlon,Nlat), intent(inout) :: AveTsurf, AveAcc, AveWind, AveSubl, AveSnowDrif, AveMelt
    character*255, intent(in) :: domain, path_forcing_mask, fname_mask, fname_shelfmask, path_forcing_averages, suffix_forcing_averages

    ! declare local variables
    integer :: ind_lon, ind_lat
    integer :: start_lon, start_lat, start_z, start_t, Nlon_file, Nlat_file, Nz_file, Nt_file, status, ncid_mask, ncid_shelf, varid_lsm, varid_lat, varid_lon, varid_ism

    print *, 'trim(path_forcing_mask)//trim(fname_mask)'
    print *, trim(path_forcing_mask)//trim(fname_mask)
    print *, ' '

    ! Load land surface mask and coordinates
    status = nf90_open(trim(path_forcing_mask)//trim(fname_mask), 0, ncid_mask)
    if (status /= nf90_noerr) call Handle_Error(status, 'nf_open_mask')
    status = nf90_inq_varid(ncid_mask, "LSM", varid_lsm)
    if (status /= nf90_noerr) call Handle_Error(status, 'nf_inq_varid_lsm')    
    status = nf90_inq_varid(ncid_mask, "lat", varid_lat)
    if (status /= nf90_noerr) call Handle_Error(status, 'nf_inq_varid_lat')
    status = nf90_inq_varid(ncid_mask, "lon", varid_lon)
    if (status /= nf90_noerr) call Handle_Error(status, 'nf_inq_varid_lon')
    if (trim(domain) == 'Greenland') then
        status = nf90_get_var(ncid_mask, varid_lsm, LSM, start=(/1, 1/), count=(/Nlon, Nlat/))
    else
        status = nf90_get_var(ncid_mask, varid_lsm, LSM, start=(/1, 1, 1, 1/), count=(/Nlon, Nlat, 1, 1/))
    endif
    if (status /= nf90_noerr) call Handle_Error(status, 'nf_get_var_lsm')
    status = nf90_get_var(ncid_mask, varid_lat, Latitude)        
    if (status /= nf90_noerr) call Handle_Error(status, 'nf_get_var_lat')
    status = nf90_get_var(ncid_mask, varid_lon, Longitude)    
    if (status /= nf90_noerr) call Handle_Error(status, 'nf_get_var_lon')
    status = nf90_close(ncid_mask)
    if (status /= nf90_noerr) call Handle_Error(status, 'nf_close_mask')

    if (trim(domain) == 'Antarctica') then
        ! Open ice shelf mask
        status = nf90_open(trim(path_forcing_mask)//trim(fname_shelfmask), 0, ncid_shelf)
        if (status /= nf90_noerr) call Handle_Error(status, 'nf_open_shelfmask')
        status = nf90_inq_varid(ncid_shelf, "ISM", varid_ism)
        if (status /= nf90_noerr) call Handle_Error(status, 'nf_inq_varid_ism')
        status = nf90_get_var(ncid_shelf, varid_ism, ISM, start=(/1, 1/), count=(/Nlon, Nlat/))
        if (status /= nf90_noerr) call Handle_Error(status, 'nf_get_var_ism')
        status = nf90_close(ncid_shelf)
        if (status /= nf90_noerr) call Handle_Error(status, 'nf_close_shelfmask')
    else
        ! No ice shelves in Greenland
        ISM(:,:) = 0
    endif

    Nlon_file = Nlon
    Nlat_file = Nlat
    Nz_file = 1
    Nt_file = 1
    start_lon = 1
    start_lat = 1
    start_z = 1
    start_t = 1

    print *, 'start_lon, start_lat, start_z, start_t, Nlon_file, Nlat_file, Nz_file, Nt_file'
    print *, start_lon, start_lat, start_z, start_t, Nlon_file, Nlat_file, Nz_file, Nt_file
    print *, ' '

    call Load_File(path_forcing_averages, "snowmelt", suffix_forcing_averages, start_lon, start_lat, start_z, start_t, Nlon_file, Nlat_file, Nz_file, Nt_file, AveMelt)
    call Load_File(path_forcing_averages, "precip", suffix_forcing_averages, start_lon, start_lat, start_z, start_t, Nlon_file, Nlat_file, Nz_file, Nt_file, AveAcc)
    call Load_File(path_forcing_averages, "ff10m", suffix_forcing_averages, start_lon, start_lat, start_z, start_t, Nlon_file, Nlat_file, Nz_file, Nt_file, AveWind)
    call Load_File(path_forcing_averages, "tskin", suffix_forcing_averages, start_lon, start_lat, start_z, start_t, Nlon_file, Nlat_file, Nz_file, Nt_file, AveTsurf)
    call Load_File(path_forcing_averages, "evap", suffix_forcing_averages, start_lon, start_lat, start_z, start_t, Nlon_file, Nlat_file, Nz_file, Nt_file, AveSubl)
    call Load_File(path_forcing_averages, "sndiv", suffix_forcing_averages, start_lon, start_lat, start_z, start_t, Nlon_file, Nlat_file, Nz_file, Nt_file, AveSnowDrif)

    ! Convert units of average forcing from [mm w.e./s] to [mm w.e./yr]
    AveAcc = (AveAcc+AveSubl-AveSnowDrif) * (365.*24.*3600.)
    !AveAcc = AveAcc * (365.*24.*3600.)
    AveMelt = AveMelt * (365.*24.*3600.)
    
end subroutine Load_Ave_Forcing


! *******************************************************


subroutine Load_TimeSeries_Forcing(SnowMelt, PreTot, PreSol,PreLiq, Sublim, SnowDrif, TempSurf, FF10m, Nt_forcing, &
    ind_lon, ind_lat, path_forcing_timeseries, prefix_forcing_timeseries, dtobs, Nlon, Nlat_timeseries, Nlon_timeseries)
    !*** Load time series forcing data from netCDF files !***

    ! declare arguments
    integer, intent(in) :: Nt_forcing, ind_lon, ind_lat, dtobs, Nlon, Nlat_timeseries, Nlon_timeseries
    double precision, dimension(Nt_forcing), intent(out) :: SnowMelt, PreTot, PreSol,PreLiq, Sublim, SnowDrif, TempSurf, FF10m
    character*255, intent(in) :: path_forcing_timeseries, prefix_forcing_timeseries

    ! declare local variables
    integer :: ind_t, ind_latfile, ind_lonfile, numb_lat_loc, numb_lon_loc, numb_lon_tot, fnumb_int, start_lon, start_lat, start_z, start_t, Nlon_file, Nlat_file, Nz_file, Nt_file
    double precision :: remove_Psol, remove_Pliq, remove_Ptot, remove_Melt
    character*255 :: fnumb, suffix_forcing_timeseries
    
    ! determine indices of the forcing file
    ind_latfile = mod(ind_lat, Nlat_timeseries)
    if (ind_latfile==0) ind_latfile = Nlat_timeseries
    ind_lonfile = mod(ind_lon, Nlon_timeseries)
    if (ind_lonfile==0) ind_lonfile = Nlon_timeseries

    ! determine number of the forcing file
    numb_lat_loc = floor((real(ind_lat)-0.001)/real(Nlat_timeseries))+1
    numb_lon_loc = floor((real(ind_lon)-0.001)/real(Nlon_timeseries))+1
    numb_lon_tot = floor((real(Nlon)-0.001)/real(Nlon_timeseries))+1
    fnumb_int = numb_lon_loc + (numb_lat_loc-1)*numb_lon_tot
    if (fnumb_int<=9) write(fnumb,'(I1)') fnumb_int
    if (fnumb_int>=10) write(fnumb,'(I2)') fnumb_int

    suffix_forcing_timeseries = trim(prefix_forcing_timeseries)//trim(fnumb)//".nc"
    
    Nlon_file = 1
    Nlat_file = 1
    Nz_file = 1
    Nt_file = Nt_forcing
    start_lon = ind_lonfile
    start_lat = ind_latfile
    start_z = 1
    start_t = 1

    print *, 'start_lon, start_lat, start_z, start_t, Nlon_file, Nlat_file, Nz_file, Nt_file'
    print *, start_lon, start_lat, start_z, start_t, Nlon_file, Nlat_file, Nz_file, Nt_file
    print *, ' '

    call Load_File(path_forcing_timeseries, "snowmelt", suffix_forcing_timeseries, start_lon, start_lat, start_z, start_t, Nlon_file, Nlat_file, Nz_file, Nt_file, SnowMelt)
    call Load_File(path_forcing_timeseries, "precip", suffix_forcing_timeseries, start_lon, start_lat, start_z, start_t, Nlon_file, Nlat_file, Nz_file, Nt_file, PreTot)
    call Load_File(path_forcing_timeseries, "snowfall", suffix_forcing_timeseries, start_lon, start_lat, start_z, start_t, Nlon_file, Nlat_file, Nz_file, Nt_file, PreSol)
    call Load_File(path_forcing_timeseries, "ff10m", suffix_forcing_timeseries, start_lon, start_lat, start_z, start_t, Nlon_file, Nlat_file, Nz_file, Nt_file, FF10m)
    call Load_File(path_forcing_timeseries, "tskin", suffix_forcing_timeseries, start_lon, start_lat, start_z, start_t, Nlon_file, Nlat_file, Nz_file, Nt_file, TempSurf)
    call Load_File(path_forcing_timeseries, "evap", suffix_forcing_timeseries, start_lon, start_lat, start_z, start_t, Nlon_file, Nlat_file, Nz_file, Nt_file, Sublim)
    call Load_File(path_forcing_timeseries, "sndiv", suffix_forcing_timeseries, start_lon, start_lat, start_z, start_t, Nlon_file, Nlat_file, Nz_file, Nt_file, SnowDrif)

    ! Convert units from [mm w.e. s-1] to [mm w.e. per input time step]
    do ind_t = 1, Nt_forcing
        SnowMelt(ind_t) = SnowMelt(ind_t) * dtobs
        PreTot(ind_t) = PreTot(ind_t) * dtobs
        PreSol(ind_t) = PreSol(ind_t) * dtobs
        Sublim(ind_t) = Sublim(ind_t) * dtobs
        SnowDrif (ind_t) = SnowDrif (ind_t) * dtobs
        PreLiq(ind_t) = PreLiq(ind_t) * dtobs
    end do

    ! Remove rounding errors in the forcing files introduced when creating input files
    remove_Psol = 0.
    remove_Pliq = 0.
    remove_Ptot = 0.
    remove_Melt = 0.

    do ind_t = 1, Nt_forcing
        if (SnowMelt(ind_t) <1.e-04) then
            remove_Melt = remove_Melt + SnowMelt(ind_t)
            SnowMelt(ind_t) = 0.     
        endif
        
        if (PreSol(ind_t) < 1.e-04) then
            remove_Psol = remove_Psol + PreSol(ind_t)
            PreSol(ind_t) = 0.
        endif
    
        if (PreTot(ind_t) < 1.e-04) then
            remove_Ptot = remove_Ptot + PreTot(ind_t)
            PreTot(ind_t) = 0.
        endif
    
        if (TempSurf(ind_t) > 267.) then
            PreLiq(ind_t) = PreTot(ind_t)-PreSol(ind_t)
            if (PreLiq(ind_t) < 1.e-04) then
                remove_Pliq = remove_Pliq + PreLiq(ind_t)
                PreLiq(ind_t) = 0.
                PreSol(ind_t) = PreTot(ind_t)
            endif
        else
            PreSol(ind_t) = PreTot(ind_t)
        endif
    end do
    
    print *, "remove_Psol: ", remove_Psol
    print *, "remove_Pliq: ", remove_Pliq
    print *, "remove_Ptot: ", remove_Ptot
    print *, "remove_Melt: ", remove_Melt
    print *, ' '

end subroutine Load_TimeSeries_Forcing
    
    
! *******************************************************


subroutine Init_From_File(ind_z_max, ind_z_surf, Rho, M, T, Depth, Mlwc, DZ, path_in_restart, fname_out_ini)
    !*** This subroutine reads the initial conditions from a netCDF ***!

    ! declare arguments
    integer, intent(in) :: ind_z_max
    integer, intent(inout) :: ind_z_surf
    double precision, dimension(ind_z_max), intent(out) :: Rho, M, T, Depth, Mlwc
    double precision, dimension(ind_z_max), intent(inout) :: DZ
    character*255, intent(in) :: path_in_restart, fname_out_ini

    ! declare local variables
    integer :: ind_z, status, ncid, ID(5), LayerID
    character*255 :: name_inifile
        
    print *, 'trim(path_in_restart)//trim(fname_out_ini)'
    print *, trim(path_in_restart)//trim(fname_out_ini)
    print *, ' '

    ! Open the snowmelt netCDF file
    status = nf90_open(trim(path_in_restart)//trim(fname_out_ini), 0, ncid)
    if (status /= nf90_noerr) call Handle_Error(status,'nf_open1')

    ! Get ID of the variables in the NetCDF files
    status = nf90_inq_varid(ncid, "dens", ID(1))
    if (status /= nf90_noerr) call Handle_Error(status, 'nf_inq_varid1')
    status = nf90_inq_varid(ncid, "temp", ID(2))
    if (status /= nf90_noerr) call Handle_Error(status, 'nf_inq_varid2')
    status = nf90_inq_varid(ncid, "mass", ID(3))
    if (status /= nf90_noerr) call Handle_Error(status, 'nf_inq_varid3')
    status = nf90_inq_varid(ncid, "depth", ID(4))
    if (status /= nf90_noerr) call Handle_Error(status, 'nf_inq_varid4')
    status = nf90_inq_varid(ncid, "lwc", ID(5))
    if (status /= nf90_noerr) call Handle_Error(status, 'nf_inq_varid5')

    ! Get dimension of the array
    status = nf90_inq_dimid(ncid, "layer", LayerID)
    if (status /= nf90_noerr) call Handle_Error(status, 'nf_inq_dimid1')
    status = nf90_inquire_dimension(ncid, LayerID, len=ind_z_surf)
    if (status /= nf90_noerr) call Handle_Error(status, 'nf_inq_dim1')
    
    ! Get all variables from the netCDF files
    status = nf90_get_var(ncid, ID(1), Rho(1:ind_z_surf), start=(/1/), count=(/ind_z_surf/))
    if (status /= nf90_noerr) call Handle_Error(status, 'nf_get_var1')
    status = nf90_get_var(ncid, ID(2), T(1:ind_z_surf), start=(/1/), count=(/ind_z_surf/))
    if (status /= nf90_noerr) call Handle_Error(status, 'nf_get_var2')
    status = nf90_get_var(ncid, ID(3), M(1:ind_z_surf), start=(/1/), count=(/ind_z_surf/))
    if (status /= nf90_noerr) call Handle_Error(status, 'nf_get_var3')
    status = nf90_get_var(ncid, ID(4), Depth(1:ind_z_surf), start=(/1/), count=(/ind_z_surf/))
    if (status /= nf90_noerr) call Handle_Error(status, 'nf_get_var4')
    status = nf90_get_var(ncid, ID(5), Mlwc(1:ind_z_surf), start=(/1/), count=(/ind_z_surf/))
    if (status /= nf90_noerr) call Handle_Error(status, 'nf_get_var5')
    
    ! Close all netCDF files    
    status = nf90_close(ncid)
    if (status /= nf90_noerr) call Handle_Error(status,'nf_close1')

    DZ(ind_z_surf) = Depth(ind_z_surf) * 2.
    do ind_z = (ind_z_surf-1), 1, -1
        DZ(ind_z) = (Depth(ind_z) - Depth(ind_z+1) - 0.5*DZ(ind_z+1)) * 2.
    end do
    
end subroutine Init_From_File


! *******************************************************


subroutine Load_File(path_forcing, variable_name, suffix_file, start_lon, start_lat, start_z, start_t, Nlon, Nlat, Nz, Nt, myoutput)
    !*** This subroutine reads the forcing data from a netCDF file !***

    ! declare arguments
    integer, intent(in) :: Nlon, Nlat, Nz, Nt, start_lon, start_lat, start_z, start_t
    double precision, dimension(Nlon,Nlat,Nz,Nt), intent(out) :: myoutput
    character*255, intent(in) :: path_forcing, suffix_file
    character(len=*), intent(in) :: variable_name

    ! declare local variables
    integer :: status=0, ncid=0, varid=0

    print *, trim(path_forcing)//trim(variable_name)//trim(suffix_file)
    status = nf90_open(trim(path_forcing)//trim(variable_name)//trim(suffix_file), 0, ncid)
    if (status /= nf90_noerr) call Handle_Error(status, 'nf_open_'//trim(variable_name))
    print *, 'opened file succesfully'
    status = nf90_inq_varid(ncid, trim(variable_name), varid)
    if (status /= nf90_noerr) call Handle_Error(status, 'nf_inq_varid'//trim(variable_name))
    print *, 'inquired variable succesfully'
    status = nf90_get_var(ncid, varid, myoutput, start=(/start_lon, start_lat, start_z, start_t/), count=(/Nlon, Nlat, Nz, Nt/))
    if (status /= nf90_noerr) call Handle_Error(status, 'nf_get_var_'//trim(variable_name))
    print *, 'obtained variable succesfully'
    status = nf90_close(ncid)
    if (status /= nf90_noerr) call Handle_Error(status, 'nf_close_'//trim(variable_name))
    print *, 'closed file succesfully'
    print *, ' '

end subroutine Load_File


! *******************************************************


subroutine Handle_Error(stat, msg)
    !*** error stop for netCDF

    integer, intent(in) :: stat
    character(len=*), intent(in) :: msg

    if (stat /= nf90_noerr) then
        print *, 'netCDF error (', msg, '): ', nf90_strerror(stat)
        stop
    endif

end subroutine Handle_Error


end module openNetCDF
