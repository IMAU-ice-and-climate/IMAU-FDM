module openNetCDF
    !*** subroutines for loading forcing data

    use netcdf, only : nf90_open,nf90_inq_varid, nf90_inq_dimid, nf90_inquire_dimension, & 
        nf90_close,nf90_get_var,nf90_noerr, nf90_strerror, nf90_max_name, nf90_inquire, &
        nf90_inquire_variable
    
    use model_settings

    implicit none
    private

    type, public :: TimeBounds
        character(5) :: start_ts_year
        character(5) :: end_ts_year
        character(30) :: model_first_timestep
        character(30) :: model_last_timestep
    end type TimeBounds

    public :: Load_Mask, Load_Ave_Forcing, Load_TimeSeries_Forcing, Handle_Error, Restart_From_Spinup, Restart_From_Run, read_avg_dimensions, read_time_bounds
    ! public :: read_avg_dimensions, Handle_Error
contains


! *******************************************************

subroutine read_avg_dimensions(netcdf_file, n_lon, n_lat)

    integer :: status, ncid, dim_id
    integer, intent(out) :: n_lon, n_lat
    character(len = *), intent(in) :: netcdf_file
    character(len=nf90_max_name) :: dim_name

    call handle_error(nf90_open(netcdf_file, 0, ncid), "open_avg_file")
    call handle_error(nf90_inq_dimid(ncid, "rlat", dim_id), "find_rlat")
    call handle_error(nf90_inquire_dimension(ncid, dim_id, dim_name, n_lat), "get_nlat")
    call handle_error(nf90_inq_dimid(ncid, "rlon", dim_id), "find_rlon")
    call handle_error(nf90_inquire_dimension(ncid, dim_id, dim_name, n_lon), "get_nlon")
    status = nf90_close(ncid)

end subroutine read_avg_dimensions

function read_time_bounds(netcdf_file) result(time_bounds)
    integer :: ncid, ndim
    integer :: nvar, nattr, unlim_dim, format_num
    integer :: idim
    integer :: bnds_varid

    integer, allocatable :: bounds(:, :)
    integer :: counts(2)
    integer :: bound_dimids(2)
    integer :: start_ts_year, end_ts_year, start_date, end_date

    type(TimeBounds), allocatable :: time_bounds

    character(len = *), intent(in) :: netcdf_file
    character(len=nf90_max_name) :: dim_name, var_name

    allocate(time_bounds)
  
    call handle_error(nf90_open("data/evap_FGRN055_era055_1939-2023_p24-007.nc", 0, ncid), "open_netcdf_file")
    call handle_error(nf90_inquire(ncid, ndim, nvar, nattr, unlim_dim, format_num), "inquire_net_file")

    ! print *, ndim, nvar, nattr, unlim_dim, format_num

    call Handle_Error(nf90_inq_varid(ncid, "date_bnds", bnds_varid), "Reading var id day_bnds")
    call Handle_Error(nf90_inquire_variable(ncid, bnds_varid, dimids=bound_dimids), "Getting dim_ids for bounds")
    do idim = 1, 2
        call Handle_Error(nf90_inquire_dimension(ncid, bound_dimids(idim), len=counts(idim)), "Reading time bound counts")
    end do
    allocate(bounds(counts(1), counts(2)))
    call Handle_Error(nf90_get_var(ncid, bnds_varid, bounds), "Reading time bounds")

    start_date = bounds(1, 1)
    end_date = bounds(1, counts(2))
    start_ts_year = start_date / 10000
    end_ts_year = bounds(1, counts(2)) / 10000
    write(time_bounds%model_first_timestep, "(i4, '-', i2.2, '-', i2.2, 'T00:00:00')") start_ts_year, mod(start_date/100, 100), mod(start_date, 100)
    write(time_bounds%model_last_timestep, "(i4, '-', i2.2, '-', i2.2, 'T21:00:00')") end_ts_year, mod(end_date/100, 100), mod(end_date, 100)
    write(time_bounds%start_ts_year, "(i4)") start_ts_year
    write(time_bounds%end_ts_year, "(i4)") end_ts_year
    call Handle_Error(nf90_close(ncid), "Closing netcdf4 file")
end function

subroutine Load_Mask(LSM, Nlat, Nlon, Latitude, Longitude, ISM, domain)
    
    integer :: status, ncid(1), ID(5)
    double precision, dimension(config%forcing_dimensions%Nlon,config%forcing_dimensions%Nlat) :: LSM, ISM, Latitude, Longitude
    character*255 :: pad

    pad = trim(reference_dir)//trim(fname_mask)

    write(log_unit, *) "Path to mask: ", trim(pad)
    write(log_unit, *) " "

    ! open ice mask
    status = nf90_open(trim(pad),0,ncid(1))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_open_icemask')

    ! import icemask, lat, lon, and ice shelves as needed
    status = nf90_inq_varid(ncid(1),"LSM",ID(1))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_inq_icemask')    
    status = nf90_inq_varid(ncid(1),"lat",ID(2))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_inq_icemask_lat')
    status = nf90_inq_varid(ncid(1),"lon",ID(3))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_inq_icemask_lon')
    
    status  = nf90_get_var(ncid(1),ID(1),LSM,start=(/1,1,1,1/), &
        count=(/config%forcing_dimensions%Nlon,config%forcing_dimensions%Nlat,1,1/))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_lsm')
    status  = nf90_get_var(ncid(1),ID(2),Latitude)        
    if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_lat')
    status  = nf90_get_var(ncid(1),ID(3),Longitude)    
    if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_lon')

    if (domain == "ANT27") then
        status = nf90_inq_varid(ncid(1),iceshelf_var,ID(5))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_inq_iceshelf')
        status  = nf90_get_var(ncid(1),ID(6),ISM,start=(/1,1/), &
            count=(/config%forcing_dimensions%Nlon,config%forcing_dimensions%Nlat/))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_iceshelfmask')
    else
        ! No Ice Shelves in Greenland / no ice shelf variable defined 
        ISM(:,:) = 0
    end if

    status = nf90_close(ncid(1))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_close_icemask')  

end subroutine Load_Mask

subroutine Load_Ave_Forcing(AveTsurf, AveAcc, AveWind, AveMelt)
    
    integer :: status,ncid(50),ID(50),i,j
    double precision, dimension(config%forcing_dimensions%Nlon,config%forcing_dimensions%Nlat) :: AveTsurf,AveAcc,AveWind,AveSubl, &
        AveSnowDrif,AveMelt

    write(log_unit, *) "Path to averages: ", trim(input_averages_dir)//"VAR"//trim(suffix_forcing_averages)
    write(log_unit, *) " "

    status = nf90_open(trim(input_averages_dir)//"snowmelt"//trim(suffix_forcing_averages),0,ncid(1))
    if(status /= nf90_noerr) call Handle_Error(status,'ave_var_open1')
    status = nf90_open(trim(input_averages_dir)//"precip"//trim(suffix_forcing_averages),0,ncid(2))
    if(status /= nf90_noerr) call Handle_Error(status,'ave_var_open2')
    status = nf90_open(trim(input_averages_dir)//"ff10m"//trim(suffix_forcing_averages),0,ncid(3))
    if(status /= nf90_noerr) call Handle_Error(status,'ave_var_open3')
    status = nf90_open(trim(input_averages_dir)//"tskin"//trim(suffix_forcing_averages),0,ncid(4))
    if(status /= nf90_noerr) call Handle_Error(status,'ave_var_open4')
    status = nf90_open(trim(input_averages_dir)//"evap"//trim(suffix_forcing_averages),0,ncid(5))
    if(status /= nf90_noerr) call Handle_Error(status,'ave_var_open5')
    status = nf90_open(trim(input_averages_dir)//"sndiv"//trim(suffix_forcing_averages),0,ncid(6))
    if(status /= nf90_noerr) call Handle_Error(status,'ave_var_open7')

    status = nf90_inq_varid(ncid(1),"snowmelt",ID(1))
    if(status /= nf90_noerr) call Handle_Error(status,'ave_var_inq_snowmelt')
    status = nf90_inq_varid(ncid(2),"precip",ID(2))
    if(status /= nf90_noerr) call Handle_Error(status,'ave_var_inq_precip')
    status = nf90_inq_varid(ncid(3),"ff10m",ID(3))
    if(status /= nf90_noerr) call Handle_Error(status,'ave_var_inq_ff10m')
    status = nf90_inq_varid(ncid(4),"tskin",ID(4))
    if(status /= nf90_noerr) call Handle_Error(status,'ave_var_inq_tskin')
    status = nf90_inq_varid(ncid(5),"evap",ID(5))
    if(status /= nf90_noerr) call Handle_Error(status,'ave_var_inq_evap')
    status = nf90_inq_varid(ncid(6),"sndiv",ID(6))
    if(status /= nf90_noerr) call Handle_Error(status,'ave_var_inq_sndiv')

    status  = nf90_get_var(ncid(1),ID(1),AveMelt,start=(/1,1,1,1/), &
        count=(/config%forcing_dimensions%Nlon,config%forcing_dimensions%Nlat,1,1/))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_avemelt')
    status  = nf90_get_var(ncid(2),ID(2),AveAcc,start=(/1,1,1,1/), &
        count=(/config%forcing_dimensions%Nlon,config%forcing_dimensions%Nlat,1,1/))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_aveacc')
    status  = nf90_get_var(ncid(3),ID(3),AveWind,start=(/1,1,1,1/), &
        count=(/config%forcing_dimensions%Nlon,config%forcing_dimensions%Nlat,1,1/))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_avewind')
    status  = nf90_get_var(ncid(4),ID(4),AveTsurf,start=(/1,1,1,1/), &
        count=(/config%forcing_dimensions%Nlon,config%forcing_dimensions%Nlat,1,1/))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_avetsurf')
    status  = nf90_get_var(ncid(5),ID(5),AveSubl,start=(/1,1,1,1/), &
        count=(/config%forcing_dimensions%Nlon,config%forcing_dimensions%Nlat,1,1/))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_avesubl')
    status  = nf90_get_var(ncid(6),ID(6),AveSnowDrif,start=(/1,1,1,1/), &
        count=(/config%forcing_dimensions%Nlon,config%forcing_dimensions%Nlat,1,1/))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_avesndiv')

    ! Close all netCDF files    

    status = nf90_close(ncid(1))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_ave_close1')
    status = nf90_close(ncid(2))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_ave_close2')
    status = nf90_close(ncid(3))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_ave_close3')
    status = nf90_close(ncid(4))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_ave_close4')
    status = nf90_close(ncid(5))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_ave_close5')
    status = nf90_close(ncid(6))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_ave_close6')

    ! Convert units from [mm w.e./s] to [mm w.e./yr]
    do i = 1, config%forcing_dimensions%Nlon
        do j = 1, config%forcing_dimensions%Nlat
            AveAcc(i,j) = (AveAcc(i,j)+AveSubl(i,j)-AveSnowDrif(i,j)) &
                * (const%seconds_per_year)
            AveMelt(i,j) = AveMelt(i,j) * (const%seconds_per_year)
        end do
    end do
    
end subroutine Load_Ave_Forcing


! *******************************************************


subroutine Load_TimeSeries_Forcing(SnowMelt, PreTot, PreSol,PreLiq, Sublim, SnowDrif, TempSurf, FF10m, &
    ind_lon, ind_lat)

    integer :: status, ind_t, ncid(50), ind_lon, ind_lat, ID(50)
    double precision :: remove_Psol,remove_Pliq,remove_Ptot,remove_Melt
    double precision, dimension(:) :: SnowMelt,PreTot,PreSol,PreLiq,Sublim,TempSurf, &
        SnowDrif,FF10m

    integer :: latfile, lonfile, fnumb_i
    character*255 :: add, fnumb

    if ( trim(project_name) == "example" ) then

        latfile = 1
        lonfile = 1
        fnumb = ""
        
        add = trim(prefix_forcing_timeseries)//trim(suffix_forcing_timeseries)
        
    else !running point

        latfile = ind_lat
        lonfile = mod(ind_lon,config%forcing_dimensions%Nlon_timeseries)

        if (lonfile==0) lonfile = config%forcing_dimensions%Nlon_timeseries
        fnumb_i = floor((real(ind_lon)-0.001)/real(config%forcing_dimensions%Nlon_timeseries))+1
        if (fnumb_i<=9) write(fnumb,'(I1)') fnumb_i
        if (fnumb_i>=10) write(fnumb,'(I2)') fnumb_i

        add = trim(prefix_forcing_timeseries)//trim(fnumb)//trim(suffix_forcing_timeseries)

    end if
    

    write(log_unit, *) "Looking for timeseries at: "
    write(log_unit, *) trim(input_timeseries_dir)//"VAR"//trim(add)
    
    ! Open the snowmelt netCDF file
    status = nf90_open(trim(input_timeseries_dir)//"snowmelt"//trim(add),0,ncid(1))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_ts_open1')
    status = nf90_open(trim(input_timeseries_dir)//"precip"//trim(add),0,ncid(2))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_ts_open2')
    status = nf90_open(trim(input_timeseries_dir)//"snowfall"//trim(add),0,ncid(3))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_ts_open3')
    status = nf90_open(trim(input_timeseries_dir)//"evap"//trim(add),0,ncid(4))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_ts_open4')
    status = nf90_open(trim(input_timeseries_dir)//"tskin"//trim(add),0,ncid(5))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_ts_open5')
    status = nf90_open(trim(input_timeseries_dir)//"sndiv"//trim(add),0,ncid(6))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_ts_open6')
    status = nf90_open(trim(input_timeseries_dir)//"ff10m"//trim(add),0,ncid(7))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_ts_open7')
    
    !Get ID of the variables in the NetCDF files
    status = nf90_inq_varid(ncid(1),"snowmelt",ID(1))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_ts_inq_varid1')
    status = nf90_inq_varid(ncid(2),"precip",ID(2))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_ts_inq_varid2')
    status = nf90_inq_varid(ncid(3),"snowfall",ID(3))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_ts_inq_varid3')
    status = nf90_inq_varid(ncid(4),"evap",ID(4))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_ts_inq_varid4')
    status = nf90_inq_varid(ncid(5),"tskin",ID(5))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_ts_inq_varid5')
    status = nf90_inq_varid(ncid(6),"sndiv",ID(6))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_ts_inq_varid6')
    status = nf90_inq_varid(ncid(7),"ff10m",ID(7))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_ts_inq_varid7')

    ! Get all variables from the netCDF files
    status  = nf90_get_var(ncid(1),ID(1),SnowMelt,start=(/lonfile,latfile,1,1/), &
        count=(/1,1,1,config%forcing_dimensions%Nt_forcing/))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var1')
    write(log_unit, *) "Read snowmelt..."
    status  = nf90_get_var(ncid(2),ID(2),PreTot,start=(/lonfile,latfile,1,1/), &
        count=(/1,1,1,config%forcing_dimensions%Nt_forcing/))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var2')
    write(log_unit, *) "Read precipitation..."
    status  = nf90_get_var(ncid(3),ID(3),PreSol,start=(/lonfile,latfile,1,1/), &
        count=(/1,1,1,config%forcing_dimensions%Nt_forcing/))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var3')
    write(log_unit, *) "Read snowfall..."
    status  = nf90_get_var(ncid(4),ID(4),Sublim,start=(/lonfile,latfile,1,1/), &
        count=(/1,1,1,config%forcing_dimensions%Nt_forcing/))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var4')
    write(log_unit, *) "Read sublimation..."
    status  = nf90_get_var(ncid(5),ID(5),TempSurf,start=(/lonfile,latfile,1,1/), &
        count=(/1,1,1,config%forcing_dimensions%Nt_forcing/))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var5')
    write(log_unit, *) "Read skin temperature..."
    status  = nf90_get_var(ncid(6),ID(6),SnowDrif,start=(/lonfile,latfile,1,1/), &
        count=(/1,1,1,config%forcing_dimensions%Nt_forcing/))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var6')
    write(log_unit, *) "Read snow drift..."
    status  = nf90_get_var(ncid(7),ID(7),FF10m,start=(/lonfile,latfile,1,1/), &
        count=(/1,1,1,config%forcing_dimensions%Nt_forcing/))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var7')
    write(log_unit, *) "Read wind speed..."
    write(log_unit, *) ' '

    ! Close all netCDF files    
    status = nf90_close(ncid(1))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_close1')
    status = nf90_close(ncid(2))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_close2')
    status = nf90_close(ncid(3))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_close3')
    status = nf90_close(ncid(4))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_close4')
    status = nf90_close(ncid(5))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_close5')
    status = nf90_close(ncid(6))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_close6')
    status = nf90_close(ncid(7))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_close7')

    remove_Psol = 0.
    remove_Pliq = 0.
    remove_Ptot = 0.
    remove_Melt = 0.

    ! from per second to per observation time period (e.g., 3 hourly)
    do ind_t = 1, config%forcing_dimensions%Nt_forcing
        SnowMelt(ind_t) = SnowMelt(ind_t) * config%forcing_dimensions%dtobs
        PreTot(ind_t) = PreTot(ind_t) * config%forcing_dimensions%dtobs
        PreSol(ind_t) = PreSol(ind_t) * config%forcing_dimensions%dtobs
        Sublim(ind_t) = Sublim(ind_t) * config%forcing_dimensions%dtobs
        SnowDrif(ind_t) = SnowDrif(ind_t) * config%forcing_dimensions%dtobs
        PreLiq(ind_t) = PreLiq(ind_t) * config%forcing_dimensions%dtobs
    end do

    do ind_t = 1, config%forcing_dimensions%Nt_forcing
        if (SnowMelt(ind_t) < config%model_choices%ts_minimum) then
            remove_Melt = remove_Melt + SnowMelt(ind_t)
            SnowMelt(ind_t) = 0.     
        endif
        
        if (PreSol(ind_t) < config%model_choices%ts_minimum) then
            remove_Psol = remove_Psol + PreSol(ind_t)
            PreSol(ind_t) = 0.
        endif
    
        if (PreTot(ind_t) < config%model_choices%ts_minimum) then
            remove_Ptot = remove_Ptot + PreTot(ind_t)
            PreTot(ind_t) = 0.
        endif
    
        if (TempSurf(ind_t) > 267.) then
            PreLiq(ind_t) = PreTot(ind_t)-PreSol(ind_t)
            if (PreLiq(ind_t) < config%model_choices%ts_minimum) then
                remove_Pliq = remove_Pliq + PreLiq(ind_t)
                PreLiq(ind_t) = 0.
                PreSol(ind_t) = PreTot(ind_t)
            endif
        else
            PreSol(ind_t) = PreTot(ind_t)
        endif
    end do
    
    write(log_unit, *) "remove_Psol: ", remove_Psol
    write(log_unit, *) "remove_Pliq: ", remove_Pliq
    write(log_unit, *) "remove_Ptot: ", remove_Ptot
    write(log_unit, *) "remove_Melt: ", remove_Melt
    write(log_unit, *) ' '

end subroutine Load_TimeSeries_Forcing
    
    
! *******************************************************


subroutine Restart_From_Spinup(Rho, M, T, Depth, Mlwc, DZ, DenRho, Refreeze)

    integer :: ind_z, status, ncid(50), ID(50), LayerID
    
    double precision, dimension(ind_z_max) :: Rho, M, T, Depth, Mlwc, DZ, DenRho, Refreeze
    
    character*255 :: pad
    
    pad = trim(restart_dir)//"spinup/"//trim(fname_restart_from_spinup)
    
    write(log_unit, *) "Path of spinup restart file: "
    write(log_unit, *) trim(pad)
    write(log_unit, *) ' '

    ! Open the restart netCDF file
    status = nf90_open(trim(pad),0,ncid(1))
    if(status /= nf90_noerr) call Handle_Error(status,'load_restart_open1')

    !Get ID of the variables in the NetCDF files
    status = nf90_inq_varid(ncid(1),"dens",ID(1))
    if(status /= nf90_noerr) call Handle_Error(status,'load_restart_inq_varid1')
    status = nf90_inq_varid(ncid(1),"temp",ID(2))
    if(status /= nf90_noerr) call Handle_Error(status,'load_restart_inq_varid2')
    status = nf90_inq_varid(ncid(1),"mass",ID(3))
    if(status /= nf90_noerr) call Handle_Error(status,'load_restart_inq_varid3')
    status = nf90_inq_varid(ncid(1),"depth",ID(4))
    if(status /= nf90_noerr) call Handle_Error(status,'load_restart_inq_varid4')
    status = nf90_inq_varid(ncid(1),"lwc",ID(5))
    if(status /= nf90_noerr) call Handle_Error(status,'load_restart_inq_varid5')

    ! Get dimension of the array
    status = nf90_inq_dimid(ncid(1),"layer",LayerID)
    if(status /= nf90_noerr) call Handle_Error(status,'load_restart_inq_dimid1')
    status = nf90_inquire_dimension(ncid(1),LayerID,len=ind_z_surf)
    if(status /= nf90_noerr) call Handle_Error(status,'load_restart_inq_dim1')
    
    ! Get all variables from the netCDF files
    status  = nf90_get_var(ncid(1),ID(1),Rho(1:ind_z_surf),start=(/1/),count=(/ind_z_surf/))
    if(status /= nf90_noerr) call Handle_Error(status,'load_restart_get_var1')
    status  = nf90_get_var(ncid(1),ID(2),T(1:ind_z_surf),start=(/1/),count=(/ind_z_surf/))
    if(status /= nf90_noerr) call Handle_Error(status,'load_restart_get_var2')
    status  = nf90_get_var(ncid(1),ID(3),M(1:ind_z_surf),start=(/1/),count=(/ind_z_surf/))
    if(status /= nf90_noerr) call Handle_Error(status,'load_restart_get_var3')
    status  = nf90_get_var(ncid(1),ID(4),Depth(1:ind_z_surf),start=(/1/),count=(/ind_z_surf/))
    if(status /= nf90_noerr) call Handle_Error(status,'load_restart_get_var4')
    status  = nf90_get_var(ncid(1),ID(5),Mlwc(1:ind_z_surf),start=(/1/),count=(/ind_z_surf/))
    if(status /= nf90_noerr) call Handle_Error(status,'load_restart_get_var5')
    
    ! Close all netCDF files    
    status = nf90_close(ncid(1))
    if(status /= nf90_noerr) call Handle_Error(status,'load_restart_close1')

    DZ(ind_z_surf) = Depth(ind_z_surf) * 2.
    do ind_z = (ind_z_surf-1), 1, -1
        DZ(ind_z) = (Depth(ind_z) - Depth(ind_z+1) - 0.5*DZ(ind_z+1)) * 2.
    end do

    ! Reset profiles that need to start at zero after the spin-up
    DenRho(:) = 0.
    Refreeze(:) = 0.
    
end subroutine Restart_From_Spinup


! *******************************************************

subroutine Restart_From_Run(prev_nt, Rho, M, T, Depth, Mlwc, DZ, Year, DenRho, Refreeze)

    integer :: prev_nt
    integer :: ind_z, status, ncid(50), ID(50), LayerID(2)

    double precision, dimension(ind_z_max) :: Rho, M, T, Depth, Mlwc, DZ, Year, DenRho, Refreeze

    character*255 :: pad

    pad = trim(restart_dir)//"run/"//trim(fname_restart_from_previous_run)
    
    write(log_unit, *) 'Path of run restart files: :'
    write(log_unit, *) trim(pad)
    write(log_unit, *) ' '

    ! Open the snowmelt netCDF file
    status = nf90_open(trim(pad),0,ncid(1))
    if(status /= nf90_noerr) call Handle_Error(status,'load_run_restart_open1')

    !Get ID of the variables in the NetCDF files
    status = nf90_inq_varid(ncid(1),"dens",ID(1))
    if(status /= nf90_noerr) call Handle_Error(status,'load_run_restart_inq_varid1')
    status = nf90_inq_varid(ncid(1),"temp",ID(2))
    if(status /= nf90_noerr) call Handle_Error(status,'load_run_restart_inq_varid2')
    status = nf90_inq_varid(ncid(1),"mass",ID(3))
    if(status /= nf90_noerr) call Handle_Error(status,'load_run_restart_inq_varid3')
    status = nf90_inq_varid(ncid(1),"depth",ID(4))
    if(status /= nf90_noerr) call Handle_Error(status,'load_run_restart_inq_varid4')
    status = nf90_inq_varid(ncid(1),"lwc",ID(5))
    if(status /= nf90_noerr) call Handle_Error(status,'load_run_restart_inq_varid5')
    status = nf90_inq_varid(ncid(1),"year",ID(6))
    if(status /= nf90_noerr) call Handle_Error(status,'load_run_restart_inq_varid6')
    status = nf90_inq_varid(ncid(1),"refreeze",ID(7))
    if(status /= nf90_noerr) call Handle_Error(status,'load_run_restart_inq_varid7')
    status = nf90_inq_varid(ncid(1),"denrho",ID(8))
    if(status /= nf90_noerr) call Handle_Error(status,'load_run_restart_inq_varid8')
    status = nf90_inq_varid(ncid(1),"prev_nt",ID(9))
    if(status /= nf90_noerr) call Handle_Error(status,'load_run_restart_inq_varid9')

    ! Get dimension of the array
    status = nf90_inq_dimid(ncid(1),"layer",LayerID(1))
    if(status /= nf90_noerr) call Handle_Error(status,'load_run_restart_inq_dimid1')
    status = nf90_inquire_dimension(ncid(1),LayerID(1),len=ind_z_surf)
    if(status /= nf90_noerr) call Handle_Error(status,'load_run_restart_inq_dim1')
    status = nf90_inq_dimid(ncid(1),"constant",LayerID(2))
    if(status /= nf90_noerr) call Handle_Error(status,'load_run_restart_inq_dimid2')
    status = nf90_inquire_dimension(ncid(1),LayerID(2))!,len=1)
    if(status /= nf90_noerr) call Handle_Error(status,'load_run_restart_inq_dim2')
    
    ! Get all variables from the netCDF files
    status  = nf90_get_var(ncid(1),ID(1),Rho(1:ind_z_surf),start=(/1/),count=(/ind_z_surf/))
    if(status /= nf90_noerr) call Handle_Error(status,'load_run_restart_get_var1')
    status  = nf90_get_var(ncid(1),ID(2),T(1:ind_z_surf),start=(/1/),count=(/ind_z_surf/))
    if(status /= nf90_noerr) call Handle_Error(status,'load_run_restart_get_var2')
    status  = nf90_get_var(ncid(1),ID(3),M(1:ind_z_surf),start=(/1/),count=(/ind_z_surf/))
    if(status /= nf90_noerr) call Handle_Error(status,'load_run_restart_get_var3')
    status  = nf90_get_var(ncid(1),ID(4),Depth(1:ind_z_surf),start=(/1/),count=(/ind_z_surf/))
    if(status /= nf90_noerr) call Handle_Error(status,'load_run_restart_get_var4')
    status  = nf90_get_var(ncid(1),ID(5),Mlwc(1:ind_z_surf),start=(/1/),count=(/ind_z_surf/))
    if(status /= nf90_noerr) call Handle_Error(status,'load_run_restart_get_var5')
    status  = nf90_get_var(ncid(1),ID(6),Year(1:ind_z_surf),start=(/1/),count=(/ind_z_surf/))
    if(status /= nf90_noerr) call Handle_Error(status,'load_run_restart_get_var6')
    status  = nf90_get_var(ncid(1),ID(7),Refreeze(1:ind_z_surf),start=(/1/),count=(/ind_z_surf/))
    if(status /= nf90_noerr) call Handle_Error(status,'load_run_restart_get_var7')
    status  = nf90_get_var(ncid(1),ID(8),DenRho(1:ind_z_surf),start=(/1/),count=(/ind_z_surf/))
    if(status /= nf90_noerr) call Handle_Error(status,'load_run_restart_get_var8')
    status  = nf90_get_var(ncid(1),ID(9),prev_nt)!,start=1,count=1)
    if(status /= nf90_noerr) call Handle_Error(status,'load_run_restart_get_var9')
    
    ! Close all netCDF files    
    status = nf90_close(ncid(1))
    if(status /= nf90_noerr) call Handle_Error(status,'load_run_restart_close1')

    DZ(ind_z_surf) = Depth(ind_z_surf) * 2.
    do ind_z = (ind_z_surf-1), 1, -1
        DZ(ind_z) = (Depth(ind_z) - Depth(ind_z+1) - 0.5*DZ(ind_z+1)) * 2.
    end do

    ! Reset profiles that need to start at zero after the spin-up
    ! DenRho(:) = 0.
    ! Refreeze(:) = 0.
    
end subroutine Restart_From_Run


! *******************************************************


subroutine Handle_Error(stat,msg)
    !*** error stop for netCDF

    integer, intent(in) :: stat
    character(len=*), intent(in) :: msg

    if(stat /= nf90_noerr) then
        write(log_unit, *) 'netCDF error (',msg,'): ', nf90_strerror(stat)
        stop
    endif

end subroutine Handle_Error


end module openNetCDF
