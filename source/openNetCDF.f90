module openNetCDF
    !*** subroutines for loading forcing data

    use netcdf, only : nf90_open,nf90_inq_varid, nf90_inq_dimid, nf90_inquire_dimension, & 
        nf90_close,nf90_get_var,nf90_noerr, nf90_strerror
    
    use model_settings

    implicit none
    private

    public :: Load_Mask, Load_Ave_Forcing, Load_TimeSeries_Forcing, Handle_Error, Restart_From_Spinup, Restart_From_Run

contains


! *******************************************************

subroutine Load_Mask(LSM, Nlat, Nlon, Latitude, Longitude, ISM, domain)
    
    integer :: status, ncid(1), ID(5), Nlat, Nlon
    double precision, dimension(Nlon,Nlat) :: LSM, ISM, Latitude, Longitude
    character*255 :: domain, pad

    pad = trim(path_forcing_mask)//trim(fname_mask)

    print *, "Path to mask: ", trim(pad)
    print *, " "

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
        count=(/Nlon,Nlat,1,1/))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_lsm')
    status  = nf90_get_var(ncid(1),ID(2),Latitude)        
    if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_lat')
    status  = nf90_get_var(ncid(1),ID(3),Longitude)    
    if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_lon')

    if (domain == "ANT27") then
        status = nf90_inq_varid(ncid(1),iceshelf_var,ID(5))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_inq_iceshelf')
        status  = nf90_get_var(ncid(1),ID(6),ISM,start=(/1,1/), &
            count=(/Nlon,Nlat/))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_iceshelfmask')
    else
        ! No Ice Shelves in Greenland / no ice shelf variable defined 
        ISM(:,:) = 0
    end if

    status = nf90_close(ncid(1))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_close_icemask')  

end subroutine Load_Mask

subroutine Load_Ave_Forcing(AveTsurf, AveAcc, AveWind, AveMelt, Nlat, Nlon)
    
    integer :: status,ncid(50),ID(50),Nlat,Nlon,i,j
    double precision, dimension(Nlon,Nlat) :: AveTsurf,AveAcc,AveWind,AveSubl, &
        AveSnowDrif,AveMelt

    print *, "Path to averages: ", trim(path_forcing_averages)//"VAR"//trim(suffix_forcing_averages)
    print *, " "

    status = nf90_open(trim(path_forcing_averages)//"snowmelt"//trim(suffix_forcing_averages),0,ncid(1))
    if(status /= nf90_noerr) call Handle_Error(status,'ave_var_open1')
    status = nf90_open(trim(path_forcing_averages)//"precip"//trim(suffix_forcing_averages),0,ncid(2))
    if(status /= nf90_noerr) call Handle_Error(status,'ave_var_open2')
    status = nf90_open(trim(path_forcing_averages)//"ff10m"//trim(suffix_forcing_averages),0,ncid(3))
    if(status /= nf90_noerr) call Handle_Error(status,'ave_var_open3')
    status = nf90_open(trim(path_forcing_averages)//"tskin"//trim(suffix_forcing_averages),0,ncid(4))
    if(status /= nf90_noerr) call Handle_Error(status,'ave_var_open4')
    status = nf90_open(trim(path_forcing_averages)//"evap"//trim(suffix_forcing_averages),0,ncid(5))
    if(status /= nf90_noerr) call Handle_Error(status,'ave_var_open5')
    status = nf90_open(trim(path_forcing_averages)//"sndiv"//trim(suffix_forcing_averages),0,ncid(6))
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
        count=(/Nlon,Nlat,1,1/))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_avemelt')
    status  = nf90_get_var(ncid(2),ID(2),AveAcc,start=(/1,1,1,1/), &
        count=(/Nlon,Nlat,1,1/))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_aveacc')
    status  = nf90_get_var(ncid(3),ID(3),AveWind,start=(/1,1,1,1/), &
        count=(/Nlon,Nlat,1,1/))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_avewind')
    status  = nf90_get_var(ncid(4),ID(4),AveTsurf,start=(/1,1,1,1/), &
        count=(/Nlon,Nlat,1,1/))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_avetsurf')
    status  = nf90_get_var(ncid(5),ID(5),AveSubl,start=(/1,1,1,1/), &
        count=(/Nlon,Nlat,1,1/))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_avesubl')
    status  = nf90_get_var(ncid(6),ID(6),AveSnowDrif,start=(/1,1,1,1/), &
        count=(/Nlon,Nlat,1,1/))
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
    do i = 1, Nlon
        do j = 1, Nlat
            AveAcc(i,j) = (AveAcc(i,j)+AveSubl(i,j)-AveSnowDrif(i,j)) &
                * (seconds_per_year)
            AveMelt(i,j) = AveMelt(i,j) * (seconds_per_year)
        end do
    end do
    
end subroutine Load_Ave_Forcing


! *******************************************************


subroutine Load_TimeSeries_Forcing(SnowMelt, PreTot, PreSol,PreLiq, Sublim, SnowDrif, TempSurf, FF10m, Nt_forcing, &
    ind_lon, ind_lat, dtobs, Nlon_timeseries)
    
    integer :: status, ind_t, ncid(50), Nt_forcing, ind_lon, ind_lat, ID(50), dtobs, Nlon_timeseries
    double precision :: remove_Psol,remove_Pliq,remove_Ptot,remove_Melt
    double precision, dimension(Nt_forcing) :: SnowMelt,PreTot,PreSol,PreLiq,Sublim,TempSurf, &
        SnowDrif,FF10m

    integer :: latfile, lonfile, fnumb_i
    character*255 :: add, fnumb

    lonfile = mod(ind_lon,Nlon_timeseries)
    if (lonfile==0) lonfile = Nlon_timeseries
    fnumb_i = floor((real(ind_lon)-0.001)/real(Nlon_timeseries))+1
    if (fnumb_i<=9) write(fnumb,'(I1)') fnumb_i
    if (fnumb_i>=10) write(fnumb,'(I2)') fnumb_i

    latfile = ind_lat

    add = trim(prefix_forcing_timeseries)//trim(fnumb)//".nc"
    
    print *, "ind_lat: ", ind_lat
    print *, "latfile: ", latfile
    print *, "ind_lon: ", ind_lon
    print *, "lonfile: ", lonfile
    print *, "fnumb: ", fnumb
    print *, "Nt_forcing: ", Nt_forcing
    print *, "Looking for timeseries at: "
    print *, trim(path_forcing_timeseries)//"VAR"//trim(add)
    print *, " "

    ! Open the snowmelt netCDF file
    status = nf90_open(trim(path_forcing_timeseries)//"snowmelt"//trim(add),0,ncid(1))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_ts_open1')
    status = nf90_open(trim(path_forcing_timeseries)//"precip"//trim(add),0,ncid(2))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_ts_open2')
    status = nf90_open(trim(path_forcing_timeseries)//"snowfall"//trim(add),0,ncid(3))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_ts_open3')
    status = nf90_open(trim(path_forcing_timeseries)//"evap"//trim(add),0,ncid(4))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_ts_open4')
    status = nf90_open(trim(path_forcing_timeseries)//"tskin"//trim(add),0,ncid(5))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_ts_open5')
    status = nf90_open(trim(path_forcing_timeseries)//"sndiv"//trim(add),0,ncid(6))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_ts_open6')
    status = nf90_open(trim(path_forcing_timeseries)//"ff10m"//trim(add),0,ncid(7))
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
        count=(/1,1,1,Nt_forcing/))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var1')
    print *, "Read snowmelt..."
    status  = nf90_get_var(ncid(2),ID(2),PreTot,start=(/lonfile,latfile,1,1/), &
        count=(/1,1,1,Nt_forcing/))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var2')
    print *, "Read precipitation..."
    status  = nf90_get_var(ncid(3),ID(3),PreSol,start=(/lonfile,latfile,1,1/), &
        count=(/1,1,1,Nt_forcing/))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var3')
    print *, "Read snowfall..."
    status  = nf90_get_var(ncid(4),ID(4),Sublim,start=(/lonfile,latfile,1,1/), &
        count=(/1,1,1,Nt_forcing/))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var4')
    print *, "Read sublimation..."
    status  = nf90_get_var(ncid(5),ID(5),TempSurf,start=(/lonfile,latfile,1,1/), &
        count=(/1,1,1,Nt_forcing/))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var5')
    print *, "Read skin temperature..."
    status  = nf90_get_var(ncid(6),ID(6),SnowDrif,start=(/lonfile,latfile,1,1/), &
        count=(/1,1,1,Nt_forcing/))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var6')
    print *, "Read snow drift..."
    status  = nf90_get_var(ncid(7),ID(7),FF10m,start=(/lonfile,latfile,1,1/), &
        count=(/1,1,1,Nt_forcing/))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var7')
    print *, "Read wind speed..."
    print *, ' '

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
    do ind_t = 1, Nt_forcing
        SnowMelt(ind_t) = SnowMelt(ind_t) * dtobs
        PreTot(ind_t) = PreTot(ind_t) * dtobs
        PreSol(ind_t) = PreSol(ind_t) * dtobs
        Sublim(ind_t) = Sublim(ind_t) * dtobs
        SnowDrif(ind_t) = SnowDrif(ind_t) * dtobs
        PreLiq(ind_t) = PreLiq(ind_t) * dtobs
    end do

    do ind_t = 1, Nt_forcing
        if (SnowMelt(ind_t) < ts_minimum) then
            remove_Melt = remove_Melt + SnowMelt(ind_t)
            SnowMelt(ind_t) = 0.     
        endif
        
        if (PreSol(ind_t) < ts_minimum) then
            remove_Psol = remove_Psol + PreSol(ind_t)
            PreSol(ind_t) = 0.
        endif
    
        if (PreTot(ind_t) < ts_minimum) then
            remove_Ptot = remove_Ptot + PreTot(ind_t)
            PreTot(ind_t) = 0.
        endif
    
        if (TempSurf(ind_t) > 267.) then
            PreLiq(ind_t) = PreTot(ind_t)-PreSol(ind_t)
            if (PreLiq(ind_t) < ts_minimum) then
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


subroutine Restart_From_Spinup(ind_z_max, ind_z_surf, Rho, M, T, Depth, Mlwc, DZ, DenRho, Refreeze)
        
    integer :: ind_z_max, ind_z_surf
    integer :: ind_z, status, ncid(50), ID(50), LayerID
    
    double precision, dimension(ind_z_max) :: Rho, M, T, Depth, Mlwc, DZ, DenRho, Refreeze
    
    character*255 :: pad
    
    pad = trim(path_restart)//trim(fname_restart_from_spinup)
    
    print *, "Path of spinup restart file: "
    print *, trim(pad)
    print *, ' '

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

subroutine Restart_From_Run(prev_nt, ind_z_max, ind_z_surf, Rho, M, T, Depth, Mlwc, DZ, Year, DenRho, Refreeze, username, &
                                point_numb, prefix_output, project_name)
        
    integer :: ind_z_max, ind_z_surf, prev_nt
    integer :: ind_z, status, ncid(50), ID(50), LayerID(2)
    
    double precision, dimension(ind_z_max) :: Rho, M, T, Depth, Mlwc, DZ, Year, DenRho, Refreeze
    
    character*255 :: pad, username, point_numb, project_name, prefix_output
    
    pad = trim(path_restart)//trim(fname_restart_from_previous_run)
    
    print *, 'Path of run restart files: :'
    print *, trim(pad)
    print *, ' '

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
        print *, 'netCDF error (',msg,'): ', nf90_strerror(stat)
        stop
    endif

end subroutine Handle_Error


end module openNetCDF
