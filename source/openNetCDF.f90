module openNetCDF
    !*** subroutines for loading forcing data

    use netcdf, only : nf90_open,nf90_inq_varid, nf90_inq_dimid, nf90_inquire_dimension, & 
        nf90_close,nf90_get_var,nf90_noerr, nf90_strerror

    implicit none
    private

    public :: Load_Ave_Forcing, Load_TimeSeries_Forcing, Handle_Error, Restart_From_Spinup, Restart_From_Run
    
contains


! *******************************************************


subroutine Load_Ave_Forcing(AveTsurf, AveAcc, AveWind, AveMelt, LSM,Nlat, Nlon, Nt_forcing, &
    nyears, Latitude, Longitude, ISM, username, domain)
    
    integer :: status,ncid(50),ID(50),Nlat,Nlon,Nt_forcing,i,j,nyears
    double precision, dimension(Nlon,Nlat) :: AveTsurf,AveAcc,AveWind,AveSubl, &
        AveSnowDrif,AveMelt,LSM,ISM,Latitude,Longitude, Icemask_GR
    
    character*255 :: add,pad,username,domain,path_dir,pad_mask

    path_dir = "/ec/res4/scratch/"
        
    if (domain == "ANT27") then    
        pad = ''//trim(path_dir)//''//trim(username)//"/backup_HPC/data/input/ant27_era_files/averages/"    
        add = "_ANT27_79-20_ave.nc"
    elseif (domain == "XPEN055") then
        pad = ''//trim(path_dir)//''//trim(username)//"/FM_Data/INPUT/XPEN055_averages/"    
        add = "_XPEN055_79-16_ave.nc"
    elseif (domain == "FGRN11") then
        pad = ''//trim(path_dir)//''//trim(username)//"/data/input/era_files/averages/"    
        add = "_FGRN11_60-79_ave.nc"
    elseif (domain == "FGRN055" .or. domain == "FGRN055_era055") then
        pad = ''//trim(path_dir)//''//trim(username)//"/FGRN055_era055/input/averages/"    
        add = "_FGRN055-era055_1960-1981_ave.nc"
    elseif (domain == "PAT055") then
        pad = ''//trim(path_dir)//''//trim(username)//"/FM_Data/INPUT/PAT055_averages/"    
        add = "_PAT055_79-12_ave.nc"
    elseif (domain == "XDML055") then
        pad = ''//trim(path_dir)//''//trim(username)//"/FM_Data/INPUT/XDML055_averages/"    
        add = "_XDML055_79-15_ave.nc"
    elseif (domain == "ASE055") then
        pad = ''//trim(path_dir)//''//trim(username)//"/FM_Data/INPUT/ASE055_averages/"    
        add = "_ASE055_79-15_ave.nc"
    elseif (domain == "DMIS055") then
        pad = ''//trim(path_dir)//''//trim(username)//"/FM_Data/INPUT/DMIS055_averages/"    
        add = "_DMIS055_79-17_ave.nc"
    else
        call Handle_Error(41,'no valid domain') 
     endif

    print *, "looking for averages at: ", pad
    
    if (domain == "ANT27") then

        status = nf90_open(trim(pad)//"/ANT27_Masks_ice_and_shelves.nc",0,ncid(1))
        !status = nf90_open(trim(pad)//"../lsm_ANT27.nc",0,ncid(1))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_open1')
        status = nf90_inq_varid(ncid(1),"IceMask",ID(1))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_inq_varid')    
        status = nf90_inq_varid(ncid(1),"lat",ID(11))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_inq_varid11')
        status = nf90_inq_varid(ncid(1),"lon",ID(12))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_inq_varid12')
        
        status  = nf90_get_var(ncid(1),ID(1),LSM,start=(/1,1,1,1/), &
            count=(/Nlon,Nlat,1,1/))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_lsm')
        status  = nf90_get_var(ncid(1),ID(11),Latitude)        
        if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_lat')
        status  = nf90_get_var(ncid(1),ID(12),Longitude)    
        if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_lon')

        ! Open Ice Shelf Mask

        status = nf90_open(trim(pad)//"/ANT27_Masks_ice_and_shelves.nc",0,ncid(6))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_open6')
        status = nf90_inq_varid(ncid(6),"IceShelve",ID(6))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_inq_varid')
        status  = nf90_get_var(ncid(6),ID(6),ISM,start=(/1,1/), &
            count=(/Nlon,Nlat/))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_ism')


        status = nf90_close(ncid(1))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_close1')
        status = nf90_close(ncid(6))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_close6')   
   
    elseif (domain == "XPEN055") then

        status = nf90_open(trim(pad)//"../Height_latlon_XPEN055.nc",0,ncid(1))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_open1')
        status = nf90_inq_varid(ncid(1),"mask2d",ID(1))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_inq_varid')    
        status = nf90_inq_varid(ncid(1),"lat",ID(11))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_inq_varid11')
        status = nf90_inq_varid(ncid(1),"lon",ID(12))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_inq_varid12')
        
        status  = nf90_get_var(ncid(1),ID(1),LSM,start=(/1,1/), &
            count=(/Nlon,Nlat/))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_lsm')
        status  = nf90_get_var(ncid(1),ID(11),Latitude)        
        if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_lat')
        status  = nf90_get_var(ncid(1),ID(12),Longitude)    
        if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_lon')

        ! Open Ice Shelf Mask

        status = nf90_inq_varid(ncid(1),"iceshelves",ID(6))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_inq_varid')
        status  = nf90_get_var(ncid(1),ID(6),ISM,start=(/1,1/), &
            count=(/Nlon,Nlat/))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_ism')

        status = nf90_close(ncid(1))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_close1')

    elseif (domain == "FGRN11") then

        status = nf90_open(trim(pad)//"../mask/FGRN11_Masks_wholedomain.nc",0,ncid(1))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_open1')
        status = nf90_inq_varid(ncid(1),"icemask",ID(1))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_inq_varid_lsm')    
        status = nf90_inq_varid(ncid(1),"lat",ID(11))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_inq_varid11')
        status = nf90_inq_varid(ncid(1),"lon",ID(12))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_inq_varid12')
        
        status  = nf90_get_var(ncid(1),ID(1),LSM,start=(/1,1/), &
            count=(/Nlon,Nlat/))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_lsm')
        status  = nf90_get_var(ncid(1),ID(11),Latitude)        
        if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_lat')
        status  = nf90_get_var(ncid(1),ID(12),Longitude)    
        if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_lon')

        ! No Ice Shelves in Greenland
        ISM(:,:) = 0

    elseif (domain == "FGRN055" .or. domain == "FGRN055_era055") then

        pad_mask = "/perm/"//trim(username)//"/code/IMAU-FDM/reference/"//trim(domain)//"/FGRN055_Masks.nc"
        
        print *, "Path to mask: ", pad_mask

        status = nf90_open(trim(pad_mask),0,ncid(1))
        
        if(status /= nf90_noerr) call Handle_Error(status,'mask_open1')
        status = nf90_inq_varid(ncid(1),"Icemask_GR",ID(1))
        if(status /= nf90_noerr) call Handle_Error(status,'mask_inq_varid_lsm')    
        status = nf90_inq_varid(ncid(1),"lat",ID(11))
        if(status /= nf90_noerr) call Handle_Error(status,'mask_inq_varid11')
        status = nf90_inq_varid(ncid(1),"lon",ID(12))
        if(status /= nf90_noerr) call Handle_Error(status,'mask_inq_varid12')
        
        status  = nf90_get_var(ncid(1),ID(1),Icemask_GR,start=(/1,1/), &
            count=(/Nlon,Nlat/))
        if(status /= nf90_noerr) call Handle_Error(status,'mask_get_var_lsm')
        status  = nf90_get_var(ncid(1),ID(11),Latitude)        
        if(status /= nf90_noerr) call Handle_Error(status,'mask_get_var_lat')
        status  = nf90_get_var(ncid(1),ID(12),Longitude)    
        if(status /= nf90_noerr) call Handle_Error(status,'mask_get_var_lon')

        ! No Ice Shelves in Greenland
        ISM(:,:) = 0

    elseif (domain == "PAT055") then

        status = nf90_open(trim(pad)//"../lsm_PAT055.nc",0,ncid(1))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_open1')
        status = nf90_inq_varid(ncid(1),"mask",ID(1))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_inq_varid')    
        status = nf90_inq_varid(ncid(1),"lat",ID(11))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_inq_varid11')
        status = nf90_inq_varid(ncid(1),"lon",ID(12))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_inq_varid12')
        
        status  = nf90_get_var(ncid(1),ID(1),LSM,start=(/1,1/), &
            count=(/Nlon,Nlat/))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_lsm')
        status  = nf90_get_var(ncid(1),ID(11),Latitude)        
        if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_lat')
        status  = nf90_get_var(ncid(1),ID(12),Longitude)    
        if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_lon')

        ! No Ice Shelves in Patagonia
        ISM(:,:) = 0

    elseif (domain == "XDML055") then

        status = nf90_open(trim(pad)//"../lsm_XDML055.nc",0,ncid(1))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_open1')
        status = nf90_inq_varid(ncid(1),"ism",ID(1))   !ice sheet mask....
        if(status /= nf90_noerr) call Handle_Error(status,'nf_inq_varid')    
        status = nf90_inq_varid(ncid(1),"lat",ID(11))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_inq_varid11')
        status = nf90_inq_varid(ncid(1),"lon",ID(12))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_inq_varid12')
        
        status  = nf90_get_var(ncid(1),ID(1),LSM,start=(/1,1/), &
            count=(/Nlon,Nlat/))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_lsm')
        status  = nf90_get_var(ncid(1),ID(11),Latitude)        
        if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_lat')
        status  = nf90_get_var(ncid(1),ID(12),Longitude)    
        if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_lon')

        ! No Ice Shelves file yet
        ISM(:,:) = 0

    elseif (domain == "ASE055") then

        status = nf90_open(trim(pad)//"../Masks_ASE055.nc",0,ncid(1))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_open1')
        status = nf90_inq_varid(ncid(1),"LSM",ID(1))   !ice sheet mask....
        if(status /= nf90_noerr) call Handle_Error(status,'nf_inq_varid')    
        status = nf90_inq_varid(ncid(1),"ISM",ID(2))   !ice shelf mask....
        if(status /= nf90_noerr) call Handle_Error(status,'nf_inq_varid')
        status = nf90_inq_varid(ncid(1),"lat",ID(11))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_inq_varid11')
        status = nf90_inq_varid(ncid(1),"lon",ID(12))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_inq_varid12')
        
        status  = nf90_get_var(ncid(1),ID(1),LSM,start=(/1,1/), &
            count=(/Nlon,Nlat/))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_lsm')
        status  = nf90_get_var(ncid(1),ID(11),Latitude)        
        if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_lat')
        status  = nf90_get_var(ncid(1),ID(12),Longitude)    
        if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_lon')

        status  = nf90_get_var(ncid(1),ID(2),ISM,start=(/1,1/), &
            count=(/Nlon,Nlat/))

    elseif (domain == "DMIS055") then

        status = nf90_open(trim(pad)//"../Masks_DMIS055.nc",0,ncid(1))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_open1')
        status = nf90_inq_varid(ncid(1),"LSM",ID(1))   !ice sheet mask....
        if(status /= nf90_noerr) call Handle_Error(status,'nf_inq_varid')    
        status = nf90_inq_varid(ncid(1),"ISM",ID(2))   !ice shelf mask....
        if(status /= nf90_noerr) call Handle_Error(status,'nf_inq_varid')
        status = nf90_inq_varid(ncid(1),"lat",ID(11))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_inq_varid11')
        status = nf90_inq_varid(ncid(1),"lon",ID(12))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_inq_varid12')
        
        status  = nf90_get_var(ncid(1),ID(1),LSM,start=(/1,1/), &
            count=(/Nlon,Nlat/))
        if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_lsm')
        status  = nf90_get_var(ncid(1),ID(11),Latitude)        
        if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_lat')
        status  = nf90_get_var(ncid(1),ID(12),Longitude)    
        if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_lon')

        status  = nf90_get_var(ncid(1),ID(2),ISM,start=(/1,1/), &
            count=(/Nlon,Nlat/))

    else

         call Handle_Error(42,'no valid domain')

    endif 

    status = nf90_open(trim(pad)//"snowmelt"//trim(add),0,ncid(1))
    if(status /= nf90_noerr) call Handle_Error(status,'ave_var_open1')
    status = nf90_open(trim(pad)//"precip"//trim(add),0,ncid(2))
    if(status /= nf90_noerr) call Handle_Error(status,'ave_var_open2')
    status = nf90_open(trim(pad)//"ff10m"//trim(add),0,ncid(3))
    if(status /= nf90_noerr) call Handle_Error(status,'ave_var_open3')
    status = nf90_open(trim(pad)//"tskin"//trim(add),0,ncid(4))
    if(status /= nf90_noerr) call Handle_Error(status,'ave_var_open4')
    status = nf90_open(trim(pad)//"evap"//trim(add),0,ncid(5))
    if(status /= nf90_noerr) call Handle_Error(status,'ave_var_open5')
    status = nf90_open(trim(pad)//"sndiv"//trim(add),0,ncid(7))
    if(status /= nf90_noerr) call Handle_Error(status,'ave_var_open7')

    status = nf90_inq_varid(ncid(1),"snowmelt",ID(1))
    if(status /= nf90_noerr) call Handle_Error(status,'ave_var_inq_varid')
    status = nf90_inq_varid(ncid(2),"precip",ID(2))
    if(status /= nf90_noerr) call Handle_Error(status,'ave_var_inq_varid')
    status = nf90_inq_varid(ncid(3),"ff10m",ID(3))
    if(status /= nf90_noerr) call Handle_Error(status,'ave_var_inq_varid')
    status = nf90_inq_varid(ncid(4),"tskin",ID(4))
    if(status /= nf90_noerr) call Handle_Error(status,'ave_var_inq_varid')
    status = nf90_inq_varid(ncid(5),"evap",ID(5))
    if(status /= nf90_noerr) call Handle_Error(status,'ave_var_inq_varid')
    status = nf90_inq_varid(ncid(7),"sndiv",ID(7))
    if(status /= nf90_noerr) call Handle_Error(status,'ave_var_inq_varid')

    status  = nf90_get_var(ncid(1),ID(1),AveMelt,start=(/1,1,1,1/), &
        count=(/Nlon,Nlat,1,1/))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_avemelt')
    status  = nf90_get_var(ncid(2),ID(2),AveAcc,start=(/1,1,1,1/), &
        count=(/Nlon,Nlat,1,1/))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_get_var_aveacc')
    status  = nf90_get_var(ncid(3),ID(3),AveWind,start=(/1,1,1,1/), &
        count=(/Nlon,Nlat,1,1/))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_get_varavewind')
    status  = nf90_get_var(ncid(4),ID(4),AveTsurf,start=(/1,1,1,1/), &
        count=(/Nlon,Nlat,1,1/))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_get_varavetsurf')
    status  = nf90_get_var(ncid(5),ID(5),AveSubl,start=(/1,1,1,1/), &
        count=(/Nlon,Nlat,1,1/))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_get_varavesubl')
    status  = nf90_get_var(ncid(7),ID(7),AveSnowDrif,start=(/1,1,1,1/), &
        count=(/Nlon,Nlat,1,1/))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_get_varavesndiv')

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
    status = nf90_close(ncid(7))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_close7')

    ! Convert units from [mm w.e./s] to [mm w.e./yr]
    do i = 1, Nlon
        do j = 1, Nlat
            AveAcc(i,j) = (AveAcc(i,j)+AveSubl(i,j)-AveSnowDrif(i,j)) &
                * (365.*24.*3600.)
            AveMelt(i,j) = AveMelt(i,j) * (365.*24.*3600.)
        end do
    end do
    
end subroutine Load_Ave_Forcing


! *******************************************************


subroutine Load_TimeSeries_Forcing(SnowMelt, PreTot, PreSol,PreLiq, Sublim, SnowDrif, TempSurf, FF10m, Nt_forcing, &
    ind_lon, ind_lat, username, domain, dtobs)
    
    integer :: status, ind_t, ncid(50), Nt_forcing, ind_lon, ind_lat, ID(50), dtobs
    double precision :: remove_Psol,remove_Pliq,remove_Ptot,remove_Melt
    double precision, dimension(Nt_forcing) :: SnowMelt,PreTot,PreSol,PreLiq,Sublim,TempSurf, &
        SnowDrif,FF10m

    integer :: latfile,lonfile,fnumb_i
    character*255 :: add,pad,fnumb,username,domain
    
    if (domain == "ANT27") then

        !latfile = mod(ind_lat,15)
        !if (latfile==0) latfile = 15
        !fnumb_i = floor((real(ind_lat)-0.001)/15.)+1
        !if (fnumb_i<=9) write(fnumb,'(I1)') fnumb_i
        !if (fnumb_i>=10) write(fnumb,'(I2)') fnumb_i

        !lonfile = ind_lon
         !from here
        lonfile = mod(ind_lon,15)
        if (lonfile.eq.0) lonfile = 15
        fnumb_i = floor((real(ind_lon)-0.001)/15.)+1
        if (fnumb_i.le.9) write(fnumb,'(I1)') fnumb_i
        if (fnumb_i.ge.10) write(fnumb,'(I2)') fnumb_i

        latfile = ind_lat

        pad = "/ec/res4/scratch/"//trim(username)//"/backup_HPC/data/input/ant27_era_files/files_2023/"    
        add = "_ANT27_79-23_p"//trim(fnumb)//".nc"

    elseif (domain == "XPEN055") then

        lonfile = mod(ind_lon,5)
        if (lonfile==0) lonfile = 5
        fnumb_i = floor((real(ind_lon)-0.001)/5.)+1
        if (fnumb_i<=9) write(fnumb,'(I1)') fnumb_i
        if (fnumb_i>=10) write(fnumb,'(I2)') fnumb_i

        latfile = ind_lat
        
        pad = "/ec/res4/scratch/"//trim(username)//"/FM_Data/INPUT/XPEN055_files/"    
        add = "_XPEN055_79-16_p"//trim(fnumb)//".nc"

    elseif (domain == "FGRN11") then

        lonfile = mod(ind_lon,5)
        if (lonfile==0) lonfile = 5
        fnumb_i = floor((real(ind_lon)-0.001)/5.)+1
        if (fnumb_i<=9) write(fnumb,'(I1)') fnumb_i
        if (fnumb_i>=10) write(fnumb,'(I2)') fnumb_i

        latfile = ind_lat

        pad = "/ec/res4/scratch/"//trim(username)//"/data/input/era_files/files/"    
        add = "_FGRN11_60-16_p"//trim(fnumb)//".nc"

    elseif (domain == "FGRN055" .OR. domain == "FGRN055_era055") then

        lonfile = mod(ind_lon,6)
        if (lonfile==0) lonfile = 6
        fnumb_i = floor((real(ind_lon)-0.001)/6.)+1
        if (fnumb_i<=9) write(fnumb,'(I1)') fnumb_i
        if (fnumb_i>=10) write(fnumb,'(I2)') fnumb_i

        latfile = ind_lat

        pad = "/ec/res4/scratch/"//trim(username)//"/FGRN055_era055/input/timeseries/"    
        add = "_FGRN055-era055_1957-2023_p"//trim(fnumb)//".nc"

    elseif (domain == "PAT055") then

        latfile = mod(ind_lat,4)
        if (latfile==0) latfile = 4
        fnumb_i = floor((real(ind_lat)-0.001)/4.)+1
        if (fnumb_i<=9) write(fnumb,'(I1)') fnumb_i
        if (fnumb_i>=10) write(fnumb,'(I2)') fnumb_i

        lonfile = ind_lon

        pad = "/ec/res4/scratch/"//trim(username)//"/FM_Data/INPUT/PAT055_files/"    
        add = "_PAT055_79-12_p"//trim(fnumb)//".nc"

    elseif (domain == "XDML055") then

        latfile = mod(ind_lat,5)
        if (latfile==0) latfile = 5
        fnumb_i = floor((real(ind_lat)-0.001)/5.)+1
        if (fnumb_i<=9) write(fnumb,'(I1)') fnumb_i
        if (fnumb_i>=10) write(fnumb,'(I2)') fnumb_i

        lonfile = ind_lon

        pad = "/ec/res4/scratch/"//trim(username)//"/FM_Data/INPUT/XDML055_files/"    
        add = "_XDML055_79-15_p"//trim(fnumb)//".nc"

    elseif (domain == "ASE055") then

        latfile = mod(ind_lat,5)
        if (latfile==0) latfile = 5
        fnumb_i = floor((real(ind_lat)-0.001)/5.)+1
        if (fnumb_i<=9) write(fnumb,'(I1)') fnumb_i
        if (fnumb_i>=10) write(fnumb,'(I2)') fnumb_i

        lonfile = ind_lon

        pad = "/ec/res4/scratch/"//trim(username)//"/FM_Data/INPUT/ASE055_files/"    
        add = "_ASE055_79-15_p"//trim(fnumb)//".nc"

    elseif (domain == "DMIS055") then

        lonfile = mod(ind_lon,6)
        if (lonfile==0) lonfile = 6
        fnumb_i = floor((real(ind_lon)-0.001)/6.)+1
        if (fnumb_i<=9) write(fnumb,'(I1)') fnumb_i
        if (fnumb_i>=10) write(fnumb,'(I2)') fnumb_i

        latfile = ind_lat

        pad = "/ec/res4/scratch/"//trim(username)//"/FM_Data/INPUT/DMIS055_files/"    
        add = "_DMIS055_79-17_p"//trim(fnumb)//".nc"

    else
        call Handle_Error(43,'no valid domain') 
    endif
    
    print *, 'ind_lat, latfile, ind_lon, lonfile, fnumb, Nt_forcing'
    print *, ind_lat, latfile, ind_lon, lonfile, trim(fnumb), Nt_forcing
    print *, 'trim(pad), trim(add)'
    print *, trim(pad), trim(add)
    print *, ' '

    ! Open the snowmelt netCDF file
    status = nf90_open(trim(pad)//"snowmelt"//trim(add),0,ncid(1))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_open1')
    status = nf90_open(trim(pad)//"precip"//trim(add),0,ncid(2))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_open2')
    status = nf90_open(trim(pad)//"snowfall"//trim(add),0,ncid(3))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_open3')
    status = nf90_open(trim(pad)//"evap"//trim(add),0,ncid(4))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_open4')
    status = nf90_open(trim(pad)//"tskin"//trim(add),0,ncid(5))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_open5')
    status = nf90_open(trim(pad)//"sndiv"//trim(add),0,ncid(6))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_open6')
    status = nf90_open(trim(pad)//"ff10m"//trim(add),0,ncid(7))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_open7')
    
    !Get ID of the variables in the NetCDF files
    status = nf90_inq_varid(ncid(1),"snowmelt",ID(1))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_inq_varid1')
    status = nf90_inq_varid(ncid(2),"precip",ID(2))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_inq_varid2')
    status = nf90_inq_varid(ncid(3),"snowfall",ID(3))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_inq_varid3')
    status = nf90_inq_varid(ncid(4),"evap",ID(4))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_inq_varid4')
    status = nf90_inq_varid(ncid(5),"tskin",ID(5))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_inq_varid5')
    status = nf90_inq_varid(ncid(6),"sndiv",ID(6))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_inq_varid6')
    status = nf90_inq_varid(ncid(7),"ff10m",ID(7))
    if(status /= nf90_noerr) call Handle_Error(status,'nf_inq_varid7')

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

    do ind_t = 1, Nt_forcing
        SnowMelt(ind_t) = SnowMelt(ind_t) * dtobs
        PreTot(ind_t) = PreTot(ind_t) * dtobs
        PreSol(ind_t) = PreSol(ind_t) * dtobs
        Sublim(ind_t) = Sublim(ind_t) * dtobs
        SnowDrif(ind_t) = SnowDrif(ind_t) * dtobs
        PreLiq(ind_t) = PreLiq(ind_t) * dtobs
    end do

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


subroutine Restart_From_Spinup(ind_z_max, ind_z_surf, Rho, M, T, Depth, Mlwc, DZ, DenRho, Refreeze, username, &
                                domain, point_numb, fname_p1, project_name)
        
    integer :: ind_z_max, ind_z_surf
    integer :: ind_z, status, ncid(50), ID(50), LayerID
    
    double precision, dimension(ind_z_max) :: Rho, M, T, Depth, Mlwc, DZ, DenRho, Refreeze
    
    character*255 :: fname, pad, username, domain, ini_fname, point_numb, project_name, fname_p1
    
    pad = "/ec/res4/scratch/"//trim(username)//"/restart/"//trim(project_name)//"/"//trim(fname_p1)//"_restart_from_spinup_"//trim(point_numb)//".nc"
    
    print *, 'pad for restart file:'
    print *, trim(pad)
    print *, ' '

    ! Open the snowmelt netCDF file
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
                                domain, point_numb, fname_p1, project_name)
        
    integer :: ind_z_max, ind_z_surf, prev_nt
    integer :: ind_z, status, ncid(50), ID(50), LayerID(2)
    
    double precision, dimension(ind_z_max) :: Rho, M, T, Depth, Mlwc, DZ, Year, DenRho, Refreeze
    
    character*255 :: fname, pad, username, domain, ini_fname, point_numb, project_name, fname_p1
    
    pad = "/ec/res4/scratch/"//trim(username)//"/restart/"//trim(project_name)//"/"//trim(fname_p1)//&
    "_restart_from_2023_run_"//trim(point_numb)//".nc"
    
    print *, 'pad for restart file:'
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
