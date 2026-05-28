module initialise_model
    !*** Subroutines for initialising the firn column before the spin-up starts ***!
    use model_settings
    implicit none
    private

    public :: Find_Grid, Interpol_Forcing, Index_Ave_Forcing, Init_Density_Prof, Init_Temp_Prof
    
contains


! *******************************************************


subroutine Find_Grid(ind_lon, ind_lat, lon_current, lat_current, Latitude, Longitude, LSM)
    ! TKTKTK: should be removed since lat/lon already found in distributor.f90

    !*** Determine the indices of the current run point on the forcing grid ***!
    
    ! declare arguments

    integer, intent(out) :: ind_lon, ind_lat
    double precision, intent(in) :: lon_current, lat_current
    double precision, dimension(config%forcing_dimensions%Nlon,config%forcing_dimensions%Nlat), intent(in) :: Latitude, Longitude, LSM

    ! declare local variables
    integer :: i, j
    double precision :: dmax, dist, lon_grid, lat_grid

    dmax = 999.
    do i = 1, config%forcing_dimensions%Nlon
        do j = 1, config%forcing_dimensions%Nlat
            dist = sqrt(((lon_current-Longitude(i,j))*cos(((lat_current+Latitude(i,j))/2.)/360.* &
                asin(1.)*4.))**2. + (lat_current-Latitude(i,j))**2.)
            if ((dist < dmax) .and. (LSM(i,j) > 0.5)) then
                dmax = dist
                ind_lon = i
                ind_lat = j
                lon_grid = Longitude(i,j)
                lat_grid = Latitude(i,j)    
            endif
        enddo
    enddo

    write(log_unit, *) " "
    write(log_unit, *) " Input Lon, Lat: ", lon_current, ", ", lat_current
    write(log_unit, *) "------------------------------------"
    write(log_unit, *) " Closest gridpoint: ", ind_lon, ", ", ind_lat
    write(log_unit, *) "     with Lon, Lat: ", lon_grid, ", ", lat_grid
    write(log_unit, *) "-----------------------------------------------------------------"
    write(log_unit, *) " "

end subroutine Find_Grid


! *******************************************************


subroutine Interpol_Forcing(TempSurf, PreSol, PreLiq, Sublim, SnowMelt, SnowDrif, FF10m, &
    TempFM, PsolFM, PliqFM, SublFM, MeltFM, DrifFM, Rho0FM, Nt_forcing, Nt_model_interpol, &
    Nt_model_tot, dtmodel)
    !*** Linearly interpolate RACMO forcing to model time step and calculate fresh snow density ***!

    ! declare arguments
    integer, intent(in) :: Nt_forcing, Nt_model_interpol, Nt_model_tot, dtmodel
    double precision, dimension(:), intent(in)  :: TempSurf, PreSol, PreLiq, Sublim, SnowMelt, SnowDrif, FF10m
    double precision, dimension(:), intent(out) :: TempFM, PsolFM, PliqFM, SublFM, MeltFM, DrifFM, Rho0FM

    ! declare local variables
    integer :: step, a, b, numSnow
    double precision :: part1, part2, TempSnow, ff10Snow
    double precision, allocatable, dimension(:) :: ff10FM

    allocate(ff10FM(Nt_model_tot))
    
    write(log_unit, *) 'Forcing before interpolation:'
    write(log_unit, *) 'TempSurf(1:10)'
    write(log_unit, *) TempSurf(1:10)
    write(log_unit, *) ' '
    write(log_unit, *) 'PreSol(1:10)'
    write(log_unit, *) PreSol(1:10)
    write(log_unit, *) ' '
    write(log_unit, *) 'PreLiq(1:10)'
    write(log_unit, *) PreLiq(1:10)
    write(log_unit, *) ' '
    write(log_unit, *) 'Sublim(1:10)'
    write(log_unit, *) Sublim(1:10)
    write(log_unit, *) ' '
    write(log_unit, *) 'SnowMelt(1:10)'
    write(log_unit, *) SnowMelt(1:10)
    write(log_unit, *) ' '
    write(log_unit, *) 'SnowDrif(1:10)'
    write(log_unit, *) SnowDrif(1:10)
    write(log_unit, *) ' '
    write(log_unit, *) 'FF10m(1:10)'
    write(log_unit, *) FF10m(1:10)
    write(log_unit, *) ' '

    do a = 1, (Nt_forcing-1)
        do  b = 1, Nt_model_interpol
            step = (a-1)*Nt_model_interpol + b
            part1 = (b-1.)/Nt_model_interpol
            part2 = 1. - part1
            TempFM(step) = part2*TempSurf(a) + part1*TempSurf(a+1)
            PSolFM(step) = (part2*PreSol(a) + part1*PreSol(a+1))/Nt_model_interpol
            PliqFM(step) = (part2*PreLiq(a) + part1*PreLiq(a+1))/Nt_model_interpol
            SublFM(step) = (part2*Sublim(a) + part1*Sublim(a+1))/Nt_model_interpol
            MeltFM(step) = (part2*SnowMelt(a) + part1*SnowMelt(a+1))/Nt_model_interpol
            DrifFM(step) = (part2*SnowDrif(a) + part1*SnowDrif(a+1))/Nt_model_interpol
            ff10FM(step) = part2*FF10m(a) + part1*FF10m(a+1)
        end do
    end do

    do b = 1, Nt_model_interpol
        step = (Nt_forcing-1)*Nt_model_interpol + b
        TempFM(step) = TempSurf(Nt_forcing)
        PSolFM(step) = PreSol(Nt_forcing)/Nt_model_interpol
        PliqFM(step) = PreLiq(Nt_forcing)/Nt_model_interpol
        SublFM(step) = Sublim(Nt_forcing)/Nt_model_interpol
        MeltFM(step) = SnowMelt(Nt_forcing)/Nt_model_interpol
        DrifFM(step) = SnowDrif(Nt_forcing)/Nt_model_interpol
        ff10FM(step) = FF10m(Nt_forcing)
    end do

    if (trim(domain) == "ANT27") then 
        numSnow = 1
    else
        numSnow = max(int(config%model_choices%dtSnow/dtmodel),1)
    end if 
    
    write(log_unit, *) 'numSnow: ', numSnow
    write(log_unit, *) ' '

    if (numSnow == 1) then
        ! TKTKTK: make choice of surface snow explicit in model.toml
        
        ! Use current temperature and wind speed for snow parameterisations
        if (trim(domain) == "FGRN11" .or. trim(domain) == "FGRN055" .or. trim(domain) == "FGRN055_era055") then
            do step = 1, Nt_model_tot
                Rho0FM(step) = 362.1 + 2.78*(TempFM(step) - const%Tmelt) ! from Fausto et al. 2018 
            end do
        else
            do step = 1, Nt_model_tot
                Rho0FM(step) = 82.97 + 0.769*TempFM(step) + 11.67*ff10FM(step)
                Rho0FM(step) = MIN(470.,Rho0FM(step))
            end do
        end if
    else
        ! Use mean temperature and wind speed, averaged over numSnow time steps
        if (trim(domain) == "FGRN11" .or. trim(domain) == "FGRN055" .or. trim(domain) == "FGRN055_era055") then
            ! Greenland
            TempSnow = sum( TempFM(1:numSnow) )/numSnow
            do step = 1, numSnow
                Rho0FM(step) = 362.1 + 2.78*(TempSnow - const%Tmelt)
            end do
            do step = (numSnow+1), Nt_model_tot
                TempSnow = sum( TempFM(step-numSnow:step) )/numSnow
                Rho0FM(step) = 362.1 + 2.78*(TempSnow - const%Tmelt)
            end do
        else
            ! Not Greenland
            TempSnow = sum( TempFM(1:numSnow) )/numSnow
            ff10Snow = sum( ff10FM(1:numSnow) )/numSnow
            do step = 1, numSnow
                Rho0FM(step) = 82.97 + 0.769*TempSnow + 11.67*ff10Snow
                Rho0FM(step) = MIN(470.,Rho0FM(step))
            end do
            do step = (numSnow+1), Nt_model_tot
                TempSnow = sum( TempFM(step-numSnow:step) )/numSnow
                ff10Snow = sum( ff10FM(step-numSnow:step) )/numSnow
                Rho0FM(step) = 82.97 + 0.769*TempSnow + 11.67*ff10Snow
                Rho0FM(step) = MIN(470.,Rho0FM(step))
            end do
        end if
    end if

    write(log_unit, *) 'TempSnow: ', TempSnow
    write(log_unit, *) ' '
    write(log_unit, *) 'Snow density (after interpolation):'
    write(log_unit, *) 'Rho0FM(1:10)'
    write(log_unit, *) Rho0FM(1:10)
    write(log_unit, *) ' '

    deallocate(ff10FM)

end subroutine Interpol_Forcing


! *******************************************************


subroutine Index_Ave_Forcing(AveTsurf, AveAcc, AveWind, AveMelt, ISM, tsav, acav, ffav, IceShelf, ind_lon, ind_lat)
    !*** Index the avarage forcing fields for the current run ***! 

    ! declare arguments
    integer, intent(in) :: ind_lon, ind_lat
    integer, intent(out) :: IceShelf
    double precision, intent(out) :: tsav, acav, ffav
    double precision, dimension(config%forcing_dimensions%Nlon,config%forcing_dimensions%Nlat), intent(in) :: AveTsurf, AveAcc, AveWind, AveMelt, ISM

    acav = AveAcc(ind_lon, ind_lat)
    if (acav < 0) acav = 0.1
    ffav = AveWind(ind_lon, ind_lat)
    tsav = AveTsurf(ind_lon, ind_lat)
    if (tsav > const%Tmelt) tsav = const%Tmelt

    IceShelf = int(ISM(ind_lon, ind_lat))

    write(log_unit, *) "  "
    write(log_unit, *) "------------------------------------"    
    write(log_unit,'(A19,1X,F8.3,1X,A1)') " Average Tsurf:     ", tsav, "K"
    write(log_unit,'(A19,1X,F8.3,1X,A13)') " Average Acc:      ", acav, "mm w.e. yr-1"
    write(log_unit,'(A19,1X,F8.3,1X,A5)') " Average Wind:      ", ffav, "m s-1"
    write(log_unit,'(A19,1X,F8.3,1X,A13)') " Total Melt:       ", AveMelt(ind_lon, ind_lat), "mm w.e. yr-1"
    write(log_unit, *) "------------------------------------"
    write(log_unit, *) "  "

end subroutine Index_Ave_Forcing


! *******************************************************


subroutine Init_Density_Prof(rho0, acav, tsav, DZ, Rho, M)
    !*** Initialise the density profile ***!
    double precision, intent(in) :: rho0, acav, tsav
    double precision, dimension(ind_z_max), intent(inout) :: DZ, Rho
    double precision, dimension(ind_z_max), intent(out) :: M

    ! declare local variables
    integer :: ind_z
    double precision :: drho, part1
        
    Rho(ind_z_surf) = rho0 ! rho0 set by Rho0FM in Interpol_Forcing
    M(ind_z_surf) = Rho(ind_z_surf) * DZ(ind_z_surf)

    do ind_z = (ind_z_surf-1), 1, -1

        part1 = (const%rhoi-Rho(ind_z+1))*exp((-const%Ec/(const%R*tsav))+(const%Eg/(const%R*tsav)))

        if (Rho(ind_z+1) <= 550.) then
            drho = 0.07*config%model_choices%dzmax*Rho(ind_z+1)*const%g*part1
        else
            drho = 0.03*config%model_choices%dzmax*Rho(ind_z+1)*const%g*part1
        endif


        Rho(ind_z) = drho + Rho(ind_z+1)
        if (Rho(ind_z) > const%rhoi) Rho(ind_z) = const%rhoi

        M(ind_z) = Rho(ind_z) * DZ(ind_z)

    enddo
    
    write(log_unit, *) 'Initial density 10 lowermost layers:'
    write(log_unit, *) Rho(1:10)
    write(log_unit, *) " "

end subroutine Init_Density_Prof

 
! *******************************************************


subroutine Init_Temp_Prof(tsav, T, Rho, Depth)
    !*** Initialise the temperature profile ***!
    double precision, intent(in) :: tsav
    double precision, dimension(ind_z_max), intent(in) :: Depth, Rho
    double precision, dimension(ind_z_max), intent(out) :: T

    ! declare local variables
    integer :: ind_z
    double precision :: ki, ci, period, om, diff
    double precision :: ampts, kice, kice_ref, kcal, kf, kair, kair_ref, theta

    ampts = 10.

    ci = 152.5 + 7.122 * tsav                       ! heat capacity, Paterson (1994)
    om = 2.*const%pi/const%seconds_per_year                               ! TKTKTK: set in constants; rads per second

    do ind_z = ind_z_surf, 1, -1                               ! temperature - depth loop

        ! TKTKTK the choice of these physics should be made explicit in model.toml file

        !ki = 0.021+2.5*(Rho(k)/1000.)**2                       ! Anderson (1976), old method
        kice = 9.828 * exp(-0.0057*T(ind_z))                    ! Paterson et al., 1994
        kice_ref = 9.828 * exp(-0.0057*270.15)                  ! Paterson et al., 1994
        kcal = 0.024 - 1.23E-4*Rho(ind_z) + 2.5E-6*Rho(ind_z)**2.
        kf = 2.107 + 0.003618*(Rho(ind_z)-const%rhoi)                     ! Calonne (2019)
        kair = (2.334E-3*T(ind_z)**(3./2.))/(164.54 + T(ind_z))           ! Reid (1966) #updated by ehc on 11/06/25 so no longer integer division (3/2)
        kair_ref = (2.334E-3*270.15**(3./2.))/(164.54 + 270.15)
        theta = 1./(1.+exp(-0.04*(Rho(ind_z)-450.)))
        ki = (1.-theta)*kice/kice_ref*kair/kair_ref*kcal + theta*kice/kice_ref*kf    ! Calonne (2019)

        Diff = ki/(Rho(ind_z)*ci)                       ! thermal diffusivity
        if (ind_z == 1) then
            T(ind_z) = tsav
        else
            if (config%initialization%beginT == 1) then
                ! Winter
                T(ind_z) = tsav - ampts*exp(-(om/(2.*Diff))**0.5*Depth(ind_z))* &
                    cos(-(om/(2.*Diff))**0.5*Depth(ind_z))
            elseif (config%initialization%beginT == 2) then
                ! Summer
                T(ind_z) = tsav + ampts*exp(-(om/(2.*Diff))**0.5*Depth(ind_z))* &
                    cos(-(om/(2.*Diff))**0.5*Depth(ind_z))
            else 
                T(ind_z) = tsav
            endif
            if (T(ind_z) > 272.15) T(ind_z) = 272.15
        endif

    enddo

    write(log_unit, *) 'Initial temperature 10 lowermost layers:'
    write(log_unit, *) T(1:10)
    write(log_unit, *) " "
    
end subroutine Init_Temp_Prof
    

end module initialise_model
