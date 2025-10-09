module firn_physics
    !*** All subroutines and functions for calculating changes in density and temperature
    
    use water_physics, only: Bucket_Method, LWrefreeze
    use model_settings

    implicit none
    private

    public :: Update_Surface, Solve_Temp_Exp, Solve_Temp_Imp, Densific, Calc_Integrated_Var
    
contains


! *******************************************************


subroutine Update_Surface(ind_z_max, ind_z_surf, dtmodel, rho0, rhoi, acav, Lh, h_surf, vice, vmelt, vacc, vsub, &
    vsnd, vfc, vbouy, Ts, PSol, PLiq, Su, Me, Sd, M, T, DZ, Rho, DenRho, Mlwc,Refreeze, &
    ImpExp, IceShelf, Msurfmelt, Mrain, Msolin, Mrunoff, Mrefreeze)
    !*** Add and remove mass to the surface layer and calculate the velocity components ***!

    ! declare arguments
    integer, intent(in) :: ind_z_max, dtmodel, ImpExp, IceShelf
    integer, intent(inout) :: ind_z_surf
    double precision, intent(in) :: rho0, rhoi, acav, Lh
    double precision, intent(inout) :: Ts, Psol, Pliq, Su, Sd, Me
    double precision, intent(inout) :: Msurfmelt, Mrain, Msolin, Mrunoff, Mrefreeze
    double precision, intent(inout) :: vmelt, vice, h_surf, vacc, vsub, vsnd, vfc, vbouy
    double precision,dimension(ind_z_max), intent(inout) :: Rho, M, T, DZ, Mlwc, Refreeze, DenRho

    ! declare local variables
    integer :: ind_z
    double precision :: mdiff, macc, mice, mrun, Mmelt, oldDZ

    ! Prepare climate input for time step
    if (ImpExp == 1) then
        if (Me > 1e-06) then
            Mmelt = Me + Pliq
            Msurfmelt = Msurfmelt + Me
            Mrain = Mrain + Pliq
        elseif (Ts > Tmelt) then
            Mmelt = Me + Pliq
            Msurfmelt = Msurfmelt + Me
            Mrain = Mrain + Pliq
        else
            Psol = Psol + Pliq
            Pliq = 0.
            Mmelt = 0.
            Me = 0.
        endif
    else
        Psol = Psol + Pliq
        Mmelt = 0.
        Pliq = 0.
        Me = 0.
    endif

    ! Calculate mass change in ind_z_surf and velocity components of the surface

    ! Firn compaction
    vfc = 0.
    do ind_z = 1, ind_z_surf
        oldDZ = DZ(ind_z)
        DZ(ind_z) = M(ind_z)/Rho(ind_z)         ! Update layer thickness for new density
        vfc = vfc - (oldDZ - DZ(ind_z))
        DenRho(ind_z) = DenRho(ind_z) - (oldDZ - DZ(ind_z))   
    enddo

    ! Add mass to the upper layer
    M(ind_z_surf) = M(ind_z_surf) + (Psol+Su-Sd)
    Msolin = Msolin + (Psol+Su-Sd)

    ! Snow accumulation
    DZ(ind_z_surf) = DZ(ind_z_surf) + (Psol/rho0)
    vacc = Psol/rho0

    ! Surface sublimation
    DZ(ind_z_surf) = DZ(ind_z_surf) + (Su/Rho(ind_z_surf))
    vsub = Su/Rho(ind_z_surf)

    ! Snow drift
    if (Sd > 0.) then
        DZ(ind_z_surf) = DZ(ind_z_surf) - (Sd/Rho(ind_z_surf))
        vsnd = -1. * (Sd/Rho(ind_z_surf))
    else
        DZ(ind_z_surf) = DZ(ind_z_surf) - (Sd/rho0)
        vsnd = -1. * (Sd/rho0)
    endif

    ! Recalculate density of upper layer
    Rho(ind_z_surf) = M(ind_z_surf)/DZ(ind_z_surf)

    ! Melt
    vmelt = -1 * Me/Rho(ind_z_surf)

    ! Vertical (downward) velocity of lowest firn layer in time step
    vice = acav*dtmodel/(seconds_per_year*rhoi)    

    ! Calculate the liquid water content
    if (ImpExp == 1) then
        call Bucket_Method(ind_z_max, ind_z_surf, rhoi, Lh, Me, Mmelt, T, M, Rho, DZ, Mlwc, Refreeze, Mrunoff, Mrefreeze)
        call LWrefreeze(ind_z_max, ind_z_surf, Lh, Mrefreeze, T, M, Rho, DZ, Mlwc, Refreeze)
    endif

    ! If ice Shelf = on, vbouy has to be calculated
    if (IceShelf == 1) then
        macc = Psol + Su + Pliq
        mice = vice * rhoi
        mrun = Mrunoff
        mdiff = macc - mice - mrun
        vbouy = -1. * (mdiff/rho_ocean)   ! density ocean water = 1027 kg m-3
    else
        vbouy = 0.
    endif

    ! Update net surface elevation change
    h_surf = h_surf - vice + vfc + vacc + vsub + vsnd + vmelt + vbouy

end subroutine Update_Surface


! *******************************************************


subroutine Solve_Temp_Exp(ind_z_max,ind_z_surf,dtmodel,Ts,T,Rho,DZ,rhoi)
    !*** Solve the heat diffusion equation explicitly ***!

    ! declare arguments
    integer, intent(in) :: ind_z_max, ind_z_surf, dtmodel
    double precision, intent(in) :: Ts, rhoi
    double precision, dimension(ind_z_max), intent(in) :: Rho, DZ
    double precision, dimension(ind_z_max), intent(inout) :: T

    ! declare local variables
    integer :: ind_z
    double precision :: ci, kis, kip, f, kint, Gn, Gs

    ! j = ind_z_surf, surface layer, fixed temperature bc
    ci = 152.5 + 7.122*T(ind_z_surf)                       ! heat capacity
    f = DZ(ind_z_surf)/(DZ(ind_z_surf-1)+DZ(ind_z_surf))
    kip = Thermal_Cond(rhoi,Rho(ind_z_surf),T(ind_z_surf))        ! calc. ki of current node
    kis = Thermal_Cond(rhoi,Rho(ind_z_surf-1),T(ind_z_surf-1))    ! calc. ki of lower node
    kint = f*kis + (1.-f)*kip                       ! calc ki at interface
    Gn = -kip*(T(ind_z_surf)-Ts)*2./DZ(ind_z_surf)                ! incoming flux from upper node
    Gs = -kint*(T(ind_z_surf-1)-T(ind_z_surf))*2./(DZ(ind_z_surf-1)+DZ(ind_z_surf)) ! outgoing flux downwards
    T(ind_z_surf) = T(ind_z_surf) + dtmodel/(Rho(ind_z_surf)*ci)*(Gn-Gs)/DZ(ind_z_surf)
    ! j = 2..ind_z_surf, interior
    do ind_z = (ind_z_surf-1), 2, -1
        ci = 152.5 + 7.122*T(ind_z)
        f = DZ(ind_z)/(DZ(ind_z-1)+DZ(ind_z))
        kip = kis
        kis = Thermal_Cond(rhoi,Rho(ind_z-1),T(ind_z-1))
        kint = f*kis + (1.-f)*kip
        Gn = Gs
        Gs = -kint*(T(ind_z-1)-T(ind_z))*2./(DZ(ind_z-1)+DZ(ind_z))
        T(ind_z) = T(ind_z) + dtmodel/(Rho(ind_z)*ci)*(Gn-Gs)/DZ(ind_z)
    enddo
    ! j = 1, bottom layer, no flux bc
    ci = 152.5 + 7.122*T(1)
    Gn = Gs
    T(1) = T(1) + dtmodel/(Rho(1)*ci)*Gn/DZ(1)

end subroutine Solve_Temp_Exp


! *******************************************************


subroutine Solve_Temp_Imp(ind_z_max, ind_z_surf, dtmodel, th, Ts, T, Rho, DZ, rhoi)
    !*** Solve the heat diffusion equation (semi-)implicitly ***!
               
    ! declare arguments
    integer, intent(in) :: ind_z_max, ind_z_surf, dtmodel
    double precision, intent(in) :: th
    double precision, dimension(ind_z_max), intent(in) :: Rho, DZ
    double precision, dimension(ind_z_max), intent(inout) :: T

    ! declare local variables
    integer :: ind_z
    double precision :: ci, kip, kin, f, rhoi, Ts, as, an, ap0, Su, Sp
    double precision, dimension(ind_z_max) :: alpha, beta, D, C, A, CC

    ! see book An Introduction To Computational Fluid Dynamics by Versteeg and
    ! Malalasekera chapter 8.1 to 8.3 for implicit scheme and chapter 7.2 to 7.5
    ! for description Thomas algorithm for solving the set of equations

    ! j = 1, bottom column, zero flux bc
    ci = 152.5 + 7.122*T(1)
    f = DZ(1)/(DZ(1)+DZ(2))
    kip = Thermal_Cond(rhoi,Rho(1),T(1))
    kin = Thermal_Cond(rhoi,Rho(2),T(2))
    an = ((1.-f)*kip+f*kin)*2./(DZ(1)+DZ(2))
    ap0 = Rho(1)*ci*DZ(1)/dtmodel
    alpha(1) = an*th
    D(1) = th*an + ap0
    C(1) = an*(1.-th)*T(2) + (ap0-(1.-th)*an)*T(1)

    ! j = 2 to ind_z_surf-1, interior
    do ind_z = 2, (ind_z_surf-1)
        ci = 152.5+7.122*T(ind_z)
        f = DZ(ind_z)/(DZ(ind_z)+DZ(ind_z+1))
        kip = kin
        kin = Thermal_Cond(rhoi,Rho(ind_z+1),T(ind_z+1))
        as = an
        an = ((1.-f)*kip+f*kin)*2./(DZ(ind_z)+DZ(ind_z+1))
        ap0 = Rho(ind_z)*ci*DZ(ind_z)/dtmodel
        beta(ind_z) = as*th
        alpha(ind_z) = an*th
        D(ind_z) = th*(as+an) + ap0
        C(ind_z) = as*(1.-th)*T(ind_z-1) + an*(1.-th)*T(ind_z+1) + (ap0-(1.-th)*(as+an))*T(ind_z)
    enddo
    
    ! j = ind_z_surf, surface layer, fixed temperature bc
    ci = 152.5 + 7.122*T(ind_z_surf)
    kip = kin
    as = an
    ap0 = Rho(ind_z_surf)*ci*DZ(ind_z_surf)/dtmodel
    Su = 2.*kip*Ts/DZ(ind_z_surf)
    Sp = -2.*kip/DZ(ind_z_surf)
    beta(ind_z_surf) = as*th
    D(ind_z_surf) = th*(as-Sp) + ap0
    C(ind_z_surf) = as*(1.-th)*T(ind_z_surf-1) + (ap0-(1.-th)*(as-Sp))*T(ind_z_surf) + Su

    ! forward substitution
    A(1) = alpha(1)/D(1)
    CC(1) = C(1)/D(1)
    do ind_z = 2, (ind_z_surf-1)
        A(ind_z) = alpha(ind_z)/(D(ind_z)-beta(ind_z)*A(ind_z-1))
        CC(ind_z) = (beta(ind_z)*CC(ind_z-1)+C(ind_z))/(D(ind_z)-beta(ind_z)*A(ind_z-1))
    enddo
    CC(ind_z_surf) = (beta(ind_z_surf)*CC(ind_z_surf-1)+C(ind_z_surf))/(D(ind_z_surf)-beta(ind_z_surf)*A(ind_z_surf-1))

    ! backward substitution
    T(ind_z_surf) = CC(ind_z_surf)
    do ind_z = (ind_z_surf-1), 1, -1
        T(ind_z) = A(ind_z)*T(ind_z+1) + CC(ind_z)
    enddo

end subroutine Solve_Temp_Imp


! *******************************************************


function Thermal_Cond(rhoi, rho, Temp) result(ki)
    !*** Calculate the thermal conducticity of a layer ***!

    ! declare arguments
    double precision, intent(in) :: rhoi, rho, Temp
    double precision :: ki

    ! declare local variables
    double precision :: kice, kice_ref, kcal, kf, kair, kair_ref, theta

    kice = 9.828 * exp(-0.0057*Temp)                    ! Paterson et al., 1994
    kice_ref = 9.828 * exp(-0.0057*270.15)              ! Paterson et al., 1994
    kcal = 0.024 - 1.23E-4*rho + 2.5E-6*rho**2.
    kf = 2.107 + 0.003618*(rho-rhoi)                    ! Calonne (2019)
    kair = (2.334E-3*Temp**(3./2.))/(164.54 + Temp)       ! Reid (1966)
    kair_ref = (2.334E-3*270.15**(3./2.))/(164.54 + 270.15)
    theta = 1./(1.+exp(-0.04*(rho-450.)))
    ki = (1.-theta)*kice/kice_ref*kair/kair_ref*kcal + theta*kice/kice_ref*kf   ! Calonne (2019)
    
end function Thermal_Cond


! *******************************************************


subroutine Densific(ind_z_max, ind_z_surf, dtmodel, R, Ec, Eg, g, rhoi, acav, Rho, T, domain)
    !*** Update the density of each layer ***!

    ! declare arguments
    integer, intent(in) :: ind_z_max, ind_z_surf, dtmodel
    double precision, intent(in) :: rhoi, R, Ec, Eg, g, acav
    double precision, dimension(ind_z_max), intent(in) :: T
    double precision, dimension(ind_z_max), intent(inout) :: Rho
    character*255 :: domain

    ! declare local variables
    integer :: ind_z
    double precision :: MO, part1, Krate
    
    do ind_z = 1, ind_z_surf

        ! low density 
        if (Rho(ind_z) <= 550.) then 
            
            if (do_MO_fit) then
                MO = 1.0 ! makes doing MO fits easier - set in model_settings.f90
            elseif ((trim(domain) == "FGRN11") .or. (trim(domain) == "FGRN055")) then
                MO = 0.6688 + 0.0048*log(acav)      ! fit after debugging heat eq.
            elseif (trim(domain)=="ANT27") then
                MO = 1.288 - 0.117*log(acav)  ! Veldhuijsen et al. 2023
            else
                print *, "Domain not recognized for MO fit, using ANT27 fit"
                MO = 1.288 - 0.117*log(acav)  ! Veldhuijsen et al. 2023
            endif

            if (MO < 0.25) MO = 0.25
            
            part1 = exp((-Ec/(R*T(ind_z)))+(Eg/(R*T(1))))
            Krate = 0.07*MO*acav*g*part1 
        
        ! high density
        else

            if (do_MO_fit) then
                MO = 1 ! makes doing MO fits easier - set in model_settings.f90
            elseif ((trim(domain) == "FGRN11") .or. (trim(domain) == "FGRN055")) then
                MO = 1.7465 - 0.2045*log(acav)      ! fit after debuggin heat eq. 
            elseif (trim(domain)=="ANT27") then
                MO = 6.387 * (acav**(-0.477))+0.195   ! Veldhuijsen et al. 2023
            else
                print *, "Domain not recognized for MO fit, using ANT27 fit"
                MO = 6.387 * (acav**(-0.477))+0.195   ! Veldhuijsen et al. 2023
            endif

            if (MO < 0.25) MO = 0.25
            
            part1 = exp((-Ec/(R*T(ind_z)))+(Eg/(R*T(1))))
            Krate = 0.03*MO*acav*g*part1
        
        endif
        
        Rho(ind_z) = Rho(ind_z) + dtmodel/(seconds_per_year)*Krate*(rhoi-Rho(ind_z))
        
        if (Rho(ind_z) > rhoi) Rho(ind_z) = rhoi
        
    enddo

end subroutine Densific


! *******************************************************


subroutine Calc_Integrated_Var(ind_z_max, ind_z_surf, rhoi, Rho, Mlwc, M, DZ, FirnAir, TotLwc, IceMass)
    !*** Calculate the integrated firn air content, ice mass and total liquid water content of the firn column ***!

    ! declare arguments
    integer, intent(in) :: ind_z_max, ind_z_surf
    double precision, intent(in) :: rhoi
    double precision, intent(inout) :: FirnAir, TotLwc, IceMass
    double precision, dimension(ind_z_max), intent(in) :: Rho, Mlwc, M, DZ

    ! declare local variables
    integer :: ind_z

    FirnAir = 0.
    TotLwc = 0.
    IceMass = 0.
    do ind_z = 1, ind_z_surf
        if (Rho(ind_z) <= 910.) FirnAir = FirnAir + DZ(ind_z)*(rhoi-Rho(ind_z))/(rhoi)
        TotLwc = TotLwc + Mlwc(ind_z)
        IceMass = IceMass + M(ind_z)
    enddo

end subroutine Calc_Integrated_Var


end module firn_physics
