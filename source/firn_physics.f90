module firn_physics
    !*** All subroutines and functions for calculating changes in density and temperature
    
    use water_physics, only: Bucket_Method, LWrefreeze

    implicit none
    private

    public :: Update_Surface, Solve_Temp_Exp, Solve_Temp_Imp, Densific, Grain_Size, Calc_Integrated_Var
    
contains


! *******************************************************


subroutine Update_Surface(ind_z_max, ind_z_surf, dtmodel, rho0, rhoi, acav, Lh, h_surf, vice, vmelt, vacc, vsub, &
    vsnd, vfc, vbouy, Ts, PSol, PLiq, Su, Me, Sd, M, T, rgrain2, DZ, Rho, DenRho, Mlwc,Refreeze, &
    ImpExp, IceShelf, Msurfmelt, Mrain, Msolin, Mrunoff, Mrefreeze)
    !*** Add and remove mass to the surface layer and calculate the velocity components ***!

    ! declare arguments
    integer, intent(in) :: ind_z_max, dtmodel, ImpExp, IceShelf
    integer, intent(inout) :: ind_z_surf
    double precision, intent(in) :: rho0, rhoi, acav, Lh
    double precision, intent(inout) :: Ts, Psol, Pliq, Su, Sd, Me
    double precision, intent(inout) :: Msurfmelt, Mrain, Msolin, Mrunoff, Mrefreeze
    double precision, intent(inout) :: vmelt, vice, h_surf, vacc, vsub, vsnd, vfc, vbouy
    double precision,dimension(ind_z_max), intent(inout) :: Rho, M, T, rgrain2, DZ, Mlwc, Refreeze, DenRho

    ! declare local variables
    integer :: ind_z
    double precision :: Madded, Mremoved, Mdiff, Macc, Mice, Mrun, Mwater, oldDZ, rgrain2_new

    ! Prepare climate input for time step
    if (ImpExp == 1) then
        if (Me > 1e-06) then
            Mwater = Me + Pliq
            Msurfmelt = Msurfmelt + Me
            Mrain = Mrain + Pliq
        elseif (Ts > 273.15) then
            Mwater = Me + Pliq
            Msurfmelt = Msurfmelt + Me
            Mrain = Mrain + Pliq
        else
            Psol = Psol + Pliq
            Pliq = 0.
            Mwater = 0.
            Me = 0.
        endif
    else
        Psol = Psol + Pliq
        Pliq = 0.
        Mwater = 0.
        Me = 0.
    endif

    ! Update surface grain size
    ! New snow has a grain size squared equal to rgrain2_new, removed firn has grain size rgrain2(ind_z_surf)
    rgrain2_new = (0.1/rho0)**2.
    if (Su > 0.) then
        if (Sd > 0.) then
            ! Riming and diverging snow drift
            Madded = Psol + Su
            Mremoved = Sd
        else
            ! Riming and converging snow drift
            Madded = Psol + Su - Sd
            Mremoved = 0.
        endif
    else
        if (Sd > 0.) then
            ! Sublimation / evaporation and diverging snow drift
            Madded = Psol
            Mremoved = -Su + Sd
        else
            ! Sublimation / evaporation and converging snow drift
            Madded = Psol - Sd
            Mremoved = -Su
        endif
    endif
    rgrain2(ind_z_surf) = (rgrain2(ind_z_surf)*(M(ind_z_surf)-Mremoved) + rgrain2_new*Madded) / (M(ind_z_surf)+Madded-Mremoved)

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

    ! Vertical (downward) velocity of lowest firn layer in time step
    vice = acav*dtmodel/(3600.*24.*365.*rhoi)    

    ! Melt
    vmelt = -1 * Me/Rho(ind_z_surf)
    M(ind_z_surf) = M(ind_z_surf) - Me                      ! Substract melted snow from upper layer
    DZ(ind_z_surf) = DZ(ind_z_surf) - (Me/Rho(ind_z_surf))  ! Recalculate the height of the upper layer

    ! Calculate the liquid water content
    if (ImpExp == 1) then
        if (Mwater > 0.) then
            call Bucket_Method(ind_z_max, ind_z_surf, rhoi, Lh, Mwater, T, rgrain2, M, Rho, DZ, Mlwc, Refreeze, Mrunoff, Mrefreeze)
        endif
        call LWrefreeze(ind_z_max, ind_z_surf, rhoi, Lh, Mrefreeze, T, rgrain2, M, Rho, DZ, Mlwc, Refreeze)
    endif

    ! If ice Shelf = on, vbouy has to be calculated
    if (IceShelf == 1) then
        Macc = Psol + Su - Sd + Pliq
        Mice = vice * rhoi
        Mrun = Mrunoff
        Mdiff = Macc - Mice - Mrun
        vbouy = -1. * (Mdiff/1027.)   ! density water = 1027 kg m-3
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
    kair = (2.334E-3*Temp**(3/2))/(164.54 + Temp)       ! Reid (1966)
    kair_ref = (2.334E-3*270.15**(3/2))/(164.54 + 270.15)
    theta = 1./(1.+exp(-0.04*(rho-450.)))
    ki = (1.-theta)*kice/kice_ref*kair/kair_ref*kcal + theta*kice/kice_ref*kf   ! Calonne (2019)
    
end function Thermal_Cond


! *******************************************************


subroutine Densific(ind_z_max, ind_z_surf, dtmodel, R, Ec, Ec2, g, rhoi, Rho, M, T, rgrain2)
    !*** Update the density of each layer ***!

    ! declare arguments
    integer, intent(in) :: ind_z_max, ind_z_surf, dtmodel
    double precision, intent(in) :: rhoi, R, Ec, Ec2, g
    double precision, dimension(ind_z_max), intent(in) :: M, T, rgrain2
    double precision, dimension(ind_z_max), intent(inout) :: Rho

    ! declare local variables
    integer :: ind_z
    double precision :: a, Rho_t, A1, A2, P, Ps, Ptot, dRho_part1, dRho_part2, theta, w1, w2, dRho, rgrain, kc

    ! fitting constants
    a = 0.06
    Rho_t = 550.
    !A1 = 9.268E-9
    !A2 = 8.903E-14
    A1 = 17.912E-9
    A2 = 8.869E-14

    !!! Update uppermost layer

    if (Rho(ind_z_surf) < rhoi) then

        rgrain = rgrain2(ind_z_surf)**0.5       ! grain size in m

        P = g*M(ind_z_surf)/2.                  ! load per unit area
        Ps = rhoi/Rho(ind_z_surf)*P             ! grain-load stress: correction for the fact that only grains, not voids, support the load
        Ptot = (rhoi/Rho(ind_z_surf)-1.)*Ps     ! effective stress

        ! Calculate the components of the densification rate
        dRho_part1 = A1*exp(-Ec/(R*T(ind_z_surf)))*Ptot/(rgrain**2.)
        dRho_part2 = A2*exp(-Ec2/(R*T(ind_z_surf)))*(Ptot**1.8)/(rgrain**1.4)
        
        ! Calculate weighing factors
        theta = 1./(1.+exp(-a*(Rho(ind_z_surf) - Rho_t)))
        w1 = 1. - theta
        w2 = theta

        ! Calculate the total densification
        dRho = (w1*dRho_part1 + w2*dRho_part2) * Rho(ind_z_surf) * dtmodel

        ! Update density
        Rho(ind_z_surf) = Rho(ind_z_surf) + dRho
        if (Rho(ind_z_surf) > rhoi) Rho(ind_z_surf) = rhoi

    end if

    !!! Update remaining layers

    do ind_z = (ind_z_surf-1), 1, -1
        
        if (Rho(ind_z) < rhoi) then

            rgrain = rgrain2(ind_z)**0.5            ! grain size in m
            
            ! Calculate the driving stress
            P = P + g*(M(ind_z+1)+M(ind_z))/2.      ! load per unit area
            Ps = rhoi/Rho(ind_z)*P                  ! grain-load stress: correction for the fact that only grains, not voids, support the load
            Ptot = (rhoi/Rho(ind_z)-1.)*Ps          ! effective stress

            ! Calculate the components of the densification rate
            dRho_part1 = A1*exp(-Ec/(R*T(ind_z)))*Ptot/(rgrain**2.)
            dRho_part2 = A2*exp(-Ec2/(R*T(ind_z)))*(Ptot**1.8)/(rgrain**1.4)

            ! Calculate weighing factors
            theta = 1./(1.+exp(-a*(Rho(ind_z) - Rho_t)))
            w1 = 1. - theta
            w2 = theta

            ! Calculate the total densification
            dRho = (w1*dRho_part1 + w2*dRho_part2) * Rho(ind_z) * dtmodel

            ! Update density
            Rho(ind_z) = Rho(ind_z) + dRho
            if (Rho(ind_z) > rhoi) Rho(ind_z) = rhoi

        end if

    enddo

end subroutine Densific


! *******************************************************


subroutine Grain_Size(ind_z_max, ind_z_surf, dtmodel, kg, Eg, R, T, rgrain2)
    !*** Update the grain size (squared) of each layer ***!

    ! declare arguments
    integer, intent(in) :: ind_z_max, ind_z_surf, dtmodel
    double precision, intent(in) :: kg, Eg, R
    double precision, dimension(ind_z_max), intent(in) :: T
    double precision, dimension(ind_z_max), intent(inout) :: rgrain2

    rgrain2(1:ind_z_surf) = rgrain2(1:ind_z_surf) + kg*exp(-Eg/(R*T(1:ind_z_surf)))*dtmodel

end subroutine Grain_Size


! *******************************************************


subroutine Calc_Integrated_Var(ind_z_max, ind_z_surf, rhoi, Rho, Mlwc, M, DZ, FirnAir, TotLwc, IceMass)
    !*** Calculate the total firn air content, ice mass and liquid water content of the firn column ***!

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
