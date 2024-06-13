module water_physics
    !*** All functions and subroutines for calculating vertical water percolation

    implicit none
    private

    public Bucket_Method, LWrefreeze
    
contains


! *******************************************************


subroutine Bucket_Method(ind_z_max, ind_z_surf, rhoi, Lh, Mwater, T, rgrain2, M, Rho, DZ, Mlwc, Refreeze, Mrunoff, Mrefreeze)
    !*** Distribute new liquid water using the bucket method ***!

    ! declare arguments
    integer, intent(in) :: ind_z_max, ind_z_surf
    double precision, intent(in) :: rhoi, Lh
    double precision, intent(inout) :: Mwater, Mrunoff, Mrefreeze
    double precision, dimension(ind_z_max), intent(inout) :: Rho, T, rgrain2, DZ, M, Mlwc, Refreeze

    ! declare local variables
    integer :: ind_z
    double precision :: cp, cp0, Efreeze, Mfreeze, Mavail, MavailMax, Madd, toomuch, rgrain2_refreeze

    cp0 = 152.5 + 7.122*273.15
    rgrain2_refreeze = (1.E-3)**2.

    ind_z = ind_z_surf
    do while ((Mwater > 0.) .and. (ind_z > 0))

        if (Rho(ind_z) < rhoi) then
            cp = 152.5 + 7.122*T(ind_z)
            Efreeze = (273.15-T(ind_z))*M(ind_z)*cp     ! Calculate the energy available for freezing in the layer
            Mfreeze = Efreeze / Lh                      ! The mass that can be frozen with that energy
            
            MavailMax = rhoi*DZ(ind_z) - M(ind_z)
    
            ! Check if all the water can be refrozen in this layer
            ! if YES, all the water is refrozen in this layer. T/M/dz recalculated
            if ((Mfreeze >= Mwater) .and. (MavailMax >= Mwater)) then
                Madd = Mwater
                Mwater = 0.
                T(ind_z) = T(ind_z) + ((Madd*Lh) / (M(ind_z)*cp0))
                rgrain2(ind_z) = (rgrain2(ind_z)*M(ind_z) + rgrain2_refreeze*Madd ) / ( M(ind_z) + Madd )
                M(ind_z) = M(ind_z) + Madd
                Mrefreeze = Mrefreeze + Madd
                Refreeze(ind_z) = Refreeze(ind_z) + Madd
                Rho(ind_z) = M(ind_z) / DZ(ind_z)
            ! if NO, check which of the two is the (most) limiting factor.
            ! If it is the available energy, the refreezable part is refrozen. T == 273.15
            ! If it is the available pore space, the pore space is filled. But T .ne. 273.15
            ! The remainder of the meltwater will percolate further
            else
                if (Mfreeze > MavailMax) then
                    Madd = MavailMax
                    rgrain2(ind_z) = (rgrain2(ind_z)*M(ind_z) + rgrain2_refreeze*Madd ) / ( M(ind_z) + Madd )
                    M(ind_z) = M(ind_z) + Madd
                    Mrefreeze = Mrefreeze + Madd
                    Refreeze(ind_z) = Refreeze(ind_z) + Madd
                    Rho(ind_z) = M(ind_z) / DZ(ind_z)
                    T(ind_z) = T(ind_z) + ((Madd*Lh) / (M(ind_z)*cp0))
                else
                    Madd = Mfreeze
                    rgrain2(ind_z) = (rgrain2(ind_z)*M(ind_z) + rgrain2_refreeze*Madd ) / ( M(ind_z) + Madd )
                    M(ind_z) = M(ind_z) + Madd
                    Mrefreeze = Mrefreeze + Madd
                    Refreeze(ind_z) = Refreeze(ind_z) + Madd
                    Rho(ind_z) = M(ind_z) / DZ(ind_z)
                    T(ind_z) = 273.15
                endif
                Mwater = Mwater - Madd
                ! Some liquid water may remain in the layer as irreducible water content
                call LWcontent(ind_z, ind_z_max, rhoi, Mwater, M, Rho, DZ, Mlwc)
            endif
        endif

        ! Check if refreezing did not cause too high densities
        if (Rho(ind_z) > rhoi) then
            Rho(ind_z) = rhoi
        endif

        ! Check if refreezing did not cause too much LWC
        Mavail = Calc_Avail_Storage(ind_z, ind_z_max, rhoi, M, Rho, DZ)
        if (Mavail < Mlwc(ind_z)) then       
            toomuch = Mlwc(ind_z) - Mavail
            Mlwc(ind_z) = Mlwc(ind_z) - toomuch
            Mwater = Mwater + toomuch
        endif

        ind_z = ind_z - 1

    enddo   ! loop over all layers
    
    ! Any remaining liquid water runs off at the bottom of the column
    if (Mwater .ne. 0.) then
        Mrunoff = Mwater
        Mwater = 0.
    endif
      
end subroutine Bucket_Method


! *******************************************************


subroutine LWcontent(ind_z, ind_z_max, rhoi, Mwater, M, Rho, DZ, Mlwc)
    !*** Update the irreducible water content of layer k ***!

    ! declare arguments
    integer, intent(in) :: ind_z, ind_z_max
    double precision, intent(in) :: rhoi
    double precision, intent(inout) :: Mwater
    double precision, dimension(ind_z_max), intent(in) :: M, Rho, DZ
    double precision, dimension(ind_z_max), intent(inout) :: Mlwc

    ! declare local variables
    double precision :: Mavail, Mporeavail, toomuch


    Mavail = Calc_Avail_Storage(ind_z, ind_z_max, rhoi, M, Rho, DZ)

    if (Mavail < Mlwc(ind_z)) then              ! Check if refreezing did not cause to much LWC
        toomuch = Mlwc(ind_z) - Mavail
        Mlwc(ind_z) = Mlwc(ind_z) - toomuch
        Mwater = Mwater + toomuch
    else
        Mporeavail = Mavail - Mlwc(ind_z)
        if (Mporeavail >= Mwater) then
            Mlwc(ind_z) = Mlwc(ind_z) + Mwater       ! Enough pore space available: all water stored
            Mwater = 0.
        else
            Mlwc(ind_z) = Mlwc(ind_z) + Mporeavail  ! Not enough pore space available: 2% pore space stored
            Mwater = Mwater - Mporeavail
            Mporeavail = 0.
        endif
    endif

end subroutine LWcontent


! *******************************************************


subroutine LWrefreeze(ind_z_max, ind_z_surf, rhoi, Lh, Mrefreeze, T, rgrain2, M, Rho, DZ, Mlwc, Refreeze)
    !*** Refreeze water if layers have cooled below the melting point after temperature update ***!

    ! declare arguments
    integer, intent(in) :: ind_z_max, ind_z_surf
    double precision, intent(in) :: rhoi, Lh
    double precision, intent(inout) :: Mrefreeze
    double precision, dimension(ind_z_max), intent(in) :: DZ
    double precision, dimension(ind_z_max), intent(inout) :: Rho, T, rgrain2, M, Mlwc, Refreeze

    ! declare local variables
    integer :: ind_z
    !double precision :: cp, cp0, Mfreeze, rgrain_refreeze
    double precision :: cp, cp0, Mfreeze, rgrain2_refreeze

    cp0 = 152.5 + 7.122*273.15
    rgrain2_refreeze = (1.E-3)**2.

    do ind_z = 1, ind_z_surf   ! loop over all layers
        if ((Mlwc(ind_z)>0.) .and. (T(ind_z).ne.273.15)) then
            cp = 152.5 + 7.122*T(ind_z)
            Mfreeze = ((273.15-T(ind_z)) * M(ind_z) * cp) / Lh  ! Available energy for refreezing (in kgs)
            if (Mfreeze >= Mlwc(ind_z)) then
                rgrain2(ind_z) = (rgrain2(ind_z)*M(ind_z) + rgrain2_refreeze*Mlwc(ind_z) ) / ( M(ind_z) + Mlwc(ind_z) )
                M(ind_z) = M(ind_z) + Mlwc(ind_z)                   ! Enough energy: all LWC is refrozen
                Mrefreeze = Mrefreeze + Mlwc(ind_z)
                Refreeze(ind_z) = Refreeze(ind_z) + Mlwc(ind_z)
                Rho(ind_z) = M(ind_z) / DZ(ind_z)
                T(ind_z) = T(ind_z) + ((Mlwc(ind_z)*Lh) / (M(ind_z)*cp0))
                Mlwc(ind_z) = 0.
            else 
                rgrain2(ind_z) = (rgrain2(ind_z)*M(ind_z) + rgrain2_refreeze*Mrefreeze ) / ( M(ind_z) + Mrefreeze )
                M(ind_z) = M(ind_z) + Mfreeze                   ! Not enough energy: all energy is used for refreezing. Rest will remain LWC
                Mrefreeze = Mrefreeze + Mfreeze
                Refreeze(ind_z) = Refreeze(ind_z) + Mfreeze
                Rho(ind_z) = M(ind_z) / DZ(ind_z)
                T(ind_z) = 273.15
                Mlwc(ind_z) = Mlwc(ind_z) - Mfreeze
            endif
            if (Rho(ind_z) > rhoi) Rho(ind_z) = rhoi
        endif
    enddo   ! loop over all layers

end subroutine LWrefreeze


! *******************************************************


function Calc_Avail_Storage(ind_z, ind_z_max, rhoi, M, Rho, DZ) result(Mavail)
    !*** Calculate how much liquid water can be stored inside the current layer ***!

    ! declare arguments
    integer, intent(in) :: ind_z, ind_z_max
    double precision, intent(in) :: rhoi
    double precision, dimension(ind_z_max), intent(in) :: M, Rho, DZ

    ! declare local variables
    double precision :: poro, maxpore, MavailCol, MavailMax, Mavail

    ! Maximum available capacity for liquid water according to Coleou, 1998
    poro = (rhoi-Rho(ind_z))/rhoi
    maxpore = 0.017 + 0.057 * (poro/(1.-poro))
    MavailCol = maxpore * M(ind_z)

    ! Maximum available capacity for liquid water for high density firn
    ! Coleou, 1998 parameterization still has 1.7% of water for the density of ice...
    MavailMax = rhoi*DZ(ind_z) - M(ind_z)

    Mavail = MIN(MavailCol, MavailMax)      ! The available pore space (in kg) in the layer
    
end function Calc_Avail_Storage


end module water_physics
