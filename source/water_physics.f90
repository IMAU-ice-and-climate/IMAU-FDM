module water_physics
    !*** All functions and subroutines for calculating vertical water percolation

    use define_paths_and_constants
    implicit none
    private

    public Bucket_Method, LWrefreeze
    
contains


! *******************************************************


subroutine Bucket_Method(ind_z_max, ind_z_surf, rhoi, Lh, Me, Mmelt, T, M, Rho, DZ, Mlwc, Refreeze, Mrunoff, Mrefreeze)
    !*** Distribute new liquid water using the bucket method ***!

    ! declare arguments
    integer, intent(in) :: ind_z_max, ind_z_surf
    double precision, intent(in) :: Me, rhoi, Lh
    double precision, intent(inout) :: Mmelt, Mrunoff, Mrefreeze
    double precision, dimension(ind_z_max), intent(inout) :: Rho, T, DZ, M, Mlwc, Refreeze

    ! declare local variables
    integer :: ind_z
    double precision :: cp, cp0, Efreeze, Mfreeze, Mavail, Madd, toomuch

    cp0 = 152.5 + 7.122*Tmelt

    M(ind_z_surf) = M(ind_z_surf) - Me                !Substract melted snow from upper layer
    DZ(ind_z_surf) = DZ(ind_z_surf) - (Me/Rho(ind_z_surf))   !Recalculate the height of the upper layer

    do ind_z = ind_z_surf, 1, -1   ! loop over all layers

        if (Rho(ind_z) < rhoi) then
            cp = 152.5 + 7.122*T(ind_z)
            Efreeze = (Tmelt-T(ind_z))*M(ind_z)*cp !Calculate the energy available for freezing in the layer
            Mfreeze = Efreeze / Lh          !the mass that can be frozen with that energy
            
            Mavail = Calc_Avail_Storage(ind_z, ind_z_max, rhoi, M, Rho, DZ)
    
            !Check if all the water can be refrozen in this layer
            !if YES, all the water is refrozen in this layer. T/M/dz recalculated
            if (Mfreeze >= Mmelt .and. Mavail >= Mmelt) then
                Madd = Mmelt
                Mmelt = 0.
                T(ind_z) = T(ind_z) + ((Madd*Lh) / (M(ind_z)*cp0))
                M(ind_z) = M(ind_z) + Madd
                Mrefreeze = Mrefreeze + Madd
                Refreeze(ind_z) = Refreeze(ind_z) + Madd
                Rho(ind_z) = M(ind_z) / DZ(ind_z)
            ! if NO, check which of the two is the (most) limiting factor.
            ! If it is the available energy, the refreezable part is refrozen. T == Tmelt
            ! If it is the available pore space, the pore space is filled. But T .ne. Tmelt
            ! The remainder of the meltwater will percolate further
            else
                if (Mfreeze > Mavail) then
                    Madd = Mavail
                    M(ind_z) = M(ind_z) + Madd
                    Mrefreeze = Mrefreeze + Madd
                    Refreeze(ind_z) = Refreeze(ind_z) + Madd
                    Rho(ind_z) = M(ind_z) / DZ(ind_z)
                    T(ind_z) = T(ind_z) + ((Madd*Lh) / (M(ind_z)*cp0))
                else
                    Madd = Mfreeze
                    M(ind_z) = M(ind_z) + Madd
                    Mrefreeze = Mrefreeze + Madd
                    Refreeze(ind_z) = Refreeze(ind_z) + Madd
                    Rho(ind_z) = M(ind_z) / DZ(ind_z)
                    T(ind_z) = Tmelt
                endif
                Mmelt = Mmelt - Madd
                ! Some liquid water may remain in the layer as irreducible water content
                call LWcontent(ind_z, ind_z_max, rhoi, Mmelt, M, Rho, DZ, Mlwc)
            endif
        endif

        ! Check if densification did not cause too much LWC
        Mavail = Calc_Avail_Storage(ind_z, ind_z_max, rhoi, M, Rho, DZ)
        if (Mavail < Mlwc(ind_z)) then       
            toomuch = Mlwc(ind_z) - Mavail
            Mlwc(ind_z) = Mlwc(ind_z) - toomuch
            Mmelt = Mmelt + toomuch
        endif

    enddo   ! loop over all layers
    
    ! Any remaining liquid water runs off at the bottom of the column
    if (Mmelt .ne. 0.) then
        Mrunoff = Mmelt
        Mmelt = 0.
    endif
      
end subroutine Bucket_Method


! *******************************************************


subroutine LWcontent(ind_z, ind_z_max, rhoi, Mmelt, M, Rho, DZ, Mlwc)
    !*** Update the irreducible water content of layer k ***!

    ! declare arguments
    integer, intent(in) :: ind_z, ind_z_max
    double precision, intent(in) :: rhoi
    double precision, intent(inout) :: Mmelt
    double precision, dimension(ind_z_max), intent(in) :: M, Rho, DZ
    double precision, dimension(ind_z_max), intent(inout) :: Mlwc

    ! declare local variables
    double precision :: Mavail, Mporeavail, toomuch


    Mavail = Calc_Avail_Storage(ind_z, ind_z_max, rhoi, M, Rho, DZ)

    if (Mavail < Mlwc(ind_z)) then              ! Check if densification did not cause to much LWC
        toomuch = Mlwc(ind_z) - Mavail
        Mlwc(ind_z) = Mlwc(ind_z) - toomuch
        Mmelt = Mmelt + toomuch
    else
        Mporeavail = Mavail - Mlwc(ind_z)
        if (Mporeavail >= Mmelt) then
            Mlwc(ind_z) = Mlwc(ind_z) + Mmelt       ! Enough pore space available: all water stored
            Mmelt = 0.
        else
            Mlwc(ind_z) = Mlwc(ind_z) + Mporeavail  ! Not enough pore space available: 2% pore space stored
            Mmelt = Mmelt - Mporeavail
            Mporeavail = 0.
        endif
    endif

end subroutine LWcontent


! *******************************************************


subroutine LWrefreeze(ind_z_max, ind_z_surf, Lh, Mrefreeze, T, M, Rho, DZ, Mlwc, Refreeze)
    !*** Refreeze water if layers have cooled below the melting point after temperature update ***!

    ! declare arguments
    integer, intent(in) :: ind_z_max, ind_z_surf
    double precision, intent(in) :: Lh
    double precision, intent(inout) :: Mrefreeze
    double precision, dimension(ind_z_max), intent(in) :: DZ
    double precision, dimension(ind_z_max), intent(inout) :: Rho, T, M, Mlwc, Refreeze

    ! declare local variables
    integer :: ind_z
    double precision :: cp, cp0, Mfreeze

    cp0 = 152.5 + 7.122*Tmelt

    do ind_z = 1, ind_z_surf   ! loop over all layers
        if (Mlwc(ind_z)>0 .and. T(ind_z).ne.Tmelt) then
            cp = 152.5 + 7.122*T(ind_z)
            Mfreeze = ((Tmelt-T(ind_z)) * M(ind_z) * cp) / Lh  ! Available energy for refreezing (in kgs)
            if (Mfreeze >= Mlwc(ind_z)) then
                M(ind_z) = M(ind_z) + Mlwc(ind_z)                   ! Enough energy: all LWC is refrozen
                Mrefreeze = Mrefreeze + Mlwc(ind_z)
                Refreeze(ind_z) = Refreeze(ind_z) + Mlwc(ind_z)
                Rho(ind_z) = M(ind_z) / DZ(ind_z)
                T(ind_z) = T(ind_z) + ((Mlwc(ind_z)*Lh) / (M(ind_z)*cp0))
                Mlwc(ind_z) = 0.
            else 
                M(ind_z) = M(ind_z) + Mfreeze                   ! Not enough energy: all energy is used for refreezing. Rest will remain LWC
                Mrefreeze = Mrefreeze + Mfreeze
                Refreeze(ind_z) = Refreeze(ind_z) + Mfreeze
                Rho(ind_z) = M(ind_z) / DZ(ind_z)
                T(ind_z) = Tmelt
                Mlwc(ind_z) = Mlwc(ind_z) - Mfreeze
            endif
        endif
    enddo   ! loop over all layers

end subroutine LWrefreeze


! *******************************************************


function Calc_Avail_Storage(ind_z, ind_z_max, rhoi, M, Rho, DZ) result(Mavail)
    !*** Calculate how much pore space can be stored inside of the current layer ***!

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

    Mavail = MIN(MavailCol, MavailMax)      !the available pore space (in kg) in the layer
    
end function Calc_Avail_Storage


end module water_physics
