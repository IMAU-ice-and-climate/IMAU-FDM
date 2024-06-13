module grid_routines
    !*** All subroutines and functions relating to increasing and decreasing the vertical resolution ***!

    implicit none
    private

    public :: Split_Layers, Merge_Layers, Delete_Layers, Add_Layers
    
contains


! *******************************************************


subroutine Split_Layers(ind_z_max, ind_z_surf, Rho, M, T, rgrain2, Mlwc, DZ, DenRho, Refreeze, Year, ind_t, Nt_model_tot, nyears)
    !*** Splits the uppermost layer if its mass exceeds a threshold value ***!
    
    ! declare arguments
    integer, intent(in) :: ind_z_max, ind_t, Nt_model_tot, nyears
    integer, intent(inout) :: ind_z_surf
    double precision, dimension(ind_z_max), intent(inout) :: Rho, M, T, rgrain2, Mlwc, DZ, DenRho, Refreeze, Year
    
    ! Create the new uppermost layer
    T(ind_z_surf+1) = T(ind_z_surf)
    rgrain2(ind_z_surf+1) = rgrain2(ind_z_surf)
    Rho(ind_z_surf+1) = Rho(ind_z_surf)
    DZ(ind_z_surf+1) = DZ(ind_z_surf)/2.
    M(ind_z_surf+1) = M(ind_z_surf)/2.
    Mlwc(ind_z_surf+1) = Mlwc(ind_z_surf)/2.
    DenRho(ind_z_surf+1) = DenRho(ind_z_surf)/2.
    Refreeze(ind_z_surf+1) = Refreeze(ind_z_surf)/2.
    Year(ind_z_surf+1) = (DBLE(ind_t) * DBLE(nyears)) / DBLE(Nt_model_tot)
    
    ! Update the old layer
    DZ(ind_z_surf) = DZ(ind_z_surf+1)
    M(ind_z_surf) = M(ind_z_surf+1)
    Mlwc(ind_z_surf) = Mlwc(ind_z_surf+1)
    DenRho(ind_z_surf) = DenRho(ind_z_surf+1)
    Refreeze(ind_z_surf) = Refreeze(ind_z_surf+1)
    
    ind_z_surf = ind_z_surf+1
    
end subroutine Split_Layers


! *******************************************************


subroutine Merge_Layers(ind_z_max, ind_z_surf, Rho, M, T, rgrain2, Mlwc, DZ, DenRho, Refreeze, Year)
    !*** Merges the two uppermost layers into a single layer if the thickness of the uppermost layer is smaller than a user specified threshold value ***!
    
    ! declare arguments
    integer, intent(in) :: ind_z_max
    integer, intent(inout) :: ind_z_surf
    double precision, dimension(ind_z_max), intent(inout) :: Rho, M, T, rgrain2, Mlwc, DZ, DenRho, Refreeze, Year
    
    ! Add the upper first layer into the second layer
    T(ind_z_surf-1) = (T(ind_z_surf-1)*M(ind_z_surf-1) + T(ind_z_surf)*M(ind_z_surf)) / (M(ind_z_surf-1)+M(ind_z_surf))
    rgrain2(ind_z_surf-1) = (rgrain2(ind_z_surf-1)*M(ind_z_surf-1) + rgrain2(ind_z_surf)*M(ind_z_surf)) / (M(ind_z_surf-1)+M(ind_z_surf))
    Year(ind_z_surf-1) = (Year(ind_z_surf-1)*M(ind_z_surf-1) + Year(ind_z_surf)*M(ind_z_surf)) / (M(ind_z_surf-1)+M(ind_z_surf))
    M(ind_z_surf-1) = M(ind_z_surf-1)+M(ind_z_surf)
    DZ(ind_z_surf-1) = DZ(ind_z_surf-1)+DZ(ind_z_surf)
    Rho(ind_z_surf-1) = M(ind_z_surf-1)/DZ(ind_z_surf-1)
    DenRho(ind_z_surf-1) = DenRho(ind_z_surf-1)+DenRho(ind_z_surf)
    Mlwc(ind_z_surf-1) = Mlwc(ind_z_surf-1) + Mlwc(ind_z_surf)
    Refreeze(ind_z_surf-1) = Refreeze(ind_z_surf-1) + Refreeze(ind_z_surf)

    ! Remove the old uppermost layer
    T(ind_z_surf) = 0.
    rgrain2(ind_z_surf) = 0.
    M(ind_z_surf) = 0.
    DZ(ind_z_surf) = 0.
    Rho(ind_z_surf) = 0.
    DenRho(ind_z_surf) = 0.
    Mlwc(ind_z_surf) = 0.
    Refreeze(ind_z_surf) = 0.
    Year(ind_z_surf) = 0.

    ind_z_surf = ind_z_surf-1

end subroutine Merge_Layers


! *******************************************************


subroutine Delete_Layers(ind_z_max, ind_z_surf, Rho, M, T, rgrain2, Mlwc, DZ, DenRho, Refreeze, Year)
    !*** Deletes the bottom 200 layers if they contain only ice to save memory and speed up computations ***!
    
    ! declare arguments
    integer, intent(in) :: ind_z_max
    integer, intent(inout) :: ind_z_surf
    double precision, dimension(ind_z_max), intent(inout) :: Rho, M, T, rgrain2, Mlwc, DZ, DenRho, Refreeze, Year

    ! declare local variables
    integer :: ind_z
    
    ! Move the upper 200 layers downwards
    do ind_z = 1, (ind_z_surf-200)
        T(ind_z) = T(ind_z+200)
        rgrain2(ind_z) = rgrain2(ind_z+200)
        M(ind_z) = M(ind_z+200)
        DZ(ind_z) = DZ(ind_z+200)
        Rho(ind_z) = Rho(ind_z+200)
        DenRho(ind_z) = DenRho(ind_z+200)
        Mlwc(ind_z) = Mlwc(ind_z+200)
        Refreeze(ind_z) = Refreeze(ind_z+200)
        Year(ind_z) = Year(ind_z+200)
    enddo
        
    ind_z_surf = ind_z_surf - 200
    
    ! Remove the old layers
    T(ind_z_surf+1:ind_z_surf+201) = 0.
    rgrain2(ind_z_surf+1:ind_z_surf+201) = 0.
    M(ind_z_surf+1:ind_z_surf+201) = 0.
    DZ(ind_z_surf+1:ind_z_surf+201) = 0.
    Rho(ind_z_surf+1:ind_z_surf+201) = 0.
    DenRho(ind_z_surf+1:ind_z_surf+201) = 0.
    Mlwc(ind_z_surf+1:ind_z_surf+201) = 0.
    Refreeze(ind_z_surf+1:ind_z_surf+201) = 0.
    Year(ind_z_surf+1:ind_z_surf+201) = 0.
    
end subroutine Delete_Layers


! *******************************************************


subroutine Add_Layers(ind_z_max, ind_z_surf, DZ_max, rhoi, Rho, M, T, rgrain2, Mlwc, DZ, DenRho, Refreeze, Year)
    !*** adds 100 layers of pure ice to the bottom of the column if the column is very thin. ***!
    
    ! declare arguments
    integer, intent(in) :: ind_z_max
    integer, intent(inout) :: ind_z_surf
    double precision, intent(in) :: DZ_max, rhoi
    double precision, dimension(ind_z_max), intent(inout) :: Rho, M, T, rgrain2, Mlwc, DZ, DenRho, Refreeze, Year

    ! declare local variables
    integer :: ind_z
    
    ! Move the old layers up a 100 indices 
    M(101:100+ind_z_surf)    = M(1:ind_z_surf)
    T(101:100+ind_z_surf)    = T(1:ind_z_surf)
    rgrain2(101:100+ind_z_surf) = rgrain2(1:ind_z_surf)
    DZ(101:100+ind_z_surf)   = DZ(1:ind_z_surf)
    Mlwc(101:100+ind_z_surf) = Mlwc(1:ind_z_surf)
    Rho(101:100+ind_z_surf)  = Rho(1:ind_z_surf)
    DenRho(101:100+ind_z_surf) = DenRho(1:ind_z_surf)
    Refreeze(101:100+ind_z_surf) = Refreeze(1:ind_z_surf)
    Year(101:100+ind_z_surf) = Year(1:ind_z_surf)
        
    ind_z_surf = ind_z_surf + 100
        
    ! Properties of the 100 new layers:
    do ind_z = 1, 100
        DZ(ind_z)  = DZ_max
        Rho(ind_z) = rhoi
        DenRho(ind_z) = 0.
        M(ind_z)   = Rho(ind_z) * DZ(ind_z)
        Mlwc(ind_z) = 0.
        T(ind_z) = T(101)
        rgrain2(ind_z) = rgrain2(101)
        Refreeze(ind_z) = 0.
        Year(ind_z) = -999.
    enddo
        
end subroutine Add_Layers


end module grid_routines
