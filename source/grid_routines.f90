module grid_routines
    !*** All subroutines and functions relating to increasing and decreasing the vertical resolution ***!

    use model_settings

    implicit none
    private

    public :: Split_Surface_Layer_By_Threshold, Merge_Surface_Layer_By_Threshold, Delete_Layers, Add_Layers, &
        Split_Surface_Layer_By_Thickness_Range, Merge_Surface_Layer_By_Thickness_Range, Split_Second_Layer_By_Threshold, &
        Merge_Second_Layer_By_Threshold, Remove_Surface_Layer
    
contains


! *******************************************************


subroutine Split_Surface_Layer_By_Threshold(ind_z_max, ind_z_surf, Rho, M, T, Mlwc, DZ, DenRho, Refreeze, Year, ind_t, Nt_model_tot, nyears)
    !*** Splits the uppermost layer if its mass exceeds a threshold value ***!
    
    ! declare arguments
    integer, intent(in) :: ind_z_max, ind_t, Nt_model_tot, nyears
    integer, intent(inout) :: ind_z_surf
    double precision, dimension(ind_z_max), intent(inout) :: Rho, M, T, Mlwc, DZ, DenRho, Refreeze, Year
    
    ! Create the new uppermost layer
    T(ind_z_surf+1) = T(ind_z_surf)
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
    
end subroutine Split_Surface_Layer_By_Threshold


! *******************************************************


subroutine Merge_Surface_Layer_By_Threshold(ind_z_max, ind_z_surf, Rho, M, T, Mlwc, DZ, DenRho, Refreeze, Year)
    !*** Merges the two uppermost layers into a single layer if the thickness of the uppermost layer is smaller than a user specified threshold value ***!
    
    ! declare arguments
    integer, intent(in) :: ind_z_max
    integer, intent(inout) :: ind_z_surf
    double precision, dimension(ind_z_max), intent(inout) :: Rho, M, T, Mlwc, DZ, DenRho, Refreeze, Year
    
    ! Add the upper first layer into the second layer
    T(ind_z_surf-1) = (T(ind_z_surf-1)*M(ind_z_surf-1) + T(ind_z_surf)*M(ind_z_surf)) / (M(ind_z_surf-1)+M(ind_z_surf))
    Year(ind_z_surf-1) = (Year(ind_z_surf-1)*M(ind_z_surf-1) + Year(ind_z_surf)*M(ind_z_surf)) / (M(ind_z_surf-1)+M(ind_z_surf))
    M(ind_z_surf-1) = M(ind_z_surf-1)+M(ind_z_surf)
    DZ(ind_z_surf-1) = DZ(ind_z_surf-1)+DZ(ind_z_surf)
    Rho(ind_z_surf-1) = M(ind_z_surf-1)/DZ(ind_z_surf-1)
    DenRho(ind_z_surf-1) = DenRho(ind_z_surf-1)+DenRho(ind_z_surf)
    Mlwc(ind_z_surf-1) = Mlwc(ind_z_surf-1) + Mlwc(ind_z_surf)
    Refreeze(ind_z_surf-1) = Refreeze(ind_z_surf-1) + Refreeze(ind_z_surf)

    ! Remove the old uppermost layer
    T(ind_z_surf) = 0.
    M(ind_z_surf) = 0.
    DZ(ind_z_surf) = 0.
    Rho(ind_z_surf) = 0.
    DenRho(ind_z_surf) = 0.
    Mlwc(ind_z_surf) = 0.
    Refreeze(ind_z_surf) = 0.
    Year(ind_z_surf) = 0.

    ind_z_surf = ind_z_surf-1

end subroutine Merge_Surface_Layer_By_Threshold


! *******************************************************


subroutine Delete_Layers(ind_z_max, ind_z_surf, Rho, M, T, Mlwc, DZ, DenRho, Refreeze, Year, rgrain2)
    !*** Deletes the bottom 200 layers if they contain only ice to save memory and speed up computations ***!
    
    ! declare arguments
    integer, intent(in) :: ind_z_max
    integer, intent(inout) :: ind_z_surf
    double precision, dimension(ind_z_max), intent(inout) :: Rho, M, T, Mlwc, DZ, DenRho, Refreeze, Year, rgrain2

    ! declare local variables
    integer :: ind_z
    
    ! Move the upper 200 layers downwards
    do ind_z = 1, (ind_z_surf-200)
        T(ind_z) = T(ind_z+200)
        M(ind_z) = M(ind_z+200)
        DZ(ind_z) = DZ(ind_z+200)
        Rho(ind_z) = Rho(ind_z+200)
        DenRho(ind_z) = DenRho(ind_z+200)
        Mlwc(ind_z) = Mlwc(ind_z+200)
        Refreeze(ind_z) = Refreeze(ind_z+200)
        Year(ind_z) = Year(ind_z+200)
        if (grainsize_veldhuijsen) then
            rgrain2(ind_z) = rgrain2(ind_z+200)
        endif
    enddo
        
    ind_z_surf = ind_z_surf - 200
    
    ! Remove the old layers
    T(ind_z_surf+1:ind_z_surf+201) = 0.
    M(ind_z_surf+1:ind_z_surf+201) = 0.
    DZ(ind_z_surf+1:ind_z_surf+201) = 0.
    Rho(ind_z_surf+1:ind_z_surf+201) = 0.
    DenRho(ind_z_surf+1:ind_z_surf+201) = 0.
    Mlwc(ind_z_surf+1:ind_z_surf+201) = 0.
    Refreeze(ind_z_surf+1:ind_z_surf+201) = 0.
    Year(ind_z_surf+1:ind_z_surf+201) = 0.
    if (grainsize_veldhuijsen) then
        rgrain2(ind_z_surf+1:ind_z_surf+201) = 0.
    endif
    
end subroutine Delete_Layers


! *******************************************************


subroutine Add_Layers(ind_z_max, ind_z_surf, dzmax, rhoi, Rho, M, T, Mlwc, DZ, DenRho, Refreeze, Year, rgrain2)
    !*** adds 100 layers of pure ice to the bottom of the column if the column is very thin. ***!
    
    ! declare arguments
    integer, intent(in) :: ind_z_max
    integer, intent(inout) :: ind_z_surf
    double precision, intent(in) :: dzmax, rhoi
    double precision, dimension(ind_z_max), intent(inout) :: Rho, M, T, Mlwc, DZ, DenRho, Refreeze, Year, rgrain2

    ! declare local variables
    integer :: ind_z
    
    ! Move the old layers up a 100 indices 
    M(101:100+ind_z_surf)    = M(1:ind_z_surf)
    T(101:100+ind_z_surf)    = T(1:ind_z_surf)
    DZ(101:100+ind_z_surf)   = DZ(1:ind_z_surf)
    Mlwc(101:100+ind_z_surf) = Mlwc(1:ind_z_surf)
    Rho(101:100+ind_z_surf)  = Rho(1:ind_z_surf)
    DenRho(101:100+ind_z_surf) = DenRho(1:ind_z_surf)
    Refreeze(101:100+ind_z_surf) = Refreeze(1:ind_z_surf)
    Year(101:100+ind_z_surf) = Year(1:ind_z_surf)
    if (grainsize_veldhuijsen) then
        rgrain2(101:100+ind_z_surf) = rgrain2(1:ind_z_surf)
    endif

    ind_z_surf = ind_z_surf + 100
        
    ! Properties of the 100 new layers:
    do ind_z = 1, 100
        DZ(ind_z)  = dzmax
        Rho(ind_z) = rhoi
        DenRho(ind_z) = 0.
        M(ind_z)   = Rho(ind_z) * DZ(ind_z)
        Mlwc(ind_z) = 0.
        T(ind_z) = T(101)
        Refreeze(ind_z) = 0.
        Year(ind_z) = -999.
        if (grainsize_veldhuijsen) then
            rgrain2(ind_z) = (0.01_8)**2
            Year(ind_z) = 999.*3600*24*365
        endif
    enddo
        
end subroutine Add_Layers


! *******************************************************

  subroutine Split_Surface_Layer_By_Thickness_Range(ind_z_max,dzmax_upper_layer,ind_z_surf,dzmax,Rho,M,T,Mlwc,DZ,DenRho,Refreeze,Year,rgrain2)
    !*** Subroutine that splits the first layer if it exceeds a range. This is build in for grainsize_veldhuijsen. The max resolution of the first layer is smaller than the rest. 

    !declare arguments
    integer, intent(in) :: ind_z_max 
    integer, intent(inout) :: ind_z_surf
    double precision, intent(in) :: dzmax, dzmax_upper_layer
    double precision, dimension(ind_z_max), intent(inout) :: Rho, M, T, Mlwc, DZ, DenRho, Refreeze, Year, rgrain2

    ! declare local variables
    integer :: DZ_to_second_layer, M_to_second_layer
  
    DZ_to_second_layer = DZ(ind_z_surf) - dzmax_upper_layer
    M_to_second_layer = DZ_to_second_layer * Rho(ind_z_surf)

    T(ind_z_surf-1) = (T(ind_z_surf-1)*M(ind_z_surf-1) + T(ind_z_surf)*M_to_second_layer)/(M(ind_z_surf-1)+M_to_second_layer)
    Rho(ind_z_surf-1) = (Rho(ind_z_surf-1)*M(ind_z_surf-1) + Rho(ind_z_surf)*M_to_second_layer)/(M(ind_z_surf-1)+M_to_second_layer)
    Year(ind_z_surf-1) = (Year(ind_z_surf-1)*M(ind_z_surf-1) + Year(ind_z_surf)*M_to_second_layer)/(M(ind_z_surf-1)+M_to_second_layer)
    rgrain2(ind_z_surf-1) = (rgrain2(ind_z_surf-1)*M(ind_z_surf-1) + rgrain2(ind_z_surf)*M_to_second_layer)/(M(ind_z_surf-1)+M_to_second_layer)
    Mlwc(ind_z_surf-1) = Mlwc(ind_z_surf-1) + (Mlwc(ind_z_surf) * DZ_to_second_layer / DZ(ind_z_surf))
    DenRho(ind_z_surf-1) = DenRho(ind_z_surf-1) + (DenRho(ind_z_surf) * DZ_to_second_layer/ DZ(ind_z_surf))
    Refreeze(ind_z_surf-1) =  Refreeze(ind_z_surf-1) + (Refreeze(ind_z_surf) * DZ_to_second_layer / DZ(ind_z_surf))
    DZ(ind_z_surf-1) = DZ(ind_z_surf-1)+ DZ_to_second_layer
    M(ind_z_surf-1) = M(ind_z_surf-1) + M_to_second_layer

    Mlwc(ind_z_surf) = Mlwc(ind_z_surf) * dzmax_upper_layer/DZ(ind_z_surf)
    DenRho(ind_z_surf) = DenRho(ind_z_surf) * dzmax_upper_layer/DZ(ind_z_surf)
    Refreeze(ind_z_surf) = Refreeze(ind_z_surf) * dzmax_upper_layer/DZ(ind_z_surf)
    DZ(ind_z_surf) = dzmax_upper_layer 
    M(ind_z_surf) = Rho(ind_z_surf) * DZ(ind_z_surf) 

    end subroutine Split_Surface_Layer_By_Thickness_Range


! *******************************************************
  subroutine Merge_Surface_Layer_By_Thickness_Range(ind_z_max,dzmax_upper_layer,ind_z_surf,Rho,M,T,Mlwc,DZ,DenRho,Refreeze,Year,rgrain2)
    !*** Merges the top two layers if the surface layer is below a threshold. This is build in for grainsize_veldhuijsen. The min resolution of the first layer is different from the rest. 

    !declare arguments
    integer, intent(in) :: ind_z_max
    integer, intent(inout) :: ind_z_surf
    double precision, intent(in) :: dzmax_upper_layer
    double precision, dimension(ind_z_max), intent(inout) :: Rho, M, T, Mlwc, DZ, DenRho, Refreeze, Year, rgrain2

    ! declare local variables
    integer :: dz_merge, m_merge

    dz_merge = dzmax_upper_layer - DZ(ind_z_surf) 
    m_merge = M(ind_z_surf-1)*dz_merge/DZ(ind_z_surf-1)

    T(ind_z_surf) = (T(ind_z_surf-1)*m_merge + T(ind_z_surf)*M(ind_z_surf)) / (m_merge+M(ind_z_surf))
    Year(ind_z_surf) =  (Year(ind_z_surf-1)*m_merge + Year(ind_z_surf)*M(ind_z_surf))/(m_merge+M(ind_z_surf))
    rgrain2(ind_z_surf) = (rgrain2(ind_z_surf-1)*m_merge + rgrain2(ind_z_surf)*M(ind_z_surf)) / (m_merge+M(ind_z_surf)) 
    Rho(Ind_z_surf) = (Rho(ind_z_surf-1)*m_merge + rgrain2(ind_z_surf)*M(ind_z_surf)) / (m_merge+M(ind_z_surf))
    DenRho(ind_z_surf) = DenRho(ind_z_surf-1)*dz_merge/DZ(ind_z_surf-1) + DenRho(ind_z_surf)
    Mlwc(ind_z_surf) = Mlwc(ind_z_surf-1)*dz_merge/DZ(ind_z_surf-1) + Mlwc(ind_z_surf)
    Refreeze(ind_z_surf) = Refreeze(ind_z_surf-1)*dz_merge/DZ(ind_z_surf-1) + Refreeze(ind_z_surf)
    M(ind_z_surf) = M(ind_z_surf)+m_merge
    DZ(ind_z_surf) = DZ(ind_z_surf)+dz_merge
   
    DenRho(ind_z_surf-1) = DenRho(ind_z_surf-1)*(M(ind_z_surf-1)-m_merge)/M(ind_z_surf-1)
    Mlwc(ind_z_surf-1) = Mlwc(ind_z_surf-1)*(M(ind_z_surf-1)-m_merge)/M(ind_z_surf-1)
    Refreeze(ind_z_surf-1) = Refreeze(ind_z_surf-1)*(M(ind_z_surf-1)-m_merge)/M(ind_z_surf-1)
    M(ind_z_surf-1)= M(ind_z_surf-1)-m_merge
    DZ(ind_z_surf-1) = DZ(ind_z_surf-1)-dz_merge
    
    end subroutine Merge_Surface_Layer_By_Thickness_Range

! *******************************************************
  subroutine Split_Second_Layer_By_Threshold(ind_z_max, ind_z_surf, Rho, M, T, Mlwc, DZ, DenRho, Refreeze, Year, rgrain2)
    !*** the uppermost layer stays 1 cm, so the second layer splits it its mass exceeds a threshold value

    ! declare arguments
    integer, intent(in) :: ind_z_max
    integer, intent(inout) :: ind_z_surf
    double precision, dimension(ind_z_max), intent(inout) :: Rho, M, T, Mlwc, DZ, DenRho, Refreeze, Year, rgrain2

    ! Create the new uppermost layer
    T(ind_z_surf+1) = T(ind_z_surf)
    Rho(ind_z_surf+1) = Rho(ind_z_surf)
    DZ(ind_z_surf+1) = DZ(ind_z_surf)
    M(ind_z_surf+1) = M(ind_z_surf)
    Mlwc(ind_z_surf+1) = Mlwc(ind_z_surf)
    DenRho(ind_z_surf+1) = DenRho(ind_z_surf)
    Refreeze(ind_z_surf+1) = Refreeze(ind_z_surf)
    rgrain2(ind_z_surf+1) = rgrain2(ind_z_surf)
    Year(ind_z_surf+1) = Year(ind_z_surf)

    ! Create the new second layer
    T(ind_z_surf) = T(ind_z_surf-1)
    Rho(ind_z_surf) = Rho(ind_z_surf-1)
    DZ(ind_z_surf) = DZ(ind_z_surf-1)/2.
    M(ind_z_surf) = M(ind_z_surf-1)/2.
    Mlwc(ind_z_surf) = Mlwc(ind_z_surf-1)/2.
    DenRho(ind_z_surf) = DenRho(ind_z_surf-1)/2.
    Refreeze(ind_z_surf) = Refreeze(ind_z_surf-1)/2.
    Year(ind_z_surf) = Year(ind_z_surf-1)
    rgrain2(ind_z_surf) = rgrain2(ind_z_surf-1)

    ! Update the third layer
    DZ(ind_z_surf-1) = DZ(ind_z_surf)
    M(ind_z_surf-1) = M(ind_z_surf)
    Mlwc(ind_z_surf-1) = Mlwc(ind_z_surf)
    DenRho(ind_z_surf-1) = DenRho(ind_z_surf)
    Refreeze(ind_z_surf-1) = Refreeze(ind_z_surf)

    ind_z_surf = ind_z_surf + 1

    end subroutine Split_Second_Layer_By_Threshold

! *******************************************************
  subroutine Merge_Second_Layer_By_Threshold(ind_z_max, ind_z_surf, Rho, M, T, Mlwc, DZ, DenRho, Refreeze, Year, rgrain2)
    !*** the uppermost layer stays 1 cm, so the second layer maerges if its mass is lower than a threshold value

    ! declare arguments
    integer, intent(in) :: ind_z_max
    integer, intent(inout) :: ind_z_surf
    double precision, dimension(ind_z_max), intent(inout) :: Rho, M, T, Mlwc, DZ, DenRho, Refreeze, Year, rgrain2

    ! Add the second layer into the third layer
    T(ind_z_surf-2) = (T(ind_z_surf-2)*M(ind_z_surf-2) + T(ind_z_surf-1)*M(ind_z_surf-1)) / (M(ind_z_surf-2)+M(ind_z_surf-1))
    Year(ind_z_surf-2) = (Year(ind_z_surf-2)*M(ind_z_surf-2) + Year(ind_z_surf-1)*M(ind_z_surf-1)) / (M(ind_z_surf-2)+M(ind_z_surf-1)) 
    rgrain2(ind_z_surf-2) =  (rgrain2(ind_z_surf-2)*M(ind_z_surf-2) + rgrain2(ind_z_surf-1)*M(ind_z_surf-1)) / (M(ind_z_surf-2)+M(ind_z_surf-1))
    M(ind_z_surf-2) = M(ind_z_surf-2)+M(ind_z_surf-1)
    DZ(ind_z_surf-2) = DZ(ind_z_surf-2)+DZ(ind_z_surf-1)
    Rho(ind_z_surf-2) = M(ind_z_surf-2)/DZ(ind_z_surf-2)
    DenRho(ind_z_surf-2) = DenRho(ind_z_surf-2)+DenRho(ind_z_surf-1)
    Mlwc(ind_z_surf-2) = Mlwc(ind_z_surf-2) + Mlwc(ind_z_surf-1)
    Refreeze(ind_z_surf-2) = Refreeze(ind_z_surf-2) + Refreeze(ind_z_surf-1)

    ! Move the surface layer one down
    T(ind_z_surf-1) = T(ind_z_surf)
    Rho(ind_z_surf-1) = Rho(ind_z_surf)
    DZ(ind_z_surf-1) = DZ(ind_z_surf)
    M(ind_z_surf-1) = M(ind_z_surf)
    Mlwc(ind_z_surf-1) = Mlwc(ind_z_surf)
    DenRho(ind_z_surf-1) = DenRho(ind_z_surf)
    Refreeze(ind_z_surf-1) = Refreeze(ind_z_surf)
    Year(ind_z_surf-1) = Year(ind_z_surf)
    rgrain2(ind_z_surf-1) = rgrain2(ind_z_surf)

    ! Remove old surface layer
    T(ind_z_surf) = 0.
    M(ind_z_surf) = 0.
    DZ(ind_z_surf) = 0.
    Rho(ind_z_surf) = 0.
    DenRho(ind_z_surf) = 0.
    Mlwc(ind_z_surf) = 0.
    Refreeze(ind_z_surf) = 0.
    Year(ind_z_surf) = 0.
    rgrain2(ind_z_surf) = 0.

    ind_z_surf = ind_z_surf-1

    end subroutine Merge_Second_Layer_By_Threshold

! *******************************************************
  subroutine Remove_Surface_Layer(ind_z_max, ind_z_surf, Rho, M, T, Mlwc, DZ, DenRho, Refreeze, Year, rgrain2)
    !*** removes surface layer when more mass is removed than present in the surface layer

    ! declare arguments
    integer, intent(in) :: ind_z_max
    integer, intent(inout) :: ind_z_surf
    double precision, dimension(ind_z_max), intent(inout) :: Rho, M, T, Mlwc, DZ, DenRho, Refreeze, Year, rgrain2

    ! Add liquid water from surface layer into new surface layer
    Mlwc(ind_z_surf-1) = Mlwc(ind_z_surf) + Mlwc(ind_z_surf - 1)
    
    ! Remove old surface layer
    T(ind_z_surf) = 0.
    M(ind_z_surf) = 0.
    DZ(ind_z_surf) = 0.
    Rho(ind_z_surf) = 0.
    DenRho(ind_z_surf) = 0.
    Mlwc(ind_z_surf) = 0.
    Refreeze(ind_z_surf) = 0.
    Year(ind_z_surf) = 0.
    rgrain2(ind_z_surf) = 0.

    ind_z_surf = ind_z_surf-1

    end subroutine Remove_Surface_Layer

end module grid_routines
