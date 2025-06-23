module time_loop
    !*** Subroutines for stepping through time ***!

    use grid_routines, only: Split_Layers, Merge_Layers, Delete_Layers, Add_Layers
    use firn_physics, only: Update_Surface, Densific, Solve_Temp_Imp, Solve_Temp_Exp, Calc_Integrated_Var
    use output, only: Accumulate_Output, To_out_1D, To_out_2D, To_out_2Ddetail
    use model_settings

    implicit none
    ! explicitly set array offset
    integer, private :: array_offset = 150

    public :: Time_Loop_SpinUp, Time_Loop_Main
    
contains


! *******************************************************


subroutine Time_Loop_SpinUp(Nt_model_tot, Nt_model_spinup, ind_z_max, ind_z_surf, dtmodel, R, Ec, Eg, g, Lh, rhoi, acav, th, dzmax, M, T, DZ, Rho, &
    DenRho, Depth, Mlwc, Refreeze, Year, TempFM, PSolFM, PLiqFM, SublFM, MeltFM, DrifFM, Rho0FM, IceShelf, &
    ImpExp, nyears, nyearsSU)
    !*** Subroutine for repeatedly repeating the spin-up until a steady state is reached ***!
        
    ! declare arguments
    integer, intent(in) :: Nt_model_tot, Nt_model_spinup, ind_z_max, dtmodel, nyears, nyearsSU, IceShelf, ImpExp
    integer, intent(inout) :: ind_z_surf
    double precision, intent(in) :: R, Ec, Eg, g, Lh, rhoi, acav, th, dzmax
    double precision, dimension(ind_z_max), intent(inout) :: Rho, M, T, Depth, Mlwc, DZ, DenRho, Refreeze, Year
    double precision, dimension(Nt_model_tot), intent(in) :: TempFM, PSolFM, PLiqFM, SublFM
    double precision, dimension(Nt_model_tot), intent(in) :: MeltFM, DrifFM, Rho0FM

    ! declare local variables
    integer :: spinup_numb, ind_z, ind_t
    double precision :: rho0, Ts, Psol, Pliq, Su, Me, Sd
    double precision :: h_surf, vice, vmelt, vacc, vsub, vsnd, vfc, vbouy, FirnAir, TotLwc, IceMass
    double precision :: Mrain = 0., Mrunoff = 0., Mrefreeze = 0., Msurfmelt = 0., Msolin = 0.
    double precision :: z_surf_old, fac_old, z_surf_error, fac_error, spinup_bound, error_bound

    spinup_numb = 0
    h_surf = 1.E3
    FirnAir = 1.E3
    z_surf_error = 1.E3
    fac_error = 1.E3

    ! do at least 3 spin-ups
    ! never do more than 70 spin-ups for Greenland or 200 for Antarctica
    ! stop when the elevation at the end changes less than 1 cm (for Greenland and 2 cm for Antarctica)
    ! and FAC changes less than 1 cm.
    if (domain .eq. "FGRN055") then
            spinup_bound = 70
            error_bound = 0.0001
    elseif (domain .eq. "ANT27") then 
           spinup_bound = 200  !perhaps this one can be reduced 
           error_bound = 0.004 
    else
            print *, "No spinup bounds available for domain: ", domain
    end if 

    do while ( ( ( (z_surf_error > error_bound) .or. (fac_error > error_bound) ) .and. (spinup_numb < spinup_bound) ) .or. (spinup_numb < 3) )   
        spinup_numb = spinup_numb + 1
        z_surf_old = h_surf
        fac_old = FirnAir
        
        h_surf = 0.
        FirnAir = 0.
        TotLwc = 0.
        IceMass = 0.

        !print *, 'After spin-up #, spinup_numb, Rho(200), T(200), Year(200), ind_z_surf, h_surf, FirnAir, IceMass'

        do ind_t = 1, Nt_model_spinup 
        
            ! Index the forcing at the current time step
            rho0 = Rho0FM(ind_t)
            Ts = TempFM(ind_t)
            Psol = PSolFM(ind_t)
            Pliq = PliqFM(ind_t)
            Su = SublFM(ind_t)
            Me = MeltFM(ind_t)
            Sd = DrifFM(ind_t)
            
            ! Calculate the densification	  
            call Densific(ind_z_max, ind_z_surf, dtmodel, R, Ec, Eg, g, rhoi, acav, Rho, T, domain)
            
            ! Re-calculate the Temp-profile (explicit or implicit)		  
            if (ImpExp == 1) call Solve_Temp_Imp(ind_z_max, ind_z_surf, dtmodel, th, Ts, T, Rho, DZ, rhoi)
            if (ImpExp == 2) call Solve_Temp_Exp(ind_z_max, ind_z_surf, dtmodel, Ts, T, Rho, DZ, rhoi)

            ! Re-caluclate DZ/M-values according to new Rho-/T-values
            call Update_Surface(ind_z_max, ind_z_surf, dtmodel, rho0, rhoi, acav, Lh, h_surf, vice, vmelt, vacc, vsub, &
                vsnd, vfc, vbouy, Ts, PSol, PLiq, Su, Me, Sd, M, T, DZ, Rho, DenRho, Mlwc,Refreeze, &
                ImpExp, IceShelf, Msurfmelt, Mrain, Msolin, Mrunoff, Mrefreeze)

            ! Check if the vertical grid is still valid
            if (DZ(ind_z_surf) > dzMAX) then
                call Split_Layers(ind_z_max, ind_z_surf, Rho, M, T, Mlwc, DZ, DenRho, Refreeze, Year, ind_t, Nt_model_tot, nyears)
            endif
            if (DZ(ind_z_surf) <= dzMAX/3.) then
                call Merge_Layers(ind_z_max, ind_z_surf, Rho, M, T, Mlwc, DZ, DenRho, Refreeze, Year)
            endif
            if ((ind_z_surf > 2500) .and. (MINVAL(Rho(1:200)) >= (rhoi-7.))) then
                call Delete_Layers(ind_z_max, ind_z_surf, Rho, M, T, Mlwc, DZ, DenRho, Refreeze, Year)
            endif
            if (ind_z_surf < 200) then
                call Add_Layers(ind_z_max, ind_z_surf, dzmax, rhoi, Rho, M, T, Mlwc, DZ, DenRho, Refreeze, Year)
            endif
    
            ! Update the depth of each layer
            Depth(ind_z_surf) = DZ(ind_z_surf)/2.
            do ind_z = (ind_z_surf-1), 1, -1
                Depth(ind_z) = Depth(ind_z+1) + DZ(ind_z+1)/2. + DZ(ind_z)/2.
            enddo
            
            ! Reset variables that accumulated during the time step (not used in spin-up)
            Msurfmelt = 0.
            Mrain = 0.
            Msolin = 0.
            Mrunoff = 0.
            Mrefreeze = 0.

        enddo  ! time loop

        Year(:) = Year(:) - DBLE(nyearsSU)
        
        ! Calculate the firn air content, ice mass and total liquid water content of the firn column
        call Calc_Integrated_Var(ind_z_max, ind_z_surf, rhoi, Rho, Mlwc, M, DZ, FirnAir, TotLwc, IceMass)
        
        print *, "After spin-up #", spinup_numb
        print *, "Rho(200) = ", Rho(200), "T(200) = ", T(200), "Year(200) = ", Year(200)
        print *, "Ind_z_surf = ", ind_z_surf, "h_surf = ", h_surf, "FAC = ", FirnAir, "IceMass = ", IceMass

        z_surf_error = (h_surf-z_surf_old)*(h_surf-z_surf_old)
        fac_error = (FirnAir-fac_old)*(FirnAir-fac_old)

    enddo  ! spin-ups

    ! Reset profiles that need to start at zero after the spin-up
    DenRho(:) = 0.
    Refreeze(:) = 0.
    
end subroutine Time_Loop_SpinUp


! *******************************************************


subroutine Time_Loop_Main(dtmodel, ImpExp, Nt_model_tot, nyears, ind_z_max, ind_z_surf, numOutputSpeed, numOutputProf, numOutputDetail, &
    outputSpeed, outputProf, outputDetail, th, R, Ec, Eg, g, Lh, dzmax, rhoi, proflayers, detlayers, detthick, acav, IceShelf, &
    TempFM, PsolFM, PliqFM, SublFM, MeltFM, DrifFM, Rho0FM, Rho, M, T, Depth, Mlwc, DZ, DenRho, Refreeze, Year, &
    domain, out_1D, out_2D_dens, out_2D_temp, out_2D_lwc, out_2D_depth, out_2D_dRho, out_2D_year, &
    out_2D_det_dens, out_2D_det_temp, out_2D_det_lwc, out_2D_det_refreeze, prev_nt, restart_type)
    !*** Subrouting for stepping through time after the spin-up is complete, meanwhile writing output to netcdf ***!
    
    ! declare arguments
    integer, intent(in) :: dtmodel, ImpExp, IceShelf, ind_z_max, Nt_model_tot, nyears, proflayers, detlayers, prev_nt
    integer, intent(in) :: numOutputSpeed, numOutputProf, numOutputDetail
    integer, intent(in) :: outputSpeed, outputProf, outputDetail
    integer, intent(inout) :: ind_z_surf
    double precision, intent(in) :: th, R, Ec, Eg, g, Lh, dzmax, rhoi, acav, detthick
    double precision,dimension(Nt_model_tot), intent(in) :: TempFM, PsolFM, PliqFM, SublFM, MeltFM, DrifFM, Rho0FM
    double precision,dimension(ind_z_max), intent(inout) :: Rho, M, T, Depth, Mlwc, DZ, DenRho, Refreeze, Year
    double precision, dimension((outputSpeed+array_offset), 18), intent(inout) :: out_1D
    double precision, dimension((outputProf+array_offset), proflayers), intent(inout) :: out_2D_dens, out_2D_temp, out_2D_lwc, out_2D_depth, out_2D_dRho, out_2D_year
    double precision, dimension((outputDetail+array_offset), detlayers), intent(inout) :: out_2D_det_dens, out_2D_det_temp, out_2D_det_lwc, out_2D_det_refreeze
    character*255, intent(in) :: domain, restart_type
    
    ! declare local variables
    integer :: ind_z, ind_t, ind_t_i
    double precision :: rho0, Ts, Psol, Pliq, Su, Me, Sd
    double precision :: vice, vmelt, vacc, vsub, vsnd, vfc, vbouy, FirnAir, TotLwc, IceMass
    double precision :: h_surf = 0., Totvice = 0., Totvfc = 0., Totvacc = 0., Totvsub = 0., Totvsnd = 0., Totvmelt = 0., Totvbouy = 0.
    double precision :: Mrunoff = 0., TotRunoff = 0., Mrefreeze = 0., TotRefreeze = 0.
    double precision :: Mrain = 0., TotRain = 0., Msurfmelt = 0., TotSurfmelt = 0., Msolin = 0., TotSolIn = 0., Rho0out = 0.
    
    ! Time integration
    print *, "Start of main time loop"

    ! start from last time step if restarting from loaded run
    if ( restart_type=="run" ) then
        ind_t_i = prev_nt
    else
        ind_t_i = 1
    endif

    do ind_t = ind_t_i, Nt_model_tot
        
        ! Index the forcing at the current time step
        rho0 = Rho0FM(ind_t)
        Ts = TempFM(ind_t)
        Psol = PSolFM(ind_t)
        Pliq = PliqFM(ind_t)
        Su = SublFM(ind_t)
        Me = MeltFM(ind_t)
        Sd = DrifFM(ind_t)
        
        ! Calculate the density profile      
        call Densific(ind_z_max, ind_z_surf, dtmodel, R, Ec, Eg, g, rhoi, acav, Rho, T, domain)
        
        ! Calculate the temperature profile (explicit or implicit)         
        if (ImpExp == 1) call Solve_Temp_Imp(ind_z_max, ind_z_surf, dtmodel, th, Ts, T, Rho, DZ, rhoi)
        if (ImpExp == 2) call Solve_Temp_Exp(ind_z_max, ind_z_surf, dtmodel, Ts, T, Rho, DZ, rhoi)
                
        ! Add/remove mass from the surface layer and the update the layer thickness
        call Update_Surface(ind_z_max, ind_z_surf, dtmodel, rho0, rhoi, acav, Lh, h_surf, vice, vmelt, vacc, vsub, &
            vsnd, vfc, vbouy, Ts, PSol, PLiq, Su, Me, Sd, M, T, DZ, Rho, DenRho, Mlwc,Refreeze, &
            ImpExp, IceShelf, Msurfmelt, Mrain, Msolin, Mrunoff, Mrefreeze)
        
        if (mod(ind_t, 200000) == 0) print *, ind_t, h_surf
        
        ! Check if the vertical grid is still valid
        if (DZ(ind_z_surf) > dzMAX) then
            call Split_Layers(ind_z_max, ind_z_surf, Rho, M, T, Mlwc, DZ, DenRho, Refreeze, Year, ind_t, Nt_model_tot, nyears)
        endif
        if (DZ(ind_z_surf) <= dzMAX/3.) then
            call Merge_Layers(ind_z_max, ind_z_surf, Rho, M, T, Mlwc, DZ, DenRho, Refreeze, Year)
        endif
        if ((ind_z_surf > 2500) .and. (MINVAL(Rho(1:200)) >= (rhoi-7.))) then
            call Delete_Layers(ind_z_max, ind_z_surf, Rho, M, T, Mlwc, DZ, DenRho, Refreeze, Year)
        endif
        if (ind_z_surf < 200) then
            call Add_Layers(ind_z_max, ind_z_surf, dzmax, rhoi, Rho, M, T, Mlwc, DZ, DenRho, Refreeze, Year)
        endif

        ! Update the depth of each layer
        Depth(ind_z_surf) = DZ(ind_z_surf)/2.
        do ind_z = (ind_z_surf-1), 1, -1
            Depth(ind_z) = Depth(ind_z+1) + DZ(ind_z+1)/2. + DZ(ind_z)/2.
        enddo

        ! Calculate the firn air content, ice mass and total liquid water content of the firn column
        call Calc_Integrated_Var(ind_z_max, ind_z_surf, rhoi, Rho, Mlwc, M, DZ, FirnAir, TotLwc, IceMass)
        
        call Accumulate_Output(dtmodel, vice, vmelt, vacc, vsub, vsnd, vfc, vbouy, Totvice, Totvacc, Totvsub, &
            Totvsnd, Totvfc, Totvmelt, Totvbouy, Mrunoff, TotRunoff, Mrefreeze, TotRefreeze, Mrain, TotRain, &
            Msurfmelt, TotSurfmelt, Msolin, TotSolIn, rho0, Rho0out)
        
        ! Save output data, if needed
        
        if (mod(ind_t, numOutputSpeed) == 0) then
            call To_out_1D(ind_t, numOutputSpeed, h_surf, Totvice, Totvfc, Totvacc, Totvsub, Totvsnd, Totvmelt, &
                Totvbouy, TotRunoff, FirnAir, TotLwc, TotRefreeze, TotRain, TotSurfmelt, TotSolIn, IceMass, Rho0out, &
                out_1D, outputSpeed)
        endif
        
        if (mod(ind_t, numOutputProf) == 0.) then
            call To_out_2D(ind_z_max, ind_z_surf, ind_t, dtmodel, numOutputProf, outputProf, proflayers, &
                Rho, T, Mlwc, Depth, DenRho, Year, out_2D_dens, out_2D_temp, out_2D_lwc, &
                out_2D_depth, out_2D_dRho, out_2D_year)
        endif
        
        if (mod(ind_t, numOutputDetail) == 0.) then
            call To_out_2Ddetail(ind_z_max, ind_z_surf, ind_t, detlayers, detthick, numOutputDetail, &
                outputDetail, Rho, T, Mlwc, Refreeze, DZ, out_2D_det_dens, out_2D_det_temp, &
                out_2D_det_lwc, out_2D_det_refreeze)
        endif
    enddo
    
    ! Finished time loop
    print *, "End of time loop"
    print *, " "
        
end subroutine Time_Loop_Main


end module time_loop
