module time_loop
    !*** Subroutines for stepping through time ***!

    use grid_routines, only: Split_Layers, Merge_Layers, Delete_Layers, Add_Layers
    use firn_physics, only: Update_Surface, Densific, Solve_Temp_Imp, Solve_Temp_Exp, Calc_Integrated_Var
    use output, only: Accumulate_Output, To_out_1D, To_out_2D, To_out_2Ddetail
    use model_settings

    implicit none

    public :: Time_Loop_SpinUp, Time_Loop_Main
    
contains


! *******************************************************


subroutine Time_Loop_SpinUp(Nt_model_tot, Nt_model_spinup, dtmodel, acav, ffav, M, T, DZ, Rho, &
    DenRho, Depth, Mlwc, Refreeze, Year, TempFM, PSolFM, PLiqFM, SublFM, MeltFM, DrifFM, Rho0FM, IceShelf)
    !*** Subroutine for repeatedly repeating the spin-up until a steady state is reached ***!

    ! declare arguments
    integer, intent(in) :: Nt_model_tot, Nt_model_spinup, dtmodel, IceShelf
    double precision, intent(in) :: acav, ffav
    double precision, dimension(ind_z_max), intent(inout) :: Rho, M, T, Depth, Mlwc, DZ, DenRho, Refreeze, Year
    double precision, dimension(:), intent(in) :: TempFM, PSolFM, PLiqFM, SublFM
    double precision, dimension(:), intent(in) :: MeltFM, DrifFM, Rho0FM

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
            write(log_unit, *) "No spinup bounds available for domain: ", domain
    end if 

    do while ( ( ( (z_surf_error > error_bound) .or. (fac_error > error_bound) ) .and. (spinup_numb < spinup_bound) ) .or. (spinup_numb < 3) )   
        spinup_numb = spinup_numb + 1
        z_surf_old = h_surf
        fac_old = FirnAir
        
        h_surf = 0.
        FirnAir = 0.
        TotLwc = 0.
        IceMass = 0.

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
            call Densific(dtmodel, acav, ffav, Rho, T, ind_t)

            ! Re-calculate the Temp-profile (explicit or implicit)
            if (config%model_choices%ImpExp == 1) call Solve_Temp_Imp(dtmodel, Ts, T, Rho, DZ)
            if (config%model_choices%ImpExp == 2) call Solve_Temp_Exp(dtmodel, Ts, T, Rho, DZ)

            ! Re-caluclate DZ/M-values according to new Rho-/T-values
            call Update_Surface(dtmodel, rho0, acav, h_surf, vice, vmelt, vacc, vsub, &
                vsnd, vfc, vbouy, Ts, PSol, PLiq, Su, Me, Sd, M, T, DZ, Rho, DenRho, Mlwc,Refreeze, &
                IceShelf, Msurfmelt, Mrain, Msolin, Mrunoff, Mrefreeze)

            ! Check if the vertical grid is still valid
            if (DZ(ind_z_surf) > config%model_choices%dzmax) then
                call Split_Layers(Rho, M, T, Mlwc, DZ, DenRho, Refreeze, Year, ind_t, Nt_model_tot)
            endif
            if (DZ(ind_z_surf) <= config%model_choices%dzmax/3.) then
                call Merge_Layers(Rho, M, T, Mlwc, DZ, DenRho, Refreeze, Year)
            endif
            if ((ind_z_surf > 2500) .and. (MINVAL(Rho(1:200)) >= (const%rhoi-7.))) then
                call Delete_Layers(Rho, M, T, Mlwc, DZ, DenRho, Refreeze, Year)
            endif
            if (ind_z_surf < 200) then
                call Add_Layers(Rho, M, T, Mlwc, DZ, DenRho, Refreeze, Year)
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

        Year(:) = Year(:) - DBLE(config%forcing_dimensions%nyears_spinup)
        
        ! Calculate the firn air content, ice mass and total liquid water content of the firn column
        call Calc_Integrated_Var(Rho, Mlwc, M, DZ, FirnAir, TotLwc, IceMass)
        
        write(log_unit, *) "After spin-up #", spinup_numb
        write(log_unit, *) "Rho(200) = ", Rho(200), "T(200) = ", T(200), "Year(200) = ", Year(200)
        write(log_unit, *) "Ind_z_surf = ", ind_z_surf, "h_surf = ", h_surf, "FAC = ", FirnAir, "IceMass = ", IceMass

        z_surf_error = (h_surf-z_surf_old)*(h_surf-z_surf_old)
        fac_error = (FirnAir-fac_old)*(FirnAir-fac_old)

    enddo  ! spin-ups

    ! Reset profiles that need to start at zero after the spin-up
    DenRho(:) = 0.
    Refreeze(:) = 0.
    
end subroutine Time_Loop_SpinUp


! *******************************************************


subroutine Time_Loop_Main(dtmodel, Nt_model_tot, numOutputSpeed, numOutputProf, numOutputDetail, &
    outputSpeed, outputProf, outputDetail, acav, ffav, IceShelf, &
    TempFM, PsolFM, PliqFM, SublFM, MeltFM, DrifFM, Rho0FM, Rho, M, T, Depth, Mlwc, DZ, DenRho, Refreeze, Year, &
    out_1D, out_2D_dens, out_2D_temp, out_2D_lwc, out_2D_depth, out_2D_dRho, out_2D_year, &
    out_2D_det_dens, out_2D_det_temp, out_2D_det_lwc, out_2D_det_refreeze, prev_nt)
    !*** Subrouting for stepping through time after the spin-up is complete, meanwhile writing output to netcdf ***!

    ! declare arguments
    integer, intent(in) :: dtmodel, IceShelf, Nt_model_tot, prev_nt
    integer, intent(in) :: numOutputSpeed, numOutputProf, numOutputDetail
    integer, intent(in) :: outputSpeed, outputProf, outputDetail
    double precision, intent(in) :: acav, ffav
    double precision, dimension(:), intent(in) :: TempFM, PsolFM, PliqFM, SublFM, MeltFM, DrifFM, Rho0FM
    double precision,dimension(ind_z_max), intent(inout) :: Rho, M, T, Depth, Mlwc, DZ, DenRho, Refreeze, Year
    double precision, dimension(:,:), intent(inout) :: out_1D
    double precision, dimension(:,:), intent(inout) :: out_2D_dens, out_2D_temp, out_2D_lwc, out_2D_depth, out_2D_dRho, out_2D_year
    double precision, dimension(:,:), intent(inout) :: out_2D_det_dens, out_2D_det_temp, out_2D_det_lwc, out_2D_det_refreeze
    
    ! declare local variables
    integer :: ind_z, ind_t, ind_t_i
    double precision :: rho0, Ts, Psol, Pliq, Su, Me, Sd
    double precision :: vice, vmelt, vacc, vsub, vsnd, vfc, vbouy, FirnAir, TotLwc, IceMass
    double precision :: h_surf = 0., Totvice = 0., Totvfc = 0., Totvacc = 0., Totvsub = 0., Totvsnd = 0., Totvmelt = 0., Totvbouy = 0.
    double precision :: Mrunoff = 0., TotRunoff = 0., Mrefreeze = 0., TotRefreeze = 0.
    double precision :: Mrain = 0., TotRain = 0., Msurfmelt = 0., TotSurfmelt = 0., Msolin = 0., TotSolIn = 0., Rho0out = 0.
    
    ! Time integration
    write(log_unit, *) "Start of main time loop"

    ! start from last time step if restarting from loaded run
    if (restart_type == "run") then
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
        call Densific(dtmodel, acav, ffav, Rho, T, ind_t)

        ! Calculate the temperature profile (explicit or implicit)
        if (config%model_choices%ImpExp == 1) call Solve_Temp_Imp(dtmodel, Ts, T, Rho, DZ)
        if (config%model_choices%ImpExp == 2) call Solve_Temp_Exp(dtmodel, Ts, T, Rho, DZ)

        ! Add/remove mass from the surface layer and the update the layer thickness
        call Update_Surface(dtmodel, rho0, acav, h_surf, vice, vmelt, vacc, vsub, &
            vsnd, vfc, vbouy, Ts, PSol, PLiq, Su, Me, Sd, M, T, DZ, Rho, DenRho, Mlwc,Refreeze, &
            IceShelf, Msurfmelt, Mrain, Msolin, Mrunoff, Mrefreeze)
        
        if (mod(ind_t, 200000) == 0) write(log_unit, *) ind_t, h_surf
        
        ! Check if the vertical grid is still valid
        if (DZ(ind_z_surf) > config%model_choices%dzmax) then
            call Split_Layers(Rho, M, T, Mlwc, DZ, DenRho, Refreeze, Year, ind_t, Nt_model_tot)
        endif
        if (DZ(ind_z_surf) <= config%model_choices%dzmax/3.) then
            call Merge_Layers(Rho, M, T, Mlwc, DZ, DenRho, Refreeze, Year)
        endif

        if ((ind_z_surf > 2500) .and. (MINVAL(Rho(1:200)) >= (const%rhoi-7.))) then
            call Delete_Layers(Rho, M, T, Mlwc, DZ, DenRho, Refreeze, Year)
        endif
        if (ind_z_surf < 200) then
            call Add_Layers(Rho, M, T, Mlwc, DZ, DenRho, Refreeze, Year)
        endif

        ! Update the depth of each layer
        Depth(ind_z_surf) = DZ(ind_z_surf)/2.
        do ind_z = (ind_z_surf-1), 1, -1
            Depth(ind_z) = Depth(ind_z+1) + DZ(ind_z+1)/2. + DZ(ind_z)/2.
        enddo

        ! Calculate the firn air content, ice mass and total liquid water content of the firn column
        call Calc_Integrated_Var(Rho, Mlwc, M, DZ, FirnAir, TotLwc, IceMass)
        
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
            call To_out_2D(ind_t, dtmodel, numOutputProf, outputProf, &
                Rho, T, Mlwc, Depth, DenRho, Year, out_2D_dens, out_2D_temp, out_2D_lwc, &
                out_2D_depth, out_2D_dRho, out_2D_year)
        endif

        if (mod(ind_t, numOutputDetail) == 0.) then
            call To_out_2Ddetail(ind_t, numOutputDetail, &
                outputDetail, Rho, T, Mlwc, Refreeze, DZ, out_2D_det_dens, out_2D_det_temp, &
                out_2D_det_lwc, out_2D_det_refreeze)
        endif
    enddo
    
    ! Finished time loop
    write(log_unit, *) "End of time loop"
    write(log_unit, *) " "
        
end subroutine Time_Loop_Main


end module time_loop
