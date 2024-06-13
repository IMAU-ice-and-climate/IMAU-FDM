module output
    !*** functions and subroutines for writing data to netcdf files

    use netcdf, only: nf90_create, nf90_def_dim, nf90_def_var, nf90_real, nf90_int, nf90_noerr, &
        nf90_enddef, nf90_put_var, nf90_close, nf90_unlimited, nf90_put_att
    use openNetCDF, only: Handle_Error

    implicit none
    private

    public :: Accumulate_Output, Write_Initial_Output, To_out_1D, To_out_2D, To_out_2Ddetail, Save_out_1D, Save_out_2D, Save_out_2Ddetail, Save_out_restart
    
contains


! *******************************************************


subroutine Accumulate_Output(dtmodel, vice, vmelt, vacc, vsub, vsnd, vfc, vbouy, Totvice, Totvacc, Totvsub, &
    Totvsnd, Totvfc, Totvmelt, Totvbouy, Mrunoff, TotRunoff, Mrefreeze, TotRefreeze, Mrain, TotRain, &
    Msurfmelt, TotSurfmelt, Msolin, TotSolIn, rho0, Rho0out)
    !*** Add up the 1D output variables that are averaged / accumulated over the output time period ***!

    ! declare arguments
    integer, intent(in) :: dtmodel
    double precision, intent(in) :: vice, vmelt, vacc, vsub, vsnd, vfc, vbouy
    double precision, intent(inout) :: Totvice, Totvfc, Totvacc, Totvsub, Totvsnd, Totvmelt, Totvbouy
    double precision, intent(inout) :: Mrunoff, TotRunoff, Mrefreeze, TotRefreeze
    double precision, intent(inout) :: Mrain, TotRain, Msurfmelt, TotSurfmelt, Msolin, TotSolIn, rho0, Rho0out

    ! All velocities in m/s
    Totvice = Totvice + (vice/dtmodel)
    Totvacc = Totvacc + (vacc/dtmodel)
    Totvsub = Totvsub + (vsub/dtmodel)
    Totvsnd = Totvsnd + (vsnd/dtmodel)
    Totvfc = Totvfc + (vfc/dtmodel)
    Totvmelt = Totvmelt + (vmelt/dtmodel)
    Totvbouy = Totvbouy + (vbouy/dtmodel)
    
    ! These variables are accumulated over time
    TotRunoff = TotRunoff + Mrunoff
    TotRefreeze = TotRefreeze + Mrefreeze
    TotRain = TotRain + Mrain
    TotSurfMelt = TotSurfMelt + Msurfmelt
    TotSolIn = TotSolIn + Msolin
    Rho0out = Rho0out + rho0

    ! Reset to 0 for the next time step
    Mrunoff = 0.
    Mrefreeze = 0.
    Mrain = 0.
    Msurfmelt = 0.
    Msolin = 0.
    
end subroutine Accumulate_Output


! *******************************************************


subroutine Write_Initial_Output(ind_z_max, ind_z_surf, Rho, M, T, Depth, Mlwc, Year, path_out_ini, fname_out_ini)
    !*** Write a netcdf file containing the firn profile after the spin-up ***!
    
    ! declare arguments
    integer, intent(in) :: ind_z_max, ind_z_surf
    double precision, dimension(ind_z_max), intent(in) :: Rho, M, T, Depth, Mlwc, Year
    character*255, intent(in) :: path_out_ini, fname_out_ini

    ! declare local arguments
    integer :: status, ncid, dimid_z, varID(10)

    ! CREATE NETCDF FILES
    status = nf90_create(trim(path_out_ini)//trim(fname_out_ini), 0, ncid)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_create')

    ! DEFINE DIMENSIONS
    status = nf90_def_dim(ncid, "layer", ind_z_surf, dimid_z)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_dim1')
    
    ! DEFINE VARIABLES
    status = nf90_def_var(ncid, "dens", nf90_real, dimid_z, varID(1))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var1')
    status = nf90_def_var(ncid, "temp", nf90_real, dimid_z, varID(2))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var2')
    status = nf90_def_var(ncid, "mass", nf90_real, dimid_z, varID(3))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var3')
    status = nf90_def_var(ncid, "depth", nf90_real, dimid_z, varID(4))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var4')
    status = nf90_def_var(ncid, "lwc", nf90_real, dimid_z, varID(5))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var5')
    status = nf90_def_var(ncid, "year", nf90_real, dimid_z, varID(6))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var6')
    
    ! END OF DEFINING FILES
    status = nf90_enddef(ncid)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_enddef')
    
    ! PUT VARIABLES
    status = nf90_put_var(ncid, varID(1), Rho(1:ind_z_surf))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var1')
    status = nf90_put_var(ncid, varID(2), T(1:ind_z_surf))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var2')
    status = nf90_put_var(ncid, varID(3), M(1:ind_z_surf))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var3')
    status = nf90_put_var(ncid, varID(4), Depth(1:ind_z_surf))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var4')
    status = nf90_put_var(ncid, varID(5), Mlwc(1:ind_z_surf))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var5')
        
    status = nf90_put_var(ncid, varID(6), Year(1:ind_z_surf))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var6')
    
    ! CLOSE NETCDF FILE
    status = nf90_close(ncid)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_close')
    
end subroutine Write_Initial_Output


! *******************************************************


subroutine To_out_1D(ind_t, numOutputSpeed, h_surf, Totvice, Totvfc, Totvacc, Totvsub, Totvsnd, Totvmelt, &
    Totvbouy, TotRunoff, FirnAir, TotLwc, TotRefreeze, TotRain, TotSurfmelt, TotSolIn, IceMass, Rho0out, &
    out_1D, outputSpeed)
    !*** Write the 1D output variables to the variable that will be converted into a netcdf file after the time loop ***!

    ! declare arguments
    integer, intent(in) :: ind_t, numOutputSpeed, outputSpeed
    double precision, intent(in) :: h_surf, FirnAir, TotLwc
    double precision, intent(inout) :: Totvice, Totvfc, Totvacc, Totvsub, Totvsnd, Totvmelt, Totvbouy
    double precision, intent(inout) :: TotRunoff, TotRefreeze, TotRain, TotSurfmelt, TotSolIn, IceMass, Rho0out
    double precision, dimension((outputSpeed+50),18), intent(out) :: out_1D

    ! declare local variables
    integer :: ind_t_out
    double precision :: factor, Totv

    ! All velocities in m/year
    ! Divide all variables by the number of time steps within the output period
    factor = (3600.*24.*365.)/numOutputSpeed
    Totvice = -1. * (Totvice * factor)
    Totvacc = Totvacc * factor
    Totvsub = Totvsub * factor
    Totvsnd = Totvsnd * factor
    Totvmelt = Totvmelt * factor
    Totvfc = Totvfc * factor
    Totvbouy = Totvbouy * factor
    Totv = Totvacc + Totvice + Totvmelt + Totvfc + Totvbouy
    Rho0out = Rho0out / numOutputSpeed
    
    ind_t_out = ind_t/numOutputSpeed

    out_1D(ind_t_out,1) = REAL(h_surf)
    out_1D(ind_t_out,2) = REAL(Totvice)
    out_1D(ind_t_out,3) = REAL(Totvacc)
    out_1D(ind_t_out,4) = REAL(Totvfc)
    out_1D(ind_t_out,5) = REAL(Totvmelt)
    out_1D(ind_t_out,6) = REAL(Totvbouy)
    out_1D(ind_t_out,7) = REAL(Totvsub)
    out_1D(ind_t_out,8) = REAL(Totvsnd)
    out_1D(ind_t_out,9) = REAL(Totv)
    out_1D(ind_t_out,10) = REAL(TotRunoff)
    out_1D(ind_t_out,11) = REAL(FirnAir)
    out_1D(ind_t_out,12) = REAL(TotLwc)
    out_1D(ind_t_out,13) = REAL(TotRefreeze)
    out_1D(ind_t_out,14) = REAL(TotRain)
    out_1D(ind_t_out,15) = REAL(TotSurfmelt)
    out_1D(ind_t_out,16) = REAL(TotSolIn)
    out_1D(ind_t_out,17) = REAL(IceMass)
    out_1D(ind_t_out,18) = REAL(Rho0out)
    
    Totvice = 0.
    Totvacc = 0.
    Totvsub = 0.
    Totvsnd = 0.
    Totvmelt = 0.
    Totvfc = 0.
    Totvbouy = 0.
    Totv = 0.
    TotRunoff = 0.
    TotRefreeze = 0.
    TotRain = 0.
    TotSurfmelt = 0.
    TotSolIn = 0.
    Rho0out = 0.
    
end subroutine To_out_1D


! *******************************************************


subroutine To_out_2D(ind_z_max, ind_z_surf, ind_t, dtmodel, numOutputProf, outputProf, proflayers, Rho, &
    T, rgrain2, Mlwc, Depth, DenRho, Year, out_2D_dens, out_2D_temp, out_2D_lwc, out_2D_depth, out_2D_dRho, &
    out_2D_year, out_2D_rgrain)
    !*** Write the 2D output variables to the variables that will be converted into a netcdf file after the time loop ***!
    
    ! declare arguments
    integer, intent(in) :: ind_z_max, ind_z_surf, ind_t, dtmodel, numOutputProf, outputProf, proflayers
    double precision, dimension(ind_z_max), intent(in) :: Rho, T, rgrain2, Mlwc, Depth, Year
    double precision, dimension(ind_z_max), intent(inout) :: DenRho
    double precision, dimension(outputProf+50,proflayers), intent(out) :: out_2D_dens, out_2D_temp, out_2D_lwc, &
        out_2D_depth, out_2D_dRho, out_2D_year, out_2D_rgrain

    ! declare local variables
    integer :: ind_t_out, ind_z_bot, ind_z_surf_out
    double precision :: factor
    
    ! Determine indices for the output
    ind_t_out = ind_t/numOutputProf
    ind_z_bot = 1
    ind_z_surf_out = ind_z_surf
    if (ind_z_surf > proflayers) then 
        ind_z_bot = ind_z_surf - proflayers + 1
        ind_z_surf_out = proflayers
    endif
    
    ! Convert to [kg m-3 year-1]
    factor = (3600.*24.*365.) / dtmodel / numOutputProf
    DenRho(:) = DenRho(:) * factor
    
    ! save output to 2D array
    out_2D_dens(ind_t_out, 1:ind_z_surf_out) = REAL(Rho(ind_z_bot:ind_z_surf))
    out_2D_temp(ind_t_out, 1:ind_z_surf_out) = REAL(T(ind_z_bot:ind_z_surf))
    out_2D_lwc(ind_t_out, 1:ind_z_surf_out) = REAL(Mlwc(ind_z_bot:ind_z_surf))
    out_2D_depth(ind_t_out, 1:ind_z_surf_out) = REAL(Depth(ind_z_bot:ind_z_surf))
    out_2D_dRho(ind_t_out, 1:ind_z_surf_out) = REAL(DenRho(ind_z_bot:ind_z_surf))
    out_2D_year(ind_t_out, 1:ind_z_surf_out) = REAL(Year(ind_z_bot:ind_z_surf))
    out_2D_rgrain(ind_t_out, 1:ind_z_surf_out) = REAL(1000.*(rgrain2(ind_z_bot:ind_z_surf)**0.5))
    
    DenRho(:) = 0.
    
end subroutine To_out_2D


! *******************************************************


subroutine To_out_2Ddetail(ind_z_max, ind_z_surf, ind_t,detlayers, detthick, numOutputDetail, outputDetail, &
    Rho, T, rgrain2, Mlwc, Refreeze, DZ, out_2D_det_dens, out_2D_det_temp, out_2D_det_lwc, out_2D_det_refreeze, out_2D_det_rgrain)
    !*** Write the 2Ddetail output variables to the variables that will be converted into a netcdf file after the time loop ***!
    
    ! declare arguments
    integer, intent(in) :: ind_z_max, ind_z_surf, ind_t, detlayers, numOutputDetail, outputDetail
    double precision, intent(in) :: detthick
    double precision, dimension(ind_z_max), intent(in) :: Rho, T, rgrain2, Mlwc, DZ
    double precision, dimension(ind_z_max), intent(inout) :: Refreeze
    double precision, dimension(outputDetail+50,detlayers), intent(out) :: out_2D_det_dens, out_2D_det_temp, out_2D_det_lwc, &
        out_2D_det_refreeze, out_2D_det_rgrain

    ! declare local arguments
    integer :: ind_t_out, ind_orig, ind_int
    double precision :: dist, part, Refreeze_Sum
    double precision, dimension(ind_z_max) :: DZ_mod
    double precision, dimension(detlayers) :: Int_Rho, Int_T, Int_rgrain, Int_Mlwc, Int_Refreeze
        
    Int_Rho(:) = 0.
    Int_T(:) = 0.
    Int_Mlwc(:) = 0.
    Int_Refreeze(:) = 0.
    Int_rgrain(:) = 0.

    DZ_mod = DZ
    dist = 0.
    ind_orig = ind_z_surf
    ind_int = 1
    do while(ind_int <= detlayers)
        if ((dist + DZ_mod(ind_orig)) < (detthick * ind_int)) then
            Int_Rho(ind_int) = Int_Rho(ind_int) + Rho(ind_orig) * (DZ_mod(ind_orig) / detthick)
            Int_T(ind_int) = Int_T(ind_int) + T(ind_orig) * (DZ_mod(ind_orig) / detthick)
            Int_Mlwc(ind_int) = Int_Mlwc(ind_int) + Mlwc(ind_orig) * (DZ_mod(ind_orig) / DZ(ind_orig))
            Int_rgrain(ind_int) = Int_rgrain(ind_int) + (1000.*(rgrain2(ind_orig)**0.5)) * (DZ_mod(ind_orig) / DZ(ind_orig))
            dist = dist + DZ_mod(ind_orig)
            ind_orig = ind_orig - 1
        else if ((dist + DZ_mod(ind_orig)) == (detthick * ind_int)) then
            Int_Rho(ind_int) = Int_Rho(ind_int) + Rho(ind_orig) * (DZ_mod(ind_orig) / detthick)
            Int_T(ind_int) = Int_T(ind_int) + T(ind_orig) * (DZ_mod(ind_orig) / detthick)
            Int_Mlwc(ind_int) = Int_Mlwc(ind_int) + Mlwc(ind_orig) * (DZ_mod(ind_orig) / DZ(ind_orig))
            Int_rgrain(ind_int) = Int_rgrain(ind_int) + (1000.*(rgrain2(ind_orig)**0.5)) * (DZ_mod(ind_orig) / DZ(ind_orig))
            dist = dist + DZ_mod(ind_orig)
            ind_orig = ind_orig - 1
            ind_int = ind_int + 1
        else
            part = (detthick * ind_int) - dist
            Int_Rho(ind_int) = Int_Rho(ind_int) + Rho(ind_orig) * (part / detthick)
            Int_T(ind_int) = Int_T(ind_int) + T(ind_orig) * (part / detthick)
            Int_Mlwc(ind_int) = Int_Mlwc(ind_int) + Mlwc(ind_orig) * (part / DZ(ind_orig))
            Int_rgrain(ind_int) = Int_rgrain(ind_int) + (1000.*(rgrain2(ind_orig)**0.5)) * (part / DZ(ind_orig))
            dist = dist + part
            DZ_mod(ind_orig) = DZ_mod(ind_orig) - part
            ind_int = ind_int + 1
        endif
    end do

    Refreeze_Sum = SUM(Refreeze)
    if (Refreeze_Sum > 1e-05) then
        DZ_mod = DZ
        dist = 0.
        ind_orig = ind_z_surf
        ind_int = 1
        do while(ind_int <= detlayers)
            if ((dist + DZ_mod(ind_orig)) < (detthick * ind_int)) then
                Int_Refreeze(ind_int) = Int_Refreeze(ind_int) + Refreeze(ind_orig) * (DZ_mod(ind_orig) / DZ(ind_orig))
                dist = dist + DZ_mod(ind_orig)
                ind_orig = ind_orig - 1
            else if ((dist + DZ_mod(ind_orig)) == (detthick * ind_int)) then
                Int_Refreeze(ind_int) = Int_Refreeze(ind_int) + Refreeze(ind_orig) * (DZ_mod(ind_orig) / DZ(ind_orig))
                dist = dist + DZ_mod(ind_orig)
                ind_orig = ind_orig - 1
                ind_int = ind_int + 1
            else
                part = (detthick * ind_int) - dist
                Int_Refreeze(ind_int) = Int_Refreeze(ind_int) + Refreeze(ind_orig) * (part / DZ(ind_orig))
                dist = dist + part
                DZ_mod(ind_orig) = DZ_mod(ind_orig) - part
                ind_int = ind_int + 1
            endif
        end do
    endif

    ind_t_out = ind_t/numOutputDetail
    
    out_2D_det_dens(ind_t_out,:) = REAL(Int_Rho)
    out_2D_det_temp(ind_t_out,:) = REAL(Int_T)
    out_2D_det_lwc(ind_t_out,:) = REAL(Int_Mlwc)
    out_2D_det_refreeze(ind_t_out,:) = REAL(Int_Refreeze)
    out_2D_det_rgrain(ind_t_out,:) = REAL(Int_rgrain)

    Refreeze(:) = 0.

end subroutine To_out_2Ddetail


! *******************************************************


subroutine Save_out_1D(outputSpeed, path_out_1d, fname_out_1d, out_1D)
    !*** Write the 1D output variables to a netcdf file !***
    
    ! declare arguments
    integer, intent(in) :: outputSpeed
    double precision, dimension((outputSpeed+50),18), intent(in) :: out_1D
    character*255, intent(in) :: path_out_1d, fname_out_1d

    ! declare local arguments
    integer :: status, ncid, dimid_t, varID(18)
    double precision :: NaN_value

    NaN_value = 9.96921e+36
    
    ! CREATE NETCDF FILES
    status = nf90_create(trim(path_out_1d)//trim(fname_out_1d), 0, ncid)

    ! DEFINE DIMENSIONS
    status = nf90_def_dim(ncid, "ind_t", outputSpeed+50, dimid_t)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_dim1')
    
    ! DEFINE VARIABLES
    status = nf90_def_var(ncid, "h_surf", nf90_real, dimid_t, varID(1))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var1')
    status = nf90_def_var(ncid, "vice", nf90_real, dimid_t, varID(2))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var2')
    status = nf90_def_var(ncid, "vacc", nf90_real, dimid_t, varID(3))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var3')
    status = nf90_def_var(ncid, "vfc", nf90_real, dimid_t, varID(4))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var4')
    status = nf90_def_var(ncid, "vmelt", nf90_real, dimid_t, varID(5))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var5')
    status = nf90_def_var(ncid, "vbouy", nf90_real, dimid_t, varID(6))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var6')
    status = nf90_def_var(ncid, "vsub", nf90_real, dimid_t, varID(7))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var7')
    status = nf90_def_var(ncid, "vsnd", nf90_real, dimid_t, varID(8))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var8')
    status = nf90_def_var(ncid, "vtotal", nf90_real, dimid_t, varID(9))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var9')
    status = nf90_def_var(ncid, "Runoff", nf90_real, dimid_t, varID(10))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var10')
    status = nf90_def_var(ncid, "FirnAir", nf90_real, dimid_t, varID(11))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var11')
    status = nf90_def_var(ncid, "TotLwc", nf90_real, dimid_t, varID(12))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var12')
    status = nf90_def_var(ncid, "refreeze", nf90_real, dimid_t, varID(13)) 
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var13') 
    status = nf90_def_var(ncid, "rain", nf90_real, dimid_t, varID(14))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var14')
    status = nf90_def_var(ncid, "surfmelt", nf90_real, dimid_t, varID(15)) 
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var15') 
    status = nf90_def_var(ncid, "solin", nf90_real, dimid_t, varID(16))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var16') 
    status = nf90_def_var(ncid, "icemass", nf90_real, dimid_t, varID(17))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var17') 
    status = nf90_def_var(ncid, "Rho0", nf90_real, dimid_t, varID(18))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var18')
    
    ! DEFINE ATTRIBUTES (unit could also be defined here)
    status = nf90_put_att(ncid, varID(1), "missing_value", NaN_value)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_att_miss_val_1')
    status = nf90_put_att(ncid, varID(2), "missing_value", NaN_value)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_att_miss_val_2')
    status = nf90_put_att(ncid, varID(3), "missing_value", NaN_value)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_att_miss_val_3')
    status = nf90_put_att(ncid, varID(4), "missing_value", NaN_value)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_att_miss_val_4')
    status = nf90_put_att(ncid, varID(5), "missing_value", NaN_value)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_att_miss_val_5')
    status = nf90_put_att(ncid, varID(6), "missing_value", NaN_value)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_att_miss_val_6')
    status = nf90_put_att(ncid, varID(7), "missing_value", NaN_value)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_att_miss_val_7')
    status = nf90_put_att(ncid, varID(8), "missing_value", NaN_value)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_att_miss_val_8')
    status = nf90_put_att(ncid, varID(9), "missing_value", NaN_value)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_att_miss_val_9')
    status = nf90_put_att(ncid, varID(10), "missing_value", NaN_value)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_att_miss_val_10')
    status = nf90_put_att(ncid, varID(11), "missing_value", NaN_value)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_att_miss_val_11')
    status = nf90_put_att(ncid, varID(12), "missing_value", NaN_value)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_att_miss_val_12')
    status = nf90_put_att(ncid, varID(13), "missing_value", NaN_value)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_att_miss_val_13')
    status = nf90_put_att(ncid, varID(14), "missing_value", NaN_value)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_att_miss_val_14')
    status = nf90_put_att(ncid, varID(15), "missing_value", NaN_value)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_att_miss_val_15')
    status = nf90_put_att(ncid, varID(16), "missing_value", NaN_value)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_att_miss_val_16')
    status = nf90_put_att(ncid, varID(17), "missing_value", NaN_value)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_att_miss_val_17')
    status = nf90_put_att(ncid, varID(18), "missing_value", NaN_value)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_att_miss_val_18')
    
    ! END OF DEFINING FILES
    status = nf90_enddef(ncid)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_enddef')

    ! SAVE DATA
    status = nf90_put_var(ncid, varID(1), out_1D(:,1), start=(/1/), count=(/(outputSpeed+50)/))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var_speed1')    
    status = nf90_put_var(ncid, varID(2), out_1D(:,2), start=(/1/), count=(/(outputSpeed+50)/))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var2')
    status = nf90_put_var(ncid, varID(3), out_1D(:,3), start=(/1/), count=(/(outputSpeed+50)/))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var_speed3')
    status = nf90_put_var(ncid, varID(4), out_1D(:,4), start=(/1/), count=(/(outputSpeed+50)/))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var_speed4')
    status = nf90_put_var(ncid, varID(5), out_1D(:,5), start=(/1/), count=(/(outputSpeed+50)/))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var_speed5')
    status = nf90_put_var(ncid, varID(6), out_1D(:,6), start=(/1/), count=(/(outputSpeed+50)/))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var_speed6')
    status = nf90_put_var(ncid, varID(7), out_1D(:,7), start=(/1/), count=(/(outputSpeed+50)/))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var_speed7')
    status = nf90_put_var(ncid, varID(8), out_1D(:,8), start=(/1/), count=(/(outputSpeed+50)/))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var_speed8')
    status = nf90_put_var(ncid, varID(9), out_1D(:,9), start=(/1/), count=(/(outputSpeed+50)/))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var_speed9')
    status = nf90_put_var(ncid, varID(10), out_1D(:,10), start=(/1/), count=(/(outputSpeed+50)/))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var_speed10')
    status = nf90_put_var(ncid, varID(11), out_1D(:,11), start=(/1/), count=(/(outputSpeed+50)/))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var_speed11')
    status = nf90_put_var(ncid, varID(12), out_1D(:,12), start=(/1/), count=(/(outputSpeed+50)/))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var_speed12')
    status = nf90_put_var(ncid, varID(13), out_1D(:,13), start=(/1/), count=(/(outputSpeed+50)/))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var_speed13')
    status = nf90_put_var(ncid, varID(14), out_1D(:,14), start=(/1/), count=(/(outputSpeed+50)/))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var_speed14')
    status = nf90_put_var(ncid, varID(15), out_1D(:,15), start=(/1/), count=(/(outputSpeed+50)/))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var_speed15')
    status = nf90_put_var(ncid, varID(16), out_1D(:,16), start=(/1/), count=(/(outputSpeed+50)/))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var_speed16')
    status = nf90_put_var(ncid, varID(17), out_1D(:,17), start=(/1/), count=(/(outputSpeed+50)/))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var_speed17')
    status = nf90_put_var(ncid, varID(18), out_1D(:,18), start=(/1/), count=(/(outputSpeed+50)/))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var_speed18')

    ! CLOSE NETCDF-FILE
    status = nf90_close(ncid)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_close')
    
end subroutine Save_out_1D


! *******************************************************


subroutine Save_out_2D(outputProf, proflayers, out_2D_dens, out_2D_temp, out_2D_lwc, out_2D_depth, out_2D_dRho, &
    out_2D_year, out_2D_rgrain, path_out_2d, fname_out_2d)
    !*** Write the 2D output variables to a netcdf file !***

    ! declare arguments
    integer, intent(in) :: outputProf, proflayers
    double precision, dimension((outputProf+50),proflayers), intent(in) :: out_2D_dens, out_2D_temp, out_2D_lwc, &
        out_2D_depth, out_2D_dRho, out_2D_year, out_2D_rgrain
    character*255, intent(in) :: path_out_2d, fname_out_2d

    ! declare local arguments
    integer :: status, ncid, dimid_t, dimid_z, varID(7)
    double precision :: NaN_value
    
    NaN_value = 9.96921e+36

    ! CREATE NETCDF FILES
    status = nf90_create(trim(path_out_2d)//trim(fname_out_2d), 0, ncid)

    ! DEFINE DIMENSIONS
    status = nf90_def_dim(ncid, "ind_t", outputProf+50, dimid_t)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_dim2')    
    status = nf90_def_dim(ncid, "layer", proflayers, dimid_z)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_dim3')

    ! DEFINE VARIABLES
    status = nf90_def_var(ncid, "dens", nf90_real, (/dimid_t, dimid_z/), varID(1))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var1')
    status = nf90_def_var(ncid, "temp", nf90_real, (/dimid_t, dimid_z/), varID(2))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var2')
    status = nf90_def_var(ncid, "year", nf90_real, (/dimid_t, dimid_z/), varID(3))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var3')
    status = nf90_def_var(ncid, "lwc", nf90_real, (/dimid_t, dimid_z/), varID(4))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var4')
    status = nf90_def_var(ncid, "depth", nf90_real, (/dimid_t, dimid_z/), varID(5))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var5')
    status = nf90_def_var(ncid, "dRho", nf90_real, (/dimid_t, dimid_z/), varID(6))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var6')
    status = nf90_def_var(ncid, "rgrain", nf90_real, (/dimid_t, dimid_z/), varID(7))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var6')

    ! DEFINE ATTRIBUTES (unit could also be defined here)
    status = nf90_put_att(ncid, varID(1), "missing_value", NaN_value)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_att_miss_val_1')
    status = nf90_put_att(ncid, varID(2), "missing_value", NaN_value)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_att_miss_val_2')
    status = nf90_put_att(ncid, varID(3), "missing_value", NaN_value)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_att_miss_val_3')
    status = nf90_put_att(ncid, varID(4), "missing_value", NaN_value)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_att_miss_val_4')
    status = nf90_put_att(ncid, varID(5), "missing_value", NaN_value)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_att_miss_val_5')
    status = nf90_put_att(ncid, varID(6), "missing_value", NaN_value)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_att_miss_val_6')
    status = nf90_put_att(ncid, varID(7), "missing_value", NaN_value)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_att_miss_val_7')

    ! END OF DEFINING FILES
    status = nf90_enddef(ncid)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_enddef')

    ! SAVE DATA
    status = nf90_put_var(ncid, varID(1), out_2D_dens, start=(/1,1/), count=(/(outputProf+50), proflayers/))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var_grid1')    
    status = nf90_put_var(ncid, varID(2), out_2D_temp, start=(/1,1/), count=(/(outputProf+50), proflayers/))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var_grid2')
    status = nf90_put_var(ncid, varID(3), out_2D_year, start=(/1,1/), count=(/(outputProf+50), proflayers/))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var_grid3')
    status = nf90_put_var(ncid, varID(4), out_2D_lwc, start=(/1,1/), count=(/(outputProf+50), proflayers/))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var_grid4')
    status = nf90_put_var(ncid, varID(5), out_2D_depth, start=(/1,1/), count=(/(outputProf+50), proflayers/))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var_grid5')
    status = nf90_put_var(ncid, varID(6), out_2D_dRho, start=(/1,1/), count=(/(outputProf+50), proflayers/))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var_grid6')
    status = nf90_put_var(ncid, varID(7), out_2D_rgrain**0.5, start=(/1,1/), count=(/(outputProf+50), proflayers/))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var_grid7')
    
    ! CLOSE NETCDF-FILE
    status = nf90_close(ncid)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_close')
    
end subroutine Save_out_2D


! *******************************************************


subroutine Save_out_2Ddetail(outputDetail, detlayers, detthick, out_2D_det_dens, out_2D_det_temp, out_2D_det_lwc, &
    out_2D_det_refreeze, out_2D_det_rgrain, path_out_2ddet, fname_out_2ddet)
    !*** Write the 2Ddetail output variables to a netcdf file !***
    
    ! declare arguments
    integer, intent(in) :: outputDetail, detlayers
    double precision, intent(in) :: detthick
    double precision, dimension((outputDetail+50),detlayers), intent(in) :: out_2D_det_dens, out_2D_det_temp, &
        out_2D_det_lwc, out_2D_det_refreeze, out_2D_det_rgrain
    character*255, intent(in) :: path_out_2ddet, fname_out_2ddet

    ! declare local arguments
    integer :: ind_z_det, status, ncid, dimid_t, dimid_z, varID(7)
    double precision :: NaN_value
    double precision, dimension(detlayers) :: DetDepth, DetDZ
    
    NaN_value = 9.96921e+36
    
    ! Calculate the thickness of depth at the interpolated depths
    DetDZ(:) = detthick
    DetDepth(1) = DetDZ(1) / 2.
    do ind_z_det = 2, detlayers
        DetDepth(ind_z_det) = DetDepth(ind_z_det-1) + DetDZ(ind_z_det)
    end do

    ! CREATE NETCDF FILES
    status = nf90_create(trim(path_out_2ddet)//trim(fname_out_2ddet), 0, ncid)    

    ! DEFINE DIMENSIONS
    status = nf90_def_dim(ncid, "ind_t", (outputDetail+50), dimid_t)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_dim4')    
    status = nf90_def_dim(ncid, "layer", detlayers, dimid_z)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_dim5')

    ! DEFINE VARIABLES
    status = nf90_def_var(ncid, "dens", nf90_real, (/dimid_t, dimid_z/), varID(1))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var1')
    status = nf90_def_var(ncid, "temp", nf90_real, (/dimid_t, dimid_z/), varID(2))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var2')
    status = nf90_def_var(ncid, "lwc", nf90_real, (/dimid_t, dimid_z/), varID(3))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var3')
    status = nf90_def_var(ncid, "depth", nf90_real, (/dimid_z/), varID(4))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var4')
    status = nf90_def_var(ncid, "dz", nf90_real, (/dimid_z/), varID(5))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var5')
    status = nf90_def_var(ncid, "refreeze", nf90_real, (/dimid_t, dimid_z/), varID(6))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var6')
    status = nf90_def_var(ncid, "rgrain", nf90_real, (/dimid_t, dimid_z/), varID(7))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var7')
    
    ! DEFINE ATTRIBUTES (unit could also be defined here)
    status = nf90_put_att(ncid, varID(1), "missing_value", NaN_value)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_att_miss_val_1')
    status = nf90_put_att(ncid, varID(2), "missing_value", NaN_value)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_att_miss_val_2')
    status = nf90_put_att(ncid, varID(3), "missing_value", NaN_value)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_att_miss_val_3')
    status = nf90_put_att(ncid, varID(6), "missing_value", NaN_value)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_att_miss_val_6')
    status = nf90_put_att(ncid, varID(7), "missing_value", NaN_value)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_att_miss_val_7')
    
    ! END OF DEFINING FILES
    status = nf90_enddef(ncid)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_enddef')
    
    ! PUT IN THE CONSTANT Z-AXIS in 2Ddetail
    status = nf90_put_var(ncid, varID(4),DetDepth)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var1_1')
    status = nf90_put_var(ncid, varID(5),DetDZ)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var2_1')
    
    ! SAVE DATA
    status = nf90_put_var(ncid, varID(1), out_2D_det_dens, start=(/1,1/), count=(/(outputDetail+50), detlayers/))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var_detail1')    
    status = nf90_put_var(ncid, varID(2), out_2D_det_temp, start=(/1,1/), count=(/(outputDetail+50), detlayers/))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var_detail2')
    status = nf90_put_var(ncid, varID(3), out_2D_det_lwc, start=(/1,1/), count=(/(outputDetail+50), detlayers/))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var_detail3')
    status = nf90_put_var(ncid, varID(6), out_2D_det_refreeze, start=(/1,1/), count=(/(outputDetail+50), detlayers/))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var_detail6')
    status = nf90_put_var(ncid, varID(7), out_2D_det_rgrain, start=(/1,1/), count=(/(outputDetail+50), detlayers/))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var_detail7')

    ! CLOSE NETCDF-FILE
    status = nf90_close(ncid)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_close')
    
end subroutine Save_out_2Ddetail


! *******************************************************


subroutine Save_out_restart(ind_t, ind_z_max, ind_z_surf, Rho, DenRho, M, T, Depth, Mlwc, Year, Refreeze, DZ, path_out_restart, fname_out_restart)
    !*** Write a netcdf file containing the firn profile and all the variables necessary for the restart functionality ***!
    
    ! declare arguments
    integer, intent(in) :: ind_z_max, ind_z_surf, ind_t
    double precision, dimension(ind_z_max), intent(in) ::  M, T, Depth, Mlwc, DenRho, Rho, Year, Refreeze, DZ
    character*255, intent(in) :: path_out_restart, fname_out_restart

    ! declare local arguments
    integer :: status, ncid, dimid_t(2), varID(11)
        
    ! CREATE NETCDF FILES
    status = nf90_create(trim(path_out_restart)//trim(fname_out_restart), 0, ncid)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_create')
    
    ! DEFINE DIMENSIONS
    status = nf90_def_dim(ncid, "layer", ind_z_max, dimid_t(1))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_dim1')
    status = nf90_def_dim(ncid, "constant", 1, dimid_t(2))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_dim2')
    
    ! DEFINE VARIABLES
    status = nf90_def_var(ncid, "dens", nf90_real, dimid_t(1), varID(1))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var1')
    status = nf90_def_var(ncid, "temp", nf90_real, dimid_t(1), varID(2))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var2')
    status = nf90_def_var(ncid, "mass", nf90_real, dimid_t(1), varID(3))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var3')
    status = nf90_def_var(ncid, "depth", nf90_real, dimid_t(1), varID(4))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var4')
    status = nf90_def_var(ncid, "lwc", nf90_real, dimid_t(1), varID(5))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var5')
    status = nf90_def_var(ncid, "year", nf90_real, dimid_t(1), varID(6))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var6')
    status = nf90_def_var(ncid, "refreeze", nf90_real, dimid_t(1), varID(7))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var7')
    status = nf90_def_var(ncid, "dz", nf90_real, dimid_t(1), varID(8))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var8')
    status = nf90_def_var(ncid, "drho", nf90_real, dimid_t(1), varID(9))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var9')
    
    status = nf90_def_var(ncid, "ind_z_surf", nf90_int, dimid_t(2), varID(10))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var10')
    status = nf90_def_var(ncid, "ind_t", nf90_int, dimid_t(2), varID(11))
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_def_var11')
    
    ! END OF DEFINING FILES
    status = nf90_enddef(ncid)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_enddef')
    
    ! PUT VARIABLES
    status = nf90_put_var(ncid, varID(1), Rho)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var1')
    status = nf90_put_var(ncid, varID(2), T)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var2')
    status = nf90_put_var(ncid, varID(3), M)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var3')
    status = nf90_put_var(ncid, varID(4), Depth)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var4')
    status = nf90_put_var(ncid, varID(5), Mlwc)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var5')
    status = nf90_put_var(ncid, varID(6), Year)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var6')
    status = nf90_put_var(ncid, varID(7), Refreeze)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var7')
    status = nf90_put_var(ncid, varID(8), DZ)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var8')
    status = nf90_put_var(ncid, varID(9), DenRho)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var9')
    
    status = nf90_put_var(ncid, varID(10), ind_z_surf)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var10')
    status = nf90_put_var(ncid, varID(11), ind_t)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_put_var11')
    
    ! CLOSE NETCDF FILE
    status = nf90_close(ncid)
    if(status /= nf90_noerr) call Handle_Error(status, 'nf_close')
    
end subroutine Save_out_restart


end module output
