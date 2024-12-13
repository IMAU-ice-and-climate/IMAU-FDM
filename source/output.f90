module output
    !*** functions and subroutines for writing data to netcdf files

    use netcdf, only: nf90_create, nf90_def_dim, nf90_def_var, nf90_real, nf90_int, nf90_noerr, &
        nf90_enddef, nf90_put_var, nf90_close, nf90_unlimited, nf90_put_att
    use openNetCDF, only: Handle_Error

    implicit none
    private

    public :: Accumulate_Output, To_out_1D, To_out_2D, To_out_2Ddetail, Save_out_1D, Save_out_2D, Save_out_2Ddetail, Save_out_spinup, Save_out_run
    
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
        T, Mlwc, Depth, DenRho, Year, out_2D_dens, out_2D_temp, out_2D_lwc, out_2D_depth, out_2D_dRho, &
        out_2D_year)
    !*** Write the 2D output variables to the variables that will be converted into a netcdf file after the time loop ***!
    
    ! declare arguments
    integer, intent(in) :: ind_z_max, ind_z_surf, ind_t, dtmodel, numOutputProf, outputProf, proflayers
    double precision, dimension(ind_z_max), intent(in) :: Rho, T, Mlwc, Depth, Year
    double precision, dimension(ind_z_max), intent(inout) :: DenRho
    double precision, dimension(outputProf+50,proflayers), intent(out) :: out_2D_dens, out_2D_temp, out_2D_lwc, &
        out_2D_depth, out_2D_dRho, out_2D_year

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
    
    DenRho(:) = 0.
    
end subroutine To_out_2D


! *******************************************************


subroutine To_out_2Ddetail(ind_z_max, ind_z_surf, ind_t,detlayers, detthick, numOutputDetail, outputDetail, &
    Rho, T, Mlwc, Refreeze, DZ, out_2D_det_dens, out_2D_det_temp, out_2D_det_lwc, out_2D_det_refreeze)
    !*** Write the 2Ddetail output variables to the variables that will be converted into a netcdf file after the time loop ***!
    
    ! declare arguments
    integer, intent(in) :: ind_z_max, ind_z_surf, ind_t, detlayers, numOutputDetail, outputDetail
    double precision, intent(in) :: detthick
    double precision, dimension(ind_z_max), intent(in) :: Rho, T, Mlwc, DZ
    double precision, dimension(ind_z_max), intent(inout) :: Refreeze
    double precision, dimension(outputDetail+50,detlayers), intent(out) :: out_2D_det_dens, out_2D_det_temp, out_2D_det_lwc, out_2D_det_refreeze

    ! declare local arguments
    integer :: ind_t_out, ind_orig, ind_int
    double precision :: dist, part, refreeze_sum
    double precision, dimension(ind_z_max) :: DZ_mod
    double precision, dimension(detlayers) :: IntRho,IntT,IntMlwc,IntRefreeze
        
    IntRho(:) = 0.
    IntT(:) = 0.
    IntMlwc(:) = 0.
    IntRefreeze(:) = 0.

    print *, "is this a memory problem?"
    print *, " detthick = ", detthick

    DZ_mod = DZ
    dist = 0.
    ind_orig = ind_z_surf
    ind_int = 1
    do while(ind_int <= detlayers)
        if ((dist + DZ_mod(ind_orig)) < (detthick * ind_int)) then
            IntRho(ind_int) = IntRho(ind_int) + Rho(ind_orig) * (DZ_mod(ind_orig) / detthick)
            IntT(ind_int) = IntT(ind_int) + T(ind_orig) * (DZ_mod(ind_orig) / detthick)
            IntMlwc(ind_int) = IntMlwc(ind_int) + Mlwc(ind_orig) * (DZ_mod(ind_orig) / DZ(ind_orig))
            dist = dist + DZ_mod(ind_orig)
            ind_orig = ind_orig - 1
        else if ((dist + DZ_mod(ind_orig)) == (detthick * ind_int)) then
            IntRho(ind_int) = IntRho(ind_int) + Rho(ind_orig) * (DZ_mod(ind_orig) / detthick)
            IntT(ind_int) = IntT(ind_int) + T(ind_orig) * (DZ_mod(ind_orig) / detthick)
            IntMlwc(ind_int) = IntMlwc(ind_int) + Mlwc(ind_orig) * (DZ_mod(ind_orig) / DZ(ind_orig))
            dist = dist + DZ_mod(ind_orig)
            ind_orig = ind_orig - 1
            ind_int = ind_int + 1
        else
            part = (detthick * ind_int) - dist
            IntRho(ind_int) = IntRho(ind_int) + Rho(ind_orig) * (part / detthick)
            IntT(ind_int) = IntT(ind_int) + T(ind_orig) * (part / detthick)
            IntMlwc(ind_int) = IntMlwc(ind_int) + Mlwc(ind_orig) * (part / DZ(ind_orig))        
            dist = dist + part
            DZ_mod(ind_orig) = DZ_mod(ind_orig) - part
            ind_int = ind_int + 1
        end if
    end do

    refreeze_sum = SUM(Refreeze)
    if (refreeze_sum > 1e-05) then
        DZ_mod = DZ
        dist = 0.
        ind_orig = ind_z_surf
        ind_int = 1
        do while(ind_int <= detlayers)
            if ((dist + DZ_mod(ind_orig)) < (detthick * ind_int)) then
                IntRefreeze(ind_int) = Intrefreeze(ind_int) + Refreeze(ind_orig) * (DZ_mod(ind_orig) / DZ(ind_orig))
                dist = dist + DZ_mod(ind_orig)
                ind_orig = ind_orig - 1
            else if ((dist + DZ_mod(ind_orig)) == (detthick * ind_int)) then
                IntRefreeze(ind_int) = Intrefreeze(ind_int) + Refreeze(ind_orig) * (DZ_mod(ind_orig) / DZ(ind_orig))
                dist = dist + DZ_mod(ind_orig)
                ind_orig = ind_orig - 1
                ind_int = ind_int + 1
            else
                part = (detthick * ind_int) - dist
                IntRefreeze(ind_int) = Intrefreeze(ind_int) + Refreeze(ind_orig) * (part / DZ(ind_orig))
                dist = dist + part
                DZ_mod(ind_orig) = DZ_mod(ind_orig) - part
                ind_int = ind_int + 1
            end if
        end do
    end if

    ind_t_out = ind_t/numOutputDetail
    
    out_2D_det_dens(ind_t_out,:) = REAL(IntRho)
    out_2D_det_temp(ind_t_out,:) = REAL(IntT)
    out_2D_det_lwc(ind_t_out,:) = REAL(IntMlwc)
    out_2D_det_refreeze(ind_t_out,:) = REAL(IntRefreeze)

    Refreeze(:) = 0.

end subroutine To_out_2Ddetail


! *******************************************************


subroutine Save_out_1D(outputSpeed, point_numb, fname_p1, username, out_1D, project_name)
    !*** Write the 1D output variables to a netcdf file !***
    
    ! declare arguments
    integer, intent(in) :: outputSpeed
    double precision, dimension((outputSpeed+50),18), intent(in) :: out_1D
    character*255, intent(in) :: point_numb, fname_p1, username, project_name

    ! declare local arguments
    integer :: status, ncid(50), IDs(50,5), varID(50,20)
    character*255 :: pad

    pad = "/ec/res4/scratch/"//trim(username)//"/"//trim(project_name)//"/output/"
    
    ncid(:) = 0
    IDs(:,:) = 0
    varID(:,:) = 0

    ! CREATE NETCDF FILES
    status = nf90_create(trim(pad)//trim(fname_p1)//"_1D_"//trim(point_numb)//".nc",0,ncid(32))

    ! DEFINE DIMENSIONS
    status = nf90_def_dim(ncid(32),"ind_t",outputSpeed+50,IDs(32,1))
    if(status /= nf90_noerr) call Handle_Error(status,'1D_def_dim1')
    
    ! DEFINE VARIABLES
    status = nf90_def_var(ncid(32),"h_surf",nf90_real,(/IDs(32,1)/),varID(32,1))
    if(status /= nf90_noerr) call Handle_Error(status,'1D_def_var1')
    status = nf90_def_var(ncid(32),"vice",nf90_real,(/IDs(32,1)/),varID(32,2))
    if(status /= nf90_noerr) call Handle_Error(status,'1D_def_var2')
    status = nf90_def_var(ncid(32),"vacc",nf90_real,(/IDs(32,1)/),varID(32,3))
    if(status /= nf90_noerr) call Handle_Error(status,'1D_def_var3')
    status = nf90_def_var(ncid(32),"vfc",nf90_real,(/IDs(32,1)/),varID(32,4))
    if(status /= nf90_noerr) call Handle_Error(status,'1D_def_var4')
    status = nf90_def_var(ncid(32),"vmelt",nf90_real,(/IDs(32,1)/),varID(32,5))
    if(status /= nf90_noerr) call Handle_Error(status,'1D_def_var5')
    status = nf90_def_var(ncid(32),"vbouy",nf90_real,(/IDs(32,1)/),varID(32,6))
    if(status /= nf90_noerr) call Handle_Error(status,'1D_def_var6')
    status = nf90_def_var(ncid(32),"vsub",nf90_real,(/IDs(32,1)/),varID(32,7))
    if(status /= nf90_noerr) call Handle_Error(status,'1D_def_var7')
    status = nf90_def_var(ncid(32),"vsnd",nf90_real,(/IDs(32,1)/),varID(32,8))
    if(status /= nf90_noerr) call Handle_Error(status,'1D_def_var8')
    status = nf90_def_var(ncid(32),"vtotal",nf90_real,(/IDs(32,1)/),varID(32,9))
    if(status /= nf90_noerr) call Handle_Error(status,'1D_def_var9')
    status = nf90_def_var(ncid(32),"Runoff",nf90_real,(/IDs(32,1)/),varID(32,10))
    if(status /= nf90_noerr) call Handle_Error(status,'1D_def_var10')
    status = nf90_def_var(ncid(32),"FirnAir",nf90_real,(/IDs(32,1)/),varID(32,11))
    if(status /= nf90_noerr) call Handle_Error(status,'1D_def_var11')
    status = nf90_def_var(ncid(32),"TotLwc",nf90_real,(/IDs(32,1)/),varID(32,12))
    if(status /= nf90_noerr) call Handle_Error(status,'1D_def_var12')
    status = nf90_def_var(ncid(32),"refreeze",nf90_real,(/IDs(32,1)/),varID(32,13)) 
    if(status /= nf90_noerr) call Handle_Error(status,'1D_def_var13') 
    status = nf90_def_var(ncid(32),"rain",nf90_real,(/IDs(32,1)/),varID(32,14))
    if(status /= nf90_noerr) call Handle_Error(status,'1D_def_var14')
    status = nf90_def_var(ncid(32),"surfmelt",nf90_real,(/IDs(32,1)/),varID(32,15)) 
    if(status /= nf90_noerr) call Handle_Error(status,'1D_def_var15') 
    status = nf90_def_var(ncid(32),"solin",nf90_real,(/IDs(32,1)/),varID(32,16))
    if(status /= nf90_noerr) call Handle_Error(status,'1D_def_var16') 
    status = nf90_def_var(ncid(32),"icemass",nf90_real,(/IDs(32,1)/),varID(32,17))
    if(status /= nf90_noerr) call Handle_Error(status,'1D_def_var17') 
    status = nf90_def_var(ncid(32),"Rho0",nf90_real,(/IDs(32,1)/),varID(32,18))
    if(status /= nf90_noerr) call Handle_Error(status,'1D_def_var18')
    
    ! DEFINE ATTRIBUTES (unit could also be defined here)
    status = nf90_put_att(ncid(32),varID(32,1),"missing_value",9.96921e+36)
    if(status /= nf90_noerr) call Handle_Error(status,'1D_def_att_miss_val_1')
    status = nf90_put_att(ncid(32),varID(32,2),"missing_value",9.96921e+36)
    if(status /= nf90_noerr) call Handle_Error(status,'1D_def_att_miss_val_2')
    status = nf90_put_att(ncid(32),varID(32,3),"missing_value",9.96921e+36)
    if(status /= nf90_noerr) call Handle_Error(status,'1D_def_att_miss_val_3')
    status = nf90_put_att(ncid(32),varID(32,4),"missing_value",9.96921e+36)
    if(status /= nf90_noerr) call Handle_Error(status,'1D_def_att_miss_val_4')
    status = nf90_put_att(ncid(32),varID(32,5),"missing_value",9.96921e+36)
    if(status /= nf90_noerr) call Handle_Error(status,'1D_def_att_miss_val_5')
    status = nf90_put_att(ncid(32),varID(32,6),"missing_value",9.96921e+36)
    if(status /= nf90_noerr) call Handle_Error(status,'1D_def_att_miss_val_6')
    status = nf90_put_att(ncid(32),varID(32,7),"missing_value",9.96921e+36)
    if(status /= nf90_noerr) call Handle_Error(status,'1D_def_att_miss_val_7')
    status = nf90_put_att(ncid(32),varID(32,8),"missing_value",9.96921e+36)
    if(status /= nf90_noerr) call Handle_Error(status,'1D_def_att_miss_val_8')
    status = nf90_put_att(ncid(32),varID(32,9),"missing_value",9.96921e+36)
    if(status /= nf90_noerr) call Handle_Error(status,'1D_def_att_miss_val_9')
    status = nf90_put_att(ncid(32),varID(32,10),"missing_value",9.96921e+36)
    if(status /= nf90_noerr) call Handle_Error(status,'1D_def_att_miss_val_10')
    status = nf90_put_att(ncid(32),varID(32,11),"missing_value",9.96921e+36)
    if(status /= nf90_noerr) call Handle_Error(status,'1D_def_att_miss_val_11')
    status = nf90_put_att(ncid(32),varID(32,12),"missing_value",9.96921e+36)
    if(status /= nf90_noerr) call Handle_Error(status,'1D_def_att_miss_val_12')
    status = nf90_put_att(ncid(32),varID(32,13),"missing_value",9.96921e+36)
    if(status /= nf90_noerr) call Handle_Error(status,'1D_def_att_miss_val_13')
    status = nf90_put_att(ncid(32),varID(32,14),"missing_value",9.96921e+36)
    if(status /= nf90_noerr) call Handle_Error(status,'1D_def_att_miss_val_14')
    status = nf90_put_att(ncid(32),varID(32,15),"missing_value",9.96921e+36)
    if(status /= nf90_noerr) call Handle_Error(status,'1D_def_att_miss_val_15')
    status = nf90_put_att(ncid(32),varID(32,16),"missing_value",9.96921e+36)
    if(status /= nf90_noerr) call Handle_Error(status,'1D_def_att_miss_val_16')
    status = nf90_put_att(ncid(32),varID(32,17),"missing_value",9.96921e+36)
    if(status /= nf90_noerr) call Handle_Error(status,'1D_def_att_miss_val_17')    
    status = nf90_put_att(ncid(32),varID(32,18),"missing_value",9.96921e+36)
    if(status /= nf90_noerr) call Handle_Error(status,'1D_def_att_miss_val_18')
    
    ! END OF DEFINING FILES
    status = nf90_enddef(ncid(32))
    if(status /= nf90_noerr) call Handle_Error(status,'1D_enddef')

    ! SAVE DATA
    status = nf90_put_var(ncid(32),varID(32,1),out_1D(:,1), &
    start=(/1/),count=(/(outputSpeed+50)/))
    if(status /= nf90_noerr) call Handle_Error(status,'1D_put_var_speed1')    
    status = nf90_put_var(ncid(32),varID(32,2),out_1D(:,2), &
    start=(/1/),count=(/(outputSpeed+50)/))
    if(status /= nf90_noerr) call Handle_Error(status,'1D_put_var2')
    status = nf90_put_var(ncid(32),varID(32,3),out_1D(:,3), &
    start=(/1/),count=(/(outputSpeed+50)/))
    if(status /= nf90_noerr) call Handle_Error(status,'1D_put_var_speed3')
    status = nf90_put_var(ncid(32),varID(32,4),out_1D(:,4), &
    start=(/1/),count=(/(outputSpeed+50)/))
    if(status /= nf90_noerr) call Handle_Error(status,'1D_put_var_speed4')
    status = nf90_put_var(ncid(32),varID(32,5),out_1D(:,5), &
    start=(/1/),count=(/(outputSpeed+50)/))
    if(status /= nf90_noerr) call Handle_Error(status,'1D_put_var_speed5')
    status = nf90_put_var(ncid(32),varID(32,6),out_1D(:,6), &
    start=(/1/),count=(/(outputSpeed+50)/))
    if(status /= nf90_noerr) call Handle_Error(status,'1D_put_var_speed6')
    status = nf90_put_var(ncid(32),varID(32,7),out_1D(:,7), &
    start=(/1/),count=(/(outputSpeed+50)/))
    if(status /= nf90_noerr) call Handle_Error(status,'1D_put_var_speed7')
    status = nf90_put_var(ncid(32),varID(32,8),out_1D(:,8), &
    start=(/1/),count=(/(outputSpeed+50)/))
    if(status /= nf90_noerr) call Handle_Error(status,'1D_put_var_speed8')
    status = nf90_put_var(ncid(32),varID(32,9),out_1D(:,9), &
    start=(/1/),count=(/(outputSpeed+50)/))
    if(status /= nf90_noerr) call Handle_Error(status,'1D_put_var_speed9')
    status = nf90_put_var(ncid(32),varID(32,10),out_1D(:,10), &
    start=(/1/),count=(/(outputSpeed+50)/))
    if(status /= nf90_noerr) call Handle_Error(status,'1D_put_var_speed10')
    status = nf90_put_var(ncid(32),varID(32,11),out_1D(:,11), &
    start=(/1/),count=(/(outputSpeed+50)/))
    if(status /= nf90_noerr) call Handle_Error(status,'1D_put_var_speed11')
    status = nf90_put_var(ncid(32),varID(32,12),out_1D(:,12), &
    start=(/1/),count=(/(outputSpeed+50)/))
    if(status /= nf90_noerr) call Handle_Error(status,'1D_put_var_speed12')
    status = nf90_put_var(ncid(32),varID(32,13),out_1D(:,13), &
    start=(/1/),count=(/(outputSpeed+50)/))
    if(status /= nf90_noerr) call Handle_Error(status,'1D_put_var_speed13')
    status = nf90_put_var(ncid(32),varID(32,14),out_1D(:,14), &
    start=(/1/),count=(/(outputSpeed+50)/))
    if(status /= nf90_noerr) call Handle_Error(status,'1D_put_var_speed14')
    status = nf90_put_var(ncid(32),varID(32,15),out_1D(:,15), &
    start=(/1/),count=(/(outputSpeed+50)/))
    if(status /= nf90_noerr) call Handle_Error(status,'1D_put_var_speed15')
    status = nf90_put_var(ncid(32),varID(32,16),out_1D(:,16), &
    start=(/1/),count=(/(outputSpeed+50)/))
    if(status /= nf90_noerr) call Handle_Error(status,'1D_put_var_speed16')
    status = nf90_put_var(ncid(32),varID(32,17),out_1D(:,17), &
    start=(/1/),count=(/(outputSpeed+50)/))
    if(status /= nf90_noerr) call Handle_Error(status,'1D_put_var_speed17')
    status = nf90_put_var(ncid(32),varID(32,18),out_1D(:,18), &
    start=(/1/),count=(/(outputSpeed+50)/))
    if(status /= nf90_noerr) call Handle_Error(status,'1D_put_var_speed18')

    
    ! CLOSE NETCDF-FILE
    status = nf90_close(ncid(32))
    if(status /= nf90_noerr) call Handle_Error(status,'1D_close')
    
end subroutine Save_out_1D


! *******************************************************


subroutine Save_out_2D(outputProf, proflayers, out_2D_dens, out_2D_temp, out_2D_lwc, out_2D_depth, out_2D_dRho, &
    out_2D_year, point_numb, fname_p1, username, project_name)
    !*** Write the 2D output variables to a netcdf file !***

    ! declare arguments
    integer, intent(in) :: outputProf, proflayers
    double precision, dimension((outputProf+50),proflayers), intent(in) :: out_2D_dens, out_2D_temp, out_2D_lwc, &
        out_2D_depth, out_2D_dRho, out_2D_year
    character*255, intent(in) :: point_numb, fname_p1, username, project_name

    ! declare local arguments
    integer :: status, ncid(50), IDs(50,5), varID(50,20)
    character*255 :: pad

    pad = "/ec/res4/scratch/"//trim(username)//"/"//trim(project_name)//"/output/"
    
    ncid(:) = 0
    IDs(:,:) = 0
    varID(:,:) = 0

    ! CREATE NETCDF FILES
    status = nf90_create(trim(pad)//trim(fname_p1)//"_2D_"//trim(point_numb)//".nc",0,ncid(31))

    ! DEFINE DIMENSIONS
    status = nf90_def_dim(ncid(31),"ind_t",outputProf+50,IDs(31,1))
    if(status /= nf90_noerr) call Handle_Error(status,'2D_def_dim2')    
    status = nf90_def_dim(ncid(31),"layer",proflayers,IDs(31,2))
    if(status /= nf90_noerr) call Handle_Error(status,'2D_def_dim3')

    ! DEFINE VARIABLES
    status = nf90_def_var(ncid(31),"dens",nf90_real,(/IDs(31,1),IDs(31,2)/), &
        varID(31,1))
    if(status /= nf90_noerr) call Handle_Error(status,'2D_def_var1')
    status = nf90_def_var(ncid(31),"temp",nf90_real,(/IDs(31,1),IDs(31,2)/), &
        varID(31,2))
    if(status /= nf90_noerr) call Handle_Error(status,'2D_def_var2')
    status = nf90_def_var(ncid(31),"year",nf90_real,(/IDs(31,1),IDs(31,2)/), &
        varID(31,3))
    if(status /= nf90_noerr) call Handle_Error(status,'2D_def_var3')
    status = nf90_def_var(ncid(31),"lwc",nf90_real,(/IDs(31,1),IDs(31,2)/), &
        varID(31,4))
    if(status /= nf90_noerr) call Handle_Error(status,'2D_def_var4')
    status = nf90_def_var(ncid(31),"depth",nf90_real,(/IDs(31,1),IDs(31,2)/), &
        varID(31,5))
    if(status /= nf90_noerr) call Handle_Error(status,'2D_def_var5')
    status = nf90_def_var(ncid(31),"dRho",nf90_real,(/IDs(31,1),IDs(31,2)/), &
        varID(31,6))
    if(status /= nf90_noerr) call Handle_Error(status,'2D_def_var6')

    ! DEFINE ATTRIBUTES (unit could also be defined here)
    status = nf90_put_att(ncid(31),varID(31,1),"missing_value",9.96921e+36)
    if(status /= nf90_noerr) call Handle_Error(status,'2D_def_att_miss_val_1')
    status = nf90_put_att(ncid(31),varID(31,2),"missing_value",9.96921e+36)
    if(status /= nf90_noerr) call Handle_Error(status,'2D_def_att_miss_val_2')
    status = nf90_put_att(ncid(31),varID(31,3),"missing_value",9.96921e+36)
    if(status /= nf90_noerr) call Handle_Error(status,'2D_def_att_miss_val_3')
    status = nf90_put_att(ncid(31),varID(31,4),"missing_value",9.96921e+36)
    if(status /= nf90_noerr) call Handle_Error(status,'2D_def_att_miss_val_4')
    status = nf90_put_att(ncid(31),varID(31,5),"missing_value",9.96921e+36)
    if(status /= nf90_noerr) call Handle_Error(status,'2D_def_att_miss_val_5')
    status = nf90_put_att(ncid(31),varID(31,6),"missing_value",9.96921e+36)
    if(status /= nf90_noerr) call Handle_Error(status,'2D_def_att_miss_val_6')

    ! END OF DEFINING FILES
    status = nf90_enddef(ncid(31))
    if(status /= nf90_noerr) call Handle_Error(status,'2D_enddef')

    ! SAVE DATA
    status = nf90_put_var(ncid(31),varID(31,1),out_2D_dens,start=(/1,1/), &
        count=(/(outputProf+50),proflayers/))
    if(status /= nf90_noerr) call Handle_Error(status,'2D_put_var_grid1')    
    status = nf90_put_var(ncid(31),varID(31,2),out_2D_temp,start=(/1,1/), &
        count=(/(outputProf+50),proflayers/))
    if(status /= nf90_noerr) call Handle_Error(status,'2D_put_var_grid2')
    status = nf90_put_var(ncid(31),varID(31,3),out_2D_year,start=(/1,1/), &
        count=(/(outputProf+50),proflayers/))
    if(status /= nf90_noerr) call Handle_Error(status,'2D_put_var_grid3')
    status = nf90_put_var(ncid(31),varID(31,4),out_2D_lwc,start=(/1,1/), &
        count=(/(outputProf+50),proflayers/))
    if(status /= nf90_noerr) call Handle_Error(status,'2D_put_var_grid4')
    status = nf90_put_var(ncid(31),varID(31,5),out_2D_depth,start=(/1,1/), &
        count=(/(outputProf+50),proflayers/))
    if(status /= nf90_noerr) call Handle_Error(status,'2D_put_var_grid5')
    status = nf90_put_var(ncid(31),varID(31,6),out_2D_dRho,start=(/1,1/), &
        count=(/(outputProf+50),proflayers/))
    if(status /= nf90_noerr) call Handle_Error(status,'2D_put_var_grid6')
    
    ! CLOSE NETCDF-FILE
    status = nf90_close(ncid(31))
    if(status /= nf90_noerr) call Handle_Error(status,'2D_close')
    
end subroutine Save_out_2D


! *******************************************************


subroutine Save_out_2Ddetail(outputDetail, detlayers, detthick, out_2D_det_dens, out_2D_det_temp, out_2D_det_lwc, &
    out_2D_det_refreeze, point_numb, fname_p1, username, project_name)
    !*** Write the 2Ddetail output variables to a netcdf file !***
    
    ! declare arguments
    integer, intent(in) :: outputDetail, detlayers
    double precision, intent(in) :: detthick
    double precision, dimension(detlayers) :: DetDepth, DetDZ
    double precision, dimension((outputDetail+50),detlayers), intent(in) :: out_2D_det_dens, out_2D_det_temp, &
        out_2D_det_lwc, out_2D_det_refreeze
    character*255, intent(in) :: point_numb, fname_p1, username, project_name

    ! declare local arguments
    integer :: dd, status, ncid(50), IDs(50,5), varID(50,20)
    character*255 :: pad
    
    pad = "/ec/res4/scratch/"//trim(username)//"/"//trim(project_name)//"/output/"

    ncid(:) = 0
    IDs(:,:) = 0
    varID(:,:) = 0

    DetDZ(1) = detthick
    DetDepth(1) = DetDZ(1) / 2.
    do dd = 2, detlayers
      DetDZ(dd) = detthick
      DetDepth(dd) = DetDepth(dd-1) + DetDZ(dd)
    end do  

    ! CREATE NETCDF FILES
    status = nf90_create(trim(pad)//trim(fname_p1)//"_2Ddetail_"//trim(point_numb)//".nc",0,ncid(33))    

    ! DEFINE DIMENSIONS
    status = nf90_def_dim(ncid(33),"ind_t",outputDetail+50,IDs(33,1))
    if(status /= nf90_noerr) call Handle_Error(status,'2D_def_dim4')    
    status = nf90_def_dim(ncid(33),"layer",detlayers,IDs(33,2))
    if(status /= nf90_noerr) call Handle_Error(status,'2D_def_dim5')

    ! DEFINE VARIABLES
    status = nf90_def_var(ncid(33),"dens",nf90_real,(/IDs(33,1),IDs(33,2)/), &
        varID(33,1))
    if(status /= nf90_noerr) call Handle_Error(status,'2D_def_var1')
    status = nf90_def_var(ncid(33),"temp",nf90_real,(/IDs(33,1),IDs(33,2)/), &
        varID(33,2))
    if(status /= nf90_noerr) call Handle_Error(status,'2D_def_var2')
    status = nf90_def_var(ncid(33),"lwc",nf90_real,(/IDs(33,1),IDs(33,2)/), &
        varID(33,3))
    if(status /= nf90_noerr) call Handle_Error(status,'2D_def_var3')
    status = nf90_def_var(ncid(33),"depth",nf90_real,(/IDs(33,2)/),varID(33,4))
    if(status /= nf90_noerr) call Handle_Error(status,'2D_def_var4')
    status = nf90_def_var(ncid(33),"dz",nf90_real,(/IDs(33,2)/),varID(33,5))
    if(status /= nf90_noerr) call Handle_Error(status,'2D_def_var5')
    status = nf90_def_var(ncid(33),"refreeze",nf90_real,(/IDs(33,1),IDs(33,2)/), & 
        varID(33,6))
    if(status /= nf90_noerr) call Handle_Error(status,'2D_def_var1')
    
    ! DEFINE ATTRIBUTES (unit could also be defined here)
    status = nf90_put_att(ncid(33),varID(33,1),"missing_value",9.96921e+36)
    if(status /= nf90_noerr) call Handle_Error(status,'2D_def_att_miss_val_1')
    status = nf90_put_att(ncid(33),varID(33,2),"missing_value",9.96921e+36)
    if(status /= nf90_noerr) call Handle_Error(status,'2D_def_att_miss_val_2')
    status = nf90_put_att(ncid(33),varID(33,3),"missing_value",9.96921e+36)
    if(status /= nf90_noerr) call Handle_Error(status,'2D_def_att_miss_val_3')
    status = nf90_put_att(ncid(33),varID(33,4),"missing_value",9.96921e+36)
    if(status /= nf90_noerr) call Handle_Error(status,'2D_def_att_miss_val_4')    
    
    ! END OF DEFINING FILES
    status = nf90_enddef(ncid(33))
    if(status /= nf90_noerr) call Handle_Error(status,'2D_enddef')
    
    ! PUT IN THE CONSTANT Z-AXIS in 2Ddetail
    status = nf90_put_var(ncid(33),varID(33,4),DetDepth)
    if(status /= nf90_noerr) call Handle_Error(status,'2D_put_var1_1')
    status = nf90_put_var(ncid(33),varID(33,5),DetDZ)
    if(status /= nf90_noerr) call Handle_Error(status,'2D_put_var2_1')
    
    ! SAVE DATA
    status = nf90_put_var(ncid(33),varID(33,1),out_2D_det_dens, &
        start=(/1,1/),count=(/(outputDetail+50),detlayers/))
    if(status /= nf90_noerr) call Handle_Error(status,'2D_put_var_detail1')    
    status = nf90_put_var(ncid(33),varID(33,2),out_2D_det_temp, &
        start=(/1,1/),count=(/(outputDetail+50),detlayers/))
    if(status /= nf90_noerr) call Handle_Error(status,'2D_put_var_detail2')
    status = nf90_put_var(ncid(33),varID(33,3),out_2D_det_lwc, &
        start=(/1,1/),count=(/(outputDetail+50),detlayers/))
    if(status /= nf90_noerr) call Handle_Error(status,'2D_put_var_detail3')
    status = nf90_put_var(ncid(33),varID(33,6),out_2D_det_refreeze, &
        start=(/1,1/),count=(/(outputDetail+50),detlayers/))
    if(status /= nf90_noerr) call Handle_Error(status,'2D_put_var_detail6')

    ! CLOSE NETCDF-FILE
    status = nf90_close(ncid(33))
    if(status /= nf90_noerr) call Handle_Error(status,'2D_close')
    
end subroutine Save_out_2Ddetail

! *******************************************************


subroutine Save_out_spinup(ind_z_max, ind_z_surf, Rho, M, T, Depth, Mlwc, Year, point_numb, fname_p1, username, project_name)
    !*** Write a netcdf file containing the firn profile after the spin-up ***!
    
    integer :: status, ncid, dimID, varID(10), ind_z_max, ind_z_surf
    double precision, dimension(ind_z_max) :: M, T, Depth, Mlwc
    double precision, dimension(ind_z_max) :: Rho, Year
    character*255 :: pad, point_numb, fname_p1, username, project_name
    
    pad = "/ec/res4/scratch/"//trim(username)//"/restart/"//trim(project_name)//"/"

    ! CREATE NETCDF FILES
    status = nf90_create(trim(pad)//trim(fname_p1)//"_restart_from_spinup_"//trim(point_numb)//".nc",0,ncid)
    if(status /= nf90_noerr) call Handle_Error(status,'restart_create')

    ! DEFINE DIMENSIONS
    ! TODO: should ind_z_surf be ind_z_max?
    status = nf90_def_dim(ncid,"layer",ind_z_surf,dimID)
    if(status /= nf90_noerr) call Handle_Error(status,'restart_def_dim1')
    
    ! DEFINE VARIABLES
    status = nf90_def_var(ncid,"dens",nf90_real,dimID,varID(1))
    if(status /= nf90_noerr) call Handle_Error(status,'restart_def_var1')
    status = nf90_def_var(ncid,"temp",nf90_real,dimID,varID(2))
    if(status /= nf90_noerr) call Handle_Error(status,'restart_def_var2')
    status = nf90_def_var(ncid,"mass",nf90_real,dimID,varID(3))
    if(status /= nf90_noerr) call Handle_Error(status,'restart_def_var3')
    status = nf90_def_var(ncid,"depth",nf90_real,dimID,varID(4))
    if(status /= nf90_noerr) call Handle_Error(status,'restart_def_var4')
    status = nf90_def_var(ncid,"lwc",nf90_real,dimID,varID(5))
    if(status /= nf90_noerr) call Handle_Error(status,'restart_def_var5')
    status = nf90_def_var(ncid,"year",nf90_real,dimID,varID(6))
    if(status /= nf90_noerr) call Handle_Error(status,'restart_def_var6')    
    
    ! END OF DEFINING FILES
    status = nf90_enddef(ncid)
    if(status /= nf90_noerr) call Handle_Error(status,'restart_enddef')
    
    ! PUT VARIABLES
    status = nf90_put_var(ncid,varID(1),Rho(1:ind_z_surf))
    if(status /= nf90_noerr) call Handle_Error(status,'restart_put_var1')    
    status = nf90_put_var(ncid,varID(2),T(1:ind_z_surf))
    if(status /= nf90_noerr) call Handle_Error(status,'restart_put_var2')
    status = nf90_put_var(ncid,varID(3),M(1:ind_z_surf))
    if(status /= nf90_noerr) call Handle_Error(status,'restart_put_var3')
    status = nf90_put_var(ncid,varID(4),Depth(1:ind_z_surf))
    if(status /= nf90_noerr) call Handle_Error(status,'restart_put_var4')
    status = nf90_put_var(ncid,varID(5),Mlwc(1:ind_z_surf))
    if(status /= nf90_noerr) call Handle_Error(status,'restart_put_var5')
        
    status = nf90_put_var(ncid,varID(6),Year(1:ind_z_surf))
    if(status /= nf90_noerr) call Handle_Error(status,'restart_put_var6')
    
    ! CLOSE NETCDF FILE
    status = nf90_close(ncid)
    if(status /= nf90_noerr) call Handle_Error(status,'restart_close')
    
end subroutine Save_out_spinup

! *******************************************************


subroutine Save_out_run(Nt_model_tot, ind_z_max, ind_z_surf, Rho, M, T, Depth, Mlwc, Year, &
    DenRho, Refreeze, point_numb, fname_p1, username, project_name)
    !*** Write a netcdf file containing the firn profile after the spin-up ***!
    
    integer :: status, ncid, dimID(2), varID(10), ind_z_max, ind_z_surf, Nt_model_tot
    double precision, dimension(ind_z_max) :: M, T, Depth, Mlwc, DenRho, Refreeze
    double precision, dimension(ind_z_max) :: Rho, Year
    character*255 :: pad, point_numb, fname_p1, username, project_name
    
    pad = "/ec/res4/scratch/"//trim(username)//"/restart/"//trim(project_name)//"/"

    ! CREATE NETCDF FILES
    status = nf90_create(trim(pad)//trim(fname_p1)//"_restart_from_2023_run_"//trim(point_numb)//".nc",0,ncid)
    if(status /= nf90_noerr) call Handle_Error(status,'restart_run_create')

    ! DEFINE DIMENSIONS
    status = nf90_def_dim(ncid,"layer",ind_z_surf,dimID(1))
    if(status /= nf90_noerr) call Handle_Error(status,'restart_run_def_dim1')
    status = nf90_def_dim(ncid,"constant",1,dimID(2))
    if(status /= nf90_noerr) call Handle_Error(status,'restart_run_def_dim1')
    
    ! DEFINE VARIABLES
    status = nf90_def_var(ncid,"dens",nf90_real,dimID(1),varID(1))
    if(status /= nf90_noerr) call Handle_Error(status,'restart_run_def_var1')
    status = nf90_def_var(ncid,"temp",nf90_real,dimID(1),varID(2))
    if(status /= nf90_noerr) call Handle_Error(status,'restart_run_def_var2')
    status = nf90_def_var(ncid,"mass",nf90_real,dimID(1),varID(3))
    if(status /= nf90_noerr) call Handle_Error(status,'restart_run_def_var3')
    status = nf90_def_var(ncid,"depth",nf90_real,dimID(1),varID(4))
    if(status /= nf90_noerr) call Handle_Error(status,'restart_run_def_var4')
    status = nf90_def_var(ncid,"lwc",nf90_real,dimID(1),varID(5))
    if(status /= nf90_noerr) call Handle_Error(status,'restart_run_def_var5')
    status = nf90_def_var(ncid,"year",nf90_real,dimID(1),varID(6))
    if(status /= nf90_noerr) call Handle_Error(status,'restart_run_def_var6')
    status = nf90_def_var(ncid,"refreeze",nf90_real,dimID(1),varID(7))
    if(status /= nf90_noerr) call Handle_Error(status,'restart_run_def_var7')
    status = nf90_def_var(ncid,"denrho",nf90_real,dimID(1),varID(8))
    if(status /= nf90_noerr) call Handle_Error(status,'restart_run_def_var8')
    status = nf90_def_var(ncid,"prev_nt",nf90_real,dimID(2),varID(9))
    if(status /= nf90_noerr) call Handle_Error(status,'restart_run_def_var9')
    
    ! END OF DEFINING FILES
    status = nf90_enddef(ncid)
    if(status /= nf90_noerr) call Handle_Error(status,'restart_run_enddef')
    
    ! PUT VARIABLES
    status = nf90_put_var(ncid,varID(1),Rho(1:ind_z_surf))
    if(status /= nf90_noerr) call Handle_Error(status,'restart_run_put_var1')    
    status = nf90_put_var(ncid,varID(2),T(1:ind_z_surf))
    if(status /= nf90_noerr) call Handle_Error(status,'restart_run_put_var2')
    status = nf90_put_var(ncid,varID(3),M(1:ind_z_surf))
    if(status /= nf90_noerr) call Handle_Error(status,'restart_run_put_var3')
    status = nf90_put_var(ncid,varID(4),Depth(1:ind_z_surf))
    if(status /= nf90_noerr) call Handle_Error(status,'restart_run_put_var4')
    status = nf90_put_var(ncid,varID(5),Mlwc(1:ind_z_surf))
    if(status /= nf90_noerr) call Handle_Error(status,'restart_run_put_var5')
    status = nf90_put_var(ncid,varID(6),Year(1:ind_z_surf))
    if(status /= nf90_noerr) call Handle_Error(status,'restart_run_put_var6')
    status = nf90_put_var(ncid,varID(7),Refreeze(1:ind_z_surf))
    if(status /= nf90_noerr) call Handle_Error(status,'restart_run_put_var7')
    status = nf90_put_var(ncid,varID(8),DenRho(1:ind_z_surf))
    if(status /= nf90_noerr) call Handle_Error(status,'restart_run_put_var8')
    status = nf90_put_var(ncid,varID(9),Nt_model_tot)
    if(status /= nf90_noerr) call Handle_Error(status,'restart_run_put_var9')
    
    ! CLOSE NETCDF FILE
    status = nf90_close(ncid)
    if(status /= nf90_noerr) call Handle_Error(status,'restart_run_close')
    
end subroutine Save_out_run

end module output
