MODULE diffusivity

USE second_Precision,  ONLY : dp ! KPP Numerical type
USE constants

IMPLICIT NONE

CONTAINS

    SUBROUTINE diffusivity_coeff(fric_veloc,temp_s,RH_s,press_s,SHTF,press,temp,LWC_cloud,PBLH,RH,u_vel,v_vel,Kzc,tr)
! Diffusivity coeff calculated as in GDAS with a non-local calc in the mixed layer and a local diffusion scheme in the free atmosphere (Hong and Pan 1996,  Monthly weather review). Free atmosphere parameterization from Hong and Pan 1996. Some of the calculations are taken from dry deposition subroutine
        IMPLICIT NONE

        REAL(dp), INTENT(in) :: fric_veloc,temp_s,RH_s,press_s,SHTF,PBLH ! surface  GDAS parameters
        REAL(dp), DIMENSION(Nz+1), INTENT(in) :: temp,press,RH,u_vel,v_vel,LWC_cloud ! vertically resolved GDAS parameters
        REAL(dp), DIMENSION(Nz+1), INTENT(out) :: Kzc ! diffusivity coeff (heat/other scalars)

        REAL(dp) :: kar,g,R_d,R_v,L,ratio,p_0,T_0,p_sat,p_v_H2O,ratio_gasconst,mix_ratio,T_v,t_s_pot,C_p_dry,C_p,dens_air, flux_horiz,u,T_star,temp_prof_sl,wind_prof,temp_prof,vel_turb_m,vel_turb_h,len_sc,mix_len,t_pot_v_der, sta_func_m,sta_func_h,Prandtl,L_Ob, C_emp
        REAL(dp), DIMENSION(Nz+1) :: T_pot,p_sat_z,p_v_H2O_z,mix_ratio_z,temp_v,t_pot_v
        REAL(dp), DIMENSION(Nz+1) :: Kzm,Ri,windgradient,z3 ! diffusivity coeff (momentum) start at first_z+dz/2
        INTEGER :: i,tr
        
        z3=(z2(1:Nz+1)+z2(2:Nz+2))/2D0
        
        kar = 0.4 ! Von Karman
        g = 9.81 ! gravitational acceleration (m/s^2)
        R_d = 287.053 ! specific gas constant for dry air (J/(kgK))
        R_v = 461.495 ! specific gas constant for water vapor (J/(kgK))
        L = 2.5D6 ! latent heat of vaporization (J/kg)
        ratio = L/R_v ! (K)
        p_0 = 0.611 ! (kPa)
        T_0 = 273.15 ! (K)
        p_sat = p_0*exp(ratio*(1D0/T_0-1D0/temp_s)) ! saturation vapor pressure of water (kPa)
        p_v_H2O = ((RH_s/1D2)*p_sat)*1D3 ! vapor pressure of water (Pa)
        ratio_gasconst = 0.622 ! ratio of gas constant for dry air to that of water vapor (g_v/g_dry)
        mix_ratio = ratio_gasconst*p_v_H2O/(press_s-p_v_H2O) ! mixing ratio
        T_v = temp_s*(1D0+0.61*mix_ratio) ! Virtual temperature (K)
        t_s_pot = T_v*(1D5/press_s)**0.286 ! potential virtual temperature
        C_p_dry = 1004.67 ! specific heat (J/kgK)
        C_p = C_p_dry*(1D0+0.84*mix_ratio) ! specific heat, moist air (J/kgK)
        C_emp = 0.34 ! Empirical constant from LES-data

        dens_air = (press_s)/(R_d*T_v) ! Air density (kg/m^3) 
        !flux_horiz = dens_air*((-UMOF/dens_air)**2.+(-VMOF/dens_air)**2.)**0.5 ! Vertical turbulent flux of horizontal momentum (see HYSPLIT model description)
       !flux_horiz = ((-UMOF)**2.+(-VMOF)**2.)**0.5 ! Vertical turbulent flux of horizontal momentum (see HYSPLIT model description)
        u = fric_veloc  ! Friction velocity (m/s), Jacobson eq 8.8 with 8.2
   !     T_star = -SHTF/(C_p*dens_air*u) ! friction temperature (K) (see Stull, page 48 and page 356) SHTF/(C_p*dens_air) = (w't_s_pot')_surf which is negative for stable conditions, 0 for neutral and positive for unstable conditions
   !     L_Ob = u**2.*t_s_pot/(kar*g*T_star) ! Monin-Obukov length scale (Jacobson eq 8.32)
!
!        !Calculation of momentum profiles in the mixed layer, identical to those in the surface layer
!        IF (0.1*PBLH/L_Ob > 0) THEN !(temp_prof_sl <= 0) THEN
!           temp_prof = 1.+5.*0.1*PBLH/L_Ob ! temp profile function in stable atmos
!           wind_prof = 1.+5.*0.1*PBLH/L_Ob ! momentum profile function in stable atmos
!        ELSE ! unstable/neutral
!           temp_prof = (1.-16.*0.1*PBLH/L_Ob)**(-0.5)
!           wind_prof = (1.-12.*0.1*PBLH/L_Ob)**(-1./3.)  
!        END IF
!        Prandtl = temp_prof/wind_prof+7.8*kar*0.1 ! Turbulent Prandtl number in mixed BL  
!        vel_turb_m = u/wind_prof ! Turbulent velocity scale for momentum
!        vel_turb_h = vel_turb_m/Prandtl ! Turbulent velocity scale for heat (assumed to hold also for other scalars) This expression comes from the definition of Pr = Kzm/Kzc
!
!       len_sc = 30. ! asymptotic length scale (m)
!       p_sat_z = p_0*exp(ratio*(1D0/T_0-1D0/temp))
!       p_v_H2O_z = ((RH/1D2)*p_sat_z)*1D3 ! vapor pressure of water in each layer
!       mix_ratio_z = ratio_gasconst*p_v_H2O_z/(press-p_v_H2O_z) ! mixing ratio in each layer
!       temp_v = temp*(1D0+0.61*mix_ratio_z)
!       t_pot_v = temp_v*(1D5/press)**0.286 

   !     Kzm(1) = 0.
   !     Kzc(1) = 0.
   !     Ri(1) = 0.
   !     windgradient(1) = 0.
   !     DO i = 2,Nz+1  
   !        IF (z(i) <= PBLH) THEN ! in mixed-layer (Troen and Mahrt (1986), Holtslag et al (1990) and Holtslag and Boville (1993)
   !            Kzm(i) = kar*vel_turb_m*z(i)*(1.-z(i)/PBLH)**2. ! Diffusivity coeff (momentum) in mixed-layer between each layer
  !             Kzc(i) = kar*vel_turb_h*z(i)*(1.-z(i)/PBLH)**2. ! Diffusivity coeff (momentum) in mixed-layer between each layer
  !          IF (Kzc(i)<1D-2) Kzc(i)=1D-2 ! Lower limit of Kzc in the PBL
  !              
  !         ELSE IF (z(i) > PBLH) THEN ! free atmosphere as in Hong and Pan 1996 since a non-local scheme is used in the mixed layer
  !             windgradient(i) = (((u_vel(i)-u_vel(i-1))/(z3(i)-z3(i-1)))**2.+((v_vel(i)-v_vel(i-1))/(z3(i)-z3(i-1)))**2.)**0.5 ! horizontal wind vector in each layer
  !             t_pot_v_der = (t_pot_v(i)-t_pot_v(i-1))/(z3(i)-z3(i-1)) ! change of virtuell pot temp between each layer
  !             mix_len = kar*z(i)/(1.+kar*z(i)/len_sc) ! mixing length (m) bewteen each layer               
  !             IF (windgradient(i) == 0.) THEN
  !                 Kzm(i) = 0.
  !                 Kzc(i) = 0.
  !             ELSE
  !                 Ri(i) = (g/(t_pot_v(i)+t_pot_v(i-1)/2.)*t_pot_v_der)/windgradient(i)**2. ! Richardsons number in each level
  !                 IF (Ri(i) > 0) THEN ! stably stratified free atmos
  !                     sta_func_m = exp(-8.5*Ri(i))+0.15/(Ri(i)+3.)
  !                     sta_func_h = exp(-8.5*Ri(i))+0.15/(Ri(i)+3.)
  !                 ELSE ! unstable/neutral conditions
  !                     sta_func_m = (1.-16.*Ri(i))**(-0.25)
  !                     sta_func_h = (1.-16.*Ri(i))**(-0.5)                
  !                 END IF
 !                  Kzm(i) = mix_len**2*sta_func_m*windgradient(i) ! diffusivity coeff in free atmos between each layer
 !                  Kzc(i) = mix_len**2*sta_func_h*windgradient(i)
 !               END IF      
 !           END IF 
!
!         END DO  
!
!         where (Kzc>5D2) Kzc=5D2 

!	 ! Parameterization of vertical diffusion and the atmospheric boundary layer height determination in the EMEP model 

        Kzc = C_emp*u*z*exp(-0.5*(z/(0.21*PBLH))**2.) ! Grisogono K(z) scheme
		where (Kzc<1D-1) Kzc=1D-1 ! Set the lowest vertical turbulent diffusivity to 0.1 m^2/s (see Pisso et al., FLEXPART v10.4 2019) 
		where (LWC_cloud>1D-2 .AND. Kzc<1D0) Kzc=1D0 ! Set the lowest vertical turbulent diffusivity to 1 m^2/s in grid cells with clouds
        ! Lower and nearly constant mixing inside canopy below the displacement height 2/3*canopy_height:
        !Kzc(2:3) = 0.1*Kzc(2:3) ! At 3 and 9 m the Kz is constant and only 10 % of the value at 18 m
        

    END SUBROUTINE diffusivity_coeff

END MODULE diffusivity
