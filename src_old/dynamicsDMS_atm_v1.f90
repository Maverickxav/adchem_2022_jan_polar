MODULE dynamicsDMS_atm
    USE constants
    USE second_Precision, ONLY : dp    ! KPP Numerical type
	USE second_Parameters !  NSPEC

  USE acdc_datatypes
    !USE cloud_chem, ONLY : cloud_aq
	!USE particle_halogen_chem, ONLY : halogen_aq
	
!    USE second_Parameters, !ONLY : NSPEC, ind_H2SO4, ind_HNO3

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: dry_dep_particles, dNdlogDp, fraction_POA_marine, &
 condensation, coagulation,coagulation_clusterin, wet_deposition, nucleation, sea_spray, &
 sea_salt_ice, dry_dep_gases, wet_deposition_gas,full_stationary_rebinning!, particle_halog_chem!, cloud_activation_processing


CONTAINS

    !-------------------------------------------------------------!
    ! Particle size distribution                                  !
    !-------------------------------------------------------------!
    SUBROUTINE dNdlogDp(d_g,dp_dry,dlogDp,modes,s,N_modes,dm,N_bins,V_bins)
      INTEGER :: ii
      INTEGER, INTENT(IN) :: modes
      Real(dp), DIMENSION(modes), INTENT(IN) :: s,N_modes,dm
      Real(dp), DIMENSION(modes,nr_bins)  :: dNdlogDp_modes
      Real(dp), DIMENSION(nr_bins), INTENT(IN) :: d_g,dp_dry,dlogDp
      Real(dp), DIMENSION(nr_bins) :: Vp
      Real(dp), DIMENSION(nr_bins), INTENT(OUT) :: N_bins,V_bins

      DO ii = 1,modes
        !dNdlogDp_modes(ii,:) = ( N_modes(ii) / sqrt(2D0*pi) / LOG10(s(ii)))*&       
        !EXP(-((LOG10(d_g) - log10(dm(ii)))**2 ) / (2*LOG10(s(ii))**2 ))
		
		 dNdlogDp_modes(ii,:) = ( N_modes(ii) / (sqrt(2D0*pi)*LOG10(s(ii))))*&       
        EXP(-((LOG10(d_g) - log10(dm(ii)))**2 ) / (2D0*LOG10(s(ii))**2 ))
      END DO
    
      !N_bins=SUM(dNdlogDp_modes, DIM=1)*dlogDp+1D-3
      N_bins=SUM(dNdlogDp_modes, DIM=1)*dlogDp
      Vp=dp_dry**3.0*pi/6.0
      V_bins=N_bins*Vp
 
    END SUBROUTINE dNdlogDp

    !-------------------------------------------------------------!
    ! Nucleation                                                  !
    !-------------------------------------------------------------!
    SUBROUTINE nucleation(cH2SO4,VOC_nucl,N_bins,V_bins,dens_p,d_p,dp_dry,c_p,c_p_nucl,MX,&
    qX,dt,temp,RH,q_ion,cHOM,CS_air,Jnucl_SAorg,Jnucl_org)

      REAL(dp), DIMENSION(NSPEC_P,nr_bins), INTENT(inout) :: c_p
      REAL(dp), DIMENSION(nr_bins), INTENT(inout) :: N_bins,V_bins,dens_p,d_p,dp_dry
      REAL(dp), INTENT(in) :: cH2SO4,VOC_nucl,dt,temp,RH,q_ion,cHOM,CS_air
      REAL(dp), DIMENSION(NSPEC_P), INTENT(in) :: MX,qX
      REAL(dp), DIMENSION(NSPEC_P), INTENT(in)  :: c_p_nucl
      REAL(dp) :: h2So4_frac,a_coeff,b_coeff,c_coeff,d_coeff,e_coeff,f_coeff,g_coeff,h_coeff, i_coeff,j_coeff
      REAL(dp), INTENT(out)  :: Jnucl_SAorg, Jnucl_org
      REAL(dp) :: Jnucl_SA, Jnucl_tot
      REAL(dp) :: vp_dry, vp_wet
      REAL(dp) :: dS, dG, dH, n_ion, alpha_ion, k_ion, a1,a2,a3,a4,a5
      INTEGER  :: j

     IF (index_bin_nucl == 0) THEN
     Jnucl_SA = 0.    
     ELSE ! paramterization for binary sulfuric acid-water nucl (Vehkamäki et al 2002)
      IF (cH2SO4 .gt. 1d4 .and. cH2SO4 .lt. 1d11) THEN
       h2so4_frac = 0.740997 - 0.00266379*temp &
    - 0.00349998*log(cH2SO4)+ 0.0000504022*temp*log(cH2SO4) &
    + 0.00201048*log(RH/1d2)- 0.000183289*temp*log(RH/1d2) &
    + 0.00157407*(log(RH/1d2))**2.-0.0000179059*temp*(log(RH/1d2))**2. &
    + 0.000184403*(log(RH/1d2))**3.-1.50345*1d-6*temp*(log(RH/1d2))**3. ! mole fraction of sulfuric acid in the critical cluster
       a_coeff = 0.14309+2.21956*temp-0.0273911*temp**2. &
    + 0.0000722811*temp**3.+5.91822/h2so4_frac
       b_coeff = 0.117489+0.462532*temp-0.0118059*temp**2 &
    + 0.0000404196*temp**3. +15.7963/h2so4_frac
       c_coeff = -0.215554-0.0810269*temp+0.00143581*temp**2. &
    - 4.7758*1d-6*temp**3.-2.91297/h2so4_frac
       d_coeff = -3.58856+0.049508*temp-0.00021382*temp**2. &
    + 3.10801*1d-7*temp**3.-0.0293333/h2so4_frac
       e_coeff = 1.14598-0.600796*temp+0.00864245*temp**2. &
    - 0.0000228947*temp**3.-8.44985/h2so4_frac
       f_coeff = 2.15855+0.0808121*temp-0.000407382*temp**2. &
    - 4.01957*1d-7*temp**3.+0.721326/h2so4_frac
       g_coeff = 1.6241-0.0160106*temp+0.0000377124*temp**2. &
    + 3.21794*1d-8*temp**3.-0.0113255/h2so4_frac
       h_coeff = 9.71682-0.115048*temp+0.000157098*temp**2. &
    + 4.00914*1d-7*temp**3.+0.71186/h2so4_frac
       i_coeff = -1.05611+0.00903378*temp-0.0000198417*temp**2. &
    + 2.46048*1d-8*temp**3.-0.0579087/h2so4_frac
       j_coeff = -0.148712+0.00283508*temp-9.24619*1d-6*temp**2. &
    + 5.00427*1d-9*temp**3.-0.0127081/h2so4_frac
       ! from above coefficients the nucl rate is given by an exponential of a third order polynomial of log(RH/100) and log(cH2so4) (in particle/cm³):
       Jnucl_SA = exp(a_coeff+b_coeff*log(RH/1d2) &
    + c_coeff*(log(RH/1d2))**2.+d_coeff*(log(RH/1d2))**3. &
    + e_coeff*log(cH2SO4)+f_coeff*log(RH/1d2)*log(cH2SO4) &
    + g_coeff*(log(RH/1d2))**2.*log(cH2SO4)+h_coeff*(log(cH2SO4))**2. &
    + i_coeff*log(RH/1d2)*(log(cH2SO4))**2.+j_coeff*(log(cH2SO4))**3.)
      ELSE
       Jnucl_SA = 0.
      END IF 
     END IF
     
     dS=-61.1 ! cal/mol/K  J/mol (Entropy change upon cluster formation between H2SO4 and a carboxylic acid (Elm et al., 2017)) 
     dG=-15100 ! cal/mol J/mol (Lowes possible Gibbs free energy of cluster formation between a di-carboxylic acid and H2SO4 (Elm et al., 2017)) 
     dH=(dG+dS*temp)*4.184 ! J/mol Entalpy change upon cluster formation between a di-carboxylic acid and H2SO4
     !Jnucl_SAorg=5D-13*exp(-dH/Rg*(1D0/temp-1D0/298D0))*VOC_nucl*cH2SO4 ! Estimated T dependent kinetic nucleation rate (J=k(T)*H2SO4x[BiOxOrg]
     Jnucl_SAorg=1D-12*exp(-dH/Rg*(1D0/temp-1D0/298D0))*VOC_nucl*cH2SO4 ! Estimated T dependent kinetic nucleation rate (J=k(T)*H2SO4x[BiOxOrg]
     
     alpha_ion=1.6D-6 ! ion-ion recombination coefficient (cm^-3 s^-1)
     Jnucl_org = 0D0 ! First guess
     n_ion = 1D0 ! First guess ion concentration (cm^-3)
     a1=0.04001D0
     a2=1.848D0
     a3=0.001366D0
     a4=1.566D0
     a5=0.1863D0
     DO j=1,5
     k_ion=CS_air+Jnucl_org/(2D0*n_ion) ! ion loss rate (s^-1)
     n_ion=((k_ion**2D0+4D0*alpha_ion*q_ion)**0.5D0-k_ion)/(2D0*alpha_ion)
     Jnucl_org=2D0*n_ion*a3*(cHOM/1D7)**(a4+a5/(cHOM/1D7))
     END DO
     
     Jnucl_org=exp(-dH/Rg*(1D0/temp-1D0/278D0))*(Jnucl_org+a1*(cHOM/1D7)**(a2+a5/(cHOM/1D7))) ! cm^-3 s^-1
     
      Jnucl_tot = Jnucl_org+Jnucl_SAorg+1d-6 ! cm^-3 s^-1
!      write(*,*) Jnucl

        N_bins(5) = N_bins(5)+Jnucl_tot*1D6*dt
        c_p(:,5) = c_p(:,5)+c_p_nucl*Jnucl_tot*dt
        vp_dry=SUM(c_p(index_dry,5)/Na*MX(index_dry)/qX(index_dry)/(N_bins(5)*1D-6)) ! m^3
        vp_wet=SUM(c_p(:,5)/Na*MX/qX/(N_bins(5)*1D-6)) ! m^3
        V_bins(5)=N_bins(5)*vp_wet
        dens_p(5)=SUM(c_p(:,5)*MX*1D6)/Na/V_bins(5) ! Total particle density
        dp_dry(5)=(vp_dry*6D0/pi)**(1D0/3D0) ! Dry particle diameters
        d_p(5)=(vp_wet*6D0/pi)**(1D0/3D0) ! Particle diameters
 
    END SUBROUTINE nucleation

    !-----------------------------------------------------------!
    ! Dry deposition for particles                              !
    !-----------------------------------------------------------!
    SUBROUTINE dry_dep_particles(landuse_index,d_p,dens_p,temp,press,RH,u,PBLH,SHTF,windspeed,tr,v_dep,vs,season_index)

      INTEGER, INTENT(in) :: landuse_index
      REAL(dp), DIMENSION(nr_bins), INTENT(in) :: d_p, dens_p ! arithmetic mean diameter and dry particle densities
      REAL(dp), INTENT(in) :: temp,press,RH,PBLH,SHTF,windspeed,u ! meteorological parameters:
        ! temp [K],pressure [Pa],RH,u-resp v-comp of mom flux [N/m^2],boundary layer height [m], sensible heat net flux at surface [W/m^2], windspeed at 10 m ([m/s] (GDAS)
        ! and time along trajectory in hour and trajectory time step
      INTEGER, INTENT(in)  :: tr,season_index
      REAL(dp), DIMENSION(nr_bins), INTENT(out) :: v_dep,vs ! deposition loss-rate coefficient for each particels size bin [m/h]
      ! empirical coefficients
      REAL(dp), DIMENSION(8,5) :: z0
      REAL(dp), DIMENSION(5,5) :: r_coll
      REAL(dp), DIMENSION(5)   :: a_landuse
      REAL(dp), DIMENSION(5)   :: j_landuse
      ! Variables used to calculate the air density
      REAL(dp) :: R_d,R_v,L,ratio,p_0,T_0,p_sat,p_v_H2O,ratio_gasconst,mix_ratio,T_v,dens_air
      ! Variables used to calculate friction velocity
      !REAL(dp) :: flux_horiz,u
      ! Variables needed to calculate the Richardson nr
      REAL(dp) :: zr,z_rough,t_s,t_s_pot,kar,g,C_p_dry,C_p,T_star,L_Ob,Rfr,Rf0
      ! Variables used to calculate aerodynamic resistance
      REAL(dp) :: y0,yr,ra,Pr,beta,gam,denominator
      ! Variables used to calculate the qausi-laminar resistance
      REAL(dp) :: v,mu,k_b,q,Cc,D,Sc,St,R1,rb,th_speed,Kn
        
      INTEGER  :: i

      ! Rougness length for diffreent land-use catagories and seasons, also table 19.2 [m]
      z0(1,:) = (/ 0.8,0.9,0.9,0.9,0.8 /) ! evergreen, needleleaf trees
      z0(2,:) = (/ 1.05,1.05,0.95,0.55,0.75 /) ! deciduous broadleaf trees
      z0(3,:) = (/ 0.1,0.1,0.05,0.02,0.05 /) ! grass
      z0(4,:) = (/ 0.04,0.04,0.04,0.04,0.04 /) ! desert
      z0(5,:) = (/ 0.1,0.1,0.1,0.1,0.1 /) ! shrubs and interrupted woodlands
      z0(6,:) = (/ 1D-5,1D-5,1D-5,1D-5,1D-5 /) ! sea, from table 8.1 in Jacobson
      z0(7,:) = (/ 0.26,0.26,0.26,0.26,0.26 /) ! small urban area, from table 8.1 in Jacobson
      z0(8,:) = (/ 1D-3,1D-3,1D-3,1D-3,1D-3 /) ! Snow and ice
      

      ! radius of collector for diffreent land-use catagories and seasons, also table 19.2 [m]
      r_coll(1,:) = (/ 2d-3,2d-3,2d-3,2d-3,2d-3 /) ! evergreen, needleleaf trees
      r_coll(2,:) = (/ 5d-3,5d-3,10d-3,10d-3,5d-3 /) ! deciduous broadleaf trees
      r_coll(3,:) = (/ 2d-3,2d-3,5d-3,5d-3,2d-3 /) ! grass
      r_coll(4,:) = (/ 9999.,9999.,9999.,9999.,9999. /) ! NaN
      r_coll(5,:) = (/ 10d-3,10d-3,10d-3,10d-3,10d-3 /) ! shrubs and interrupted woodlands
    
      ! coefficients based on land use categories, also table 19.2 
      a_landuse = (/ 1.,0.8,1.2,50.,1.3 /)
      !j_landuse = (/ 0.56,0.56,0.54,0.54,0.54 /)
      
      j_landuse = 2D0/3D0 ! E.-Y Nho-Kim et al 2004, Atmos Envi
    
      ! Calculation of air density using the Clasius-Clapeyron eq to calculate the saturation vapor pressure 
      R_d = 287.053 ! specific gas constant for dry air (J/(kgK))
      R_v = 461.495 ! specific gas constant for water vapor (J/(kgK))
      L = 2.5D6 ! latent heat of vaporization (J/kg)
      ratio = L/R_v ! (K)
      p_0 = 0.611 ! (kPa)
      T_0 = 273.15 ! (K)
      p_sat = p_0*exp(ratio*(1D0/T_0-1D0/temp)) ! saturation vapor pressure of water (kPa)
      p_v_H2O = ((RH/1D2)*p_sat)*1D3 ! vapor pressure of water (Pa)
      ratio_gasconst = 0.622 ! ratio of gas constant for dry air to that of water vapor (g_v/g_dry)
      mix_ratio = ratio_gasconst*p_v_H2O/(press-p_v_H2O) ! mixing ratio
      T_v = temp*(1D0+0.61*mix_ratio) ! Virtual temperature (K)
      dens_air = (press)/(R_d*T_v) ! Air density (kg/m^3)
  
      ! Calculation of the friction velocity using vertical turbulent flux of
      ! horizontal momentum (as computed in HYSPLIT)
      !flux_horiz = dens_air*((-UMOF/dens_air)**2.+(-VMOF/dens_air)**2.)**0.5 ! Vertical turbulent flux of horizontal momentum
      !flux_horiz = ((-UMOF)**2.+(-VMOF)**2.)**0.5 ! Vertical turbulent flux of horizontal momentum (see HYSPLIT model description)
      !u = (flux_horiz/dens_air)**0.5;            ! Friction velocity (m/s)
        
      ! Calculation of Rfr = z/L (L is the Obukhov length) computed from the friction values to be consistent with model derived flux fields
      !z_rough = z0(landuse_index,season_index)
      zr = 10. ! reference height where dry dep occur. ?? 0.1.*PBLH; % estimated height of the surface layer (bottom 10 % of PBL) (m)
!     IF (zr>1D2) THEN ! if using zr = 0.1*PBLH
!       zr=1D2
!     END IF
  
      t_s = T_v ! temperature 2 m agl taken as temperature at z_s (K)
      t_s_pot = T_v*(1D5/press)**0.286
      kar = 0.4 ! Von Karman constant
      g = 9.81 ! gravitational acceleration (m/s^2)
      C_p_dry = 1004.67 ! specific heat (J/kgK)
      C_p = C_p_dry*(1D0+0.84*mix_ratio) ! specific heat, moist air (J/kgK)
      T_star = -SHTF/(C_p*dens_air*u) ! friction temperature (K)
      L_Ob = u**2.*t_s_pot/(kar*g*T_star)

      Rfr = zr/L_Ob  ! z/L; if > 0 stable, if < 0 unstable, if = 0 neutral
      !Rf0 = z_rough/L_Ob  ! Richards number ground
        
      ! Calculation of aerodynamic resistance (Seinfeld and Pandis, 2006, p 907), y0 and yr becomes imagenary when (1-15*Richards number) becomes negative!!
!     y0=(1D0-15.*Rf0)**0.25
!     yr=(1D0-15.*Rfr)**0.25
!     IF (Rfr >= 0) THEN
!       ra=(1D0/(kar*u))*(log(zr/z0)+4.7*(Rfr-Rf0)) ! aerodynamic resistance stable atmosphere
!     ELSE 
!       ra=(1D0/(kar*u))*(log(zr/z0)+log((y0**2+1D0)*(y0+1D0)**2/((yr**2+1D0)*(yr+1D0)**2))+2.*(atan(yr)-atan(y0))) ! aerodynamic resistance unstable atmosphere
!     END IF

      ! Calculation of aerodynamic resistance (Jacobson  8.4.2.3 and p.667 eq 20.12) 
      Pr = 0.95 ! turbulent Prandtl number (when kar = 0.4 (Hogstrom, 1988))
      beta = 7.8 ! when kar = 0.4 (Hogstrom, 1988)
      gam = 11.6 ! when kar = 0.4 (Hogstrom, 1988)
      denominator = kar*u

  ! Caclulate ra
  z_rough = z0(landuse_index,season_index) !D/(kar*u) ! surface roughness length of particle ( characterize the ability of elements protruding from the surface to absorb particles) (m) eq 8.13
 IF (Rfr > 1D-3) THEN
   ra = (Pr*log(zr/z_rough)+beta/L_Ob*(zr-z_rough))/denominator ! stable
 ELSE IF (Rfr<1D-3 .AND. Rfr>-1D-3) THEN
   ra = Pr*log(zr/z_rough)/denominator ! neutral
 ELSE
   ra = (Pr*(log(((1.-gam*zr/L_Ob)**0.5-1.)/((1.-gam*zr/L_Ob)**0.5+1.))-log(((1.-gam*z_rough/L_Ob)**0.5-1.)/((1.-gam*z_rough/L_Ob)**0.5+1.))))/(kar*u) ! unstable
 END IF
 
      ! Quasi-laminar resistance (rb) and deposition rate (v_dep)
      DO i = 1,nr_bins
        mu = 1.8325D-5*(416.16/(temp+120.))*(temp/296.16)**1.5 ! dynamic viscosity of air (kg/(ms)) eq 4.54 Jacobson
        th_speed = SQRT(8.*kb*temp/(pi*4.8096D-26)) ! average thermal speed of an air molecule (m/s), eq 2.3 Jacobson
        q = 2.*mu/(dens_air*th_speed) ! mean free path of an air molecule eq 15.24 (m)
        Kn = q/(d_p(i)/2.) ! Knudsen nr for air, eq 15.23         
        v = mu/dens_air ! (m²/s) kinematic viscosity of air 
        Cc = 1D0+Kn*(1.249+0.42*EXP(-0.87/Kn)) ! Cunningham slip-flow correction, when reynold nr is less than 0.01 (particle diameter less than 20 um) the particles are smaller than mean free path and flow needs to be corrected
        D = kb*temp*Cc/(3.*pi*mu*d_p(i)) ! Particle diffusion coefficient (m²/s) eq 15.29
        Sc = v/D ! Schmidt number
        vs(i) = d_p(i)**2.*(dens_p(i)-dens_air)*Cc/(9.*mu) ! sedimentation velosity, eq 20.4 (m/s)
        IF (landuse_index == 6 .OR. landuse_index == 8) THEN ! Slinn and Slinn 1980
!   rb = u**2./(kar*windspeed)*(Sc**(-0.5)+10.**(-3./(vs*u**2./(v*g))))
         rb = 1./(u*(Sc**(-0.5)+10.**(-3./(vs(i)*u**2./(g*v))))) ! E.-Y Nho-Kim et al 2004, Atmos Envi
          v_dep(i) = 1D0/(ra+rb+ra*rb*vs(i))+vs(i) ! over ocean
        ELSE IF (landuse_index == 7) THEN
          v_dep(i) = 1./ra+vs(i) ! no dry dep via rb urban areas
        ELSE 
          St = vs(i)*u/(g*r_coll(landuse_index,season_index)) ! Stokes number vegetation
          R1 = EXP(-(St)**(0.5)) ! fraction of particles, once in contact, that sticks to the surface
          rb = 1D0/(3.*u*R1*(Sc**(-j_landuse(landuse_index))+ & 
          (St/(a_landuse(landuse_index)+St))**2+0.5*(d_p(i)/r_coll(landuse_index,season_index))**2))
          v_dep(i) =1D0/(ra+rb+ra*rb*vs(i))+vs(i) ! m/s
        END IF    
      END DO
    END SUBROUTINE dry_dep_particles
	

    !-------------------------------------------------------------!
    ! Wet deposition for particles                                !
    !-------------------------------------------------------------!
    SUBROUTINE wet_deposition(d_p,rain,v_wet)
    ! Below cloud scavenging according to parameterization from Laakso et al 2013

      REAL(dp), DIMENSION(nr_bins), INTENT(in) :: d_p
      REAL(dp), INTENT(in)        :: rain ! rainfall intensity (mm/h)

      REAL(dp), DIMENSION(nr_bins), INTENT(out) :: v_wet

      REAL(dp), DIMENSION(nr_bins) :: x, d_p_wet
      REAL(dp)      :: a,b,c,d,e,f
 
      d_p_wet = d_p
      WHERE (d_p_wet > 1D-5) d_p_wet = 1D-5 ! Parameterization not valid for large d_p
      IF (rain > 0.1) THEN
        a = 274.35758
        b = 332839.59273
        c = 226656.57259
        d = 58005.91340
        e = 6588.38582
        f = 0.244984
        x = a+b/log10(d_p_wet)**4.+c/log10(d_p_wet)**3.+d/log10(d_p_wet)**2.+e/log10(d_p_wet)+f*rain**0.5
        v_wet = 10.**x ! wet deposition coefficients^-1;
      ELSE
        v_wet = d_p*0.    
      END IF 

    END SUBROUTINE wet_deposition
    
    ! Wet deposition for gases (based on the EMEP EMEP MSC-W chemical transport model– technical description)
    SUBROUTINE wet_deposition_gas(rain,LWC,conc,dt)
    REAL(dp), INTENT(in)        :: rain,dt ! rainfall intensity (mm/h)
    REAL(dp), DIMENSION(Nz), INTENT(in)        :: LWC ! %
    REAL(dp), DIMENSION(21) :: Scaveng_gas
    REAL(dp), DIMENSION(21) :: Win,Wsub
    REAL(dp) :: Precipitation, hs
	REAL(dp), DIMENSION(Nz,NSPEC), INTENT(inout) :: conc
	INTEGER :: j
    Precipitation=rain/3.6D3 ! Presipitation rate in kg m^-2 s^-1
    hs=1000 ! characteristic scavenging depth (assumed to be 1000 m)
    !     SO2  HNO3 HONO NH3 H2O2 HCHO H2SO4 MSA HCl DMSO DMSO2 MSIA HIO3 HOCl HBr HOBr HI HOI INO2 INO3 DMA
    Win=(/0.3, 1.4, 1.4, 1.4, 1.4, 0.1, 1.4, 1.4, 1.4, 0.1, 0.1, 0.1, 1.4, 0.1, 1.4, 0.1, 1.4, 1.4, 1.4, 1.4, 1.4/)*1D6
    Wsub=(/0.15, 0.5, 0.5, 0.5, 0.5, 0.03, 0.5, 0.5, 0.5, 0.03, 0.03, 0.03, 0.5, 0.03, 0.5, 0.03, 0.5, 0.5, 0.5, 0.5, 0.5/)*1D6
    
	DO j=1,Nz
    IF (LWC(j)>=1D-2) THEN
    Scaveng_gas=Win*Precipitation/(hs*qH2O) ! s^-1
    ELSE
    Scaveng_gas=Wsub*Precipitation/(hs*qH2O) ! s^-1
    END IF
    
	conc(j,ind_SO2)=conc(j,ind_SO2)*EXP(-Scaveng_gas(1)*dt)
	conc(j,ind_HNO3)=conc(j,ind_HNO3)*EXP(-Scaveng_gas(2)*dt)
	conc(j,ind_HONO)=conc(j,ind_HONO)*EXP(-Scaveng_gas(3)*dt)
	conc(j,ind_NH3)=conc(j,ind_NH3)*EXP(-Scaveng_gas(4)*dt)
	conc(j,ind_H2O2)=conc(j,ind_H2O2)*EXP(-Scaveng_gas(5)*dt)
	conc(j,ind_HCHO)=conc(j,ind_HCHO)*EXP(-Scaveng_gas(6)*dt)
	conc(j,ind_H2SO4)=conc(j,ind_H2SO4)*EXP(-Scaveng_gas(7)*dt)
	conc(j,ind_MSA)=conc(j,ind_MSA)*EXP(-Scaveng_gas(8)*dt)
	conc(j,ind_HCl)=conc(j,ind_HCl)*EXP(-Scaveng_gas(9)*dt)
	conc(j,ind_DMSO)=conc(j,ind_DMSO)*EXP(-Scaveng_gas(10)*dt)
	conc(j,ind_DMSO2)=conc(j,ind_DMSO2)*EXP(-Scaveng_gas(11)*dt)
	conc(j,ind_MSIA)=conc(j,ind_MSIA)*EXP(-Scaveng_gas(12)*dt)
	conc(j,ind_HIO3)=conc(j,ind_HIO3)*EXP(-Scaveng_gas(13)*dt)
	conc(j,ind_HOCL)=conc(j,ind_HOCL)*EXP(-Scaveng_gas(14)*dt)
    conc(j,ind_HBr)=conc(j,ind_HBr)*EXP(-Scaveng_gas(15)*dt)
    conc(j,ind_HOBr)=conc(j,ind_HOBr)*EXP(-Scaveng_gas(16)*dt)
    conc(j,ind_HI)=conc(j,ind_HI)*EXP(-Scaveng_gas(17)*dt)
    conc(j,ind_HOI)=conc(j,ind_HOI)*EXP(-Scaveng_gas(18)*dt)
    conc(j,ind_INO2)=conc(j,ind_INO2)*EXP(-Scaveng_gas(19)*dt)
    conc(j,ind_INO3)=conc(j,ind_INO3)*EXP(-Scaveng_gas(20)*dt)
	conc(j,ind_DMA)=conc(j,ind_DMA)*EXP(-Scaveng_gas(21)*dt)
	END DO
	
    END SUBROUTINE wet_deposition_gas
    
    ! IN cloud SO2 and H2O2 dissolution and S(VI) formation:
    
    ! SUBROUTINE cloud_activation_processing(c_p,N_bins,SO2,H2O2,DMS,DMSO,MSIA,HPMTF,OH,O3,&
	! cDMSaq,cDMSOaq,cMSIAaq,cO3aq,cOHaq,cHPMTFaq,cH2O2aq,cSO2aq,T,S_c,S_super,I_nSVI,I_nCH3SO3H,LWC,dt)
    ! REAL(dp), DIMENSION(Nlayers,NSPEC_P,nr_bins), INTENT(in) :: c_p
    ! REAL(dp), DIMENSION(nr_bins), INTENT(in) :: N_bins
    ! REAL(dp), DIMENSION(nr_bins), INTENT(out) :: S_c
    ! REAL(dp), INTENT(out) :: S_super, I_nSVI, I_nCH3SO3H
    ! REAL(dp), INTENT(in) :: T,LWC,dt
    ! REAL(dp), DIMENSION(nr_bins) :: ns, B_Kohler, vp0,dp0, alpha1, alpha2, c_H
    ! REAL(dp), INTENT(inout) :: SO2, H2O2, DMS, DMSO, MSIA,HPMTF,OH,O3
	! REAL(dp), INTENT(inout) :: cDMSaq,cDMSOaq,cMSIAaq,cHPMTFaq,cH2O2aq,cSO2aq,cO3aq,cOHaq
	! REAL(dp)                :: M_DMSaq,M_DMSOaq,M_MSIAaq,M_HPMTFaq,M_H2O2aq,M_SO2aq,M_O3aq,M_OHaq
	
    ! REAL(dp) :: CCN, dp_drops, A_Kohler, surf_tens, m_drops, cH,r_drops,L_cloud_vol_frac,&
	! speedOH,mOH,aOH,DOH,pOH,dOHaqdt,dpOHdt,speedO3,mO3,aO3,DO3,pO3,dO3aqdt,dpO3dt,dpDMSdt,dpDMSOdt,&
	! dpMSIAdt,dDMSaqdt,dDMSOaqdt,dMSIAaqdt,dt_aq,Lcloud_DMS,Lcloud_DMSO,Lcloud_MSIA,Pcloud_DMSO,Pcloud_MSIA,&
	! Lcloud_OH,Lcloud_O3,Pcloud_MSA
    ! REAL(dp) :: Ks1,Ks2,H_SO2,H_H2O2,H_S,H_DMS,H_DMSO,H_MSIA,H_O3,H_OH,kSIV,K_SIV,xHSO3,&
    ! k_DMS_O3,k_DMS_OH,k_DMSO_OH,k_MSIA_O3,k_MSIA_OH,pSO2,pH2O2,pDMS,pDMSO,pMSIA,pHPMTF,&
	! kmtOH,kmtO3,kmtDMS,kmtDMSO,kmtMSIA,DDMS,DDMSO,DMSIA,aDMS,aDMSO,aMSIA,speedDMS,speedDMSO,speedMSIA,mDMS,mDMSO,mMSIA
	! REAL(dp), DIMENSION(NSPECAQ) :: conc_aq
    ! INTEGER :: i,ii
	
    ! S_super = 2D-3 ! Assumed maximum supersaturation in the clouds (fraction)
    ! cH = 1D-4 ! Assumed cloud droplet acidity (mol H+ dm^-3, pH = 4)
    
    ! surf_tens=(76.1-0.155*(T-273.15))*1D-3 ! Surface tension of water (kg s^-2) 
    
	! c_H=SUM(SUM(c_p(:,1:3,:), DIM=1),DIM=1)+SUM(c_p(:,9,:), DIM=1)-SUM(SUM(c_p(:,4:5,:), DIM=1),DIM=1) ! Approximate H+ concentration
	! where (c_H<0D0) c_H=0D0 ! Possible if some S(VI) is in the form of SO4^(2-)
	! ns=(SUM(SUM(c_p(:,1:5,:), DIM=1),DIM=1)+SUM(SUM(c_p(:,9:NSPEC_P,:),DIM=1),DIM=1)+c_H)/(N_bins*1D-6)/Na ! Moles of soluble material in each single particle inorganic ions + all organic compounds but not soot
	    
    ! vp0=SUM(c_p(:,6,:)/Na*MEC/qEC)/(N_bins*1D-6) ! Equivalent volume of insoluble core (considered to be soot only)
    
    ! dp0=(vp0*6D0/pi)**(1D0/3D0) ! Equivalent diameter of insoluble core (considered to be soot only)
     
    ! A_Kohler=4*MH2O*surf_tens/(Rg*T*qH2O)
    ! B_Kohler=6*ns*MH2O/(pi*qH2O)
    
    ! ! Calculate the critical supersaturation for droplet activation according to Kokkola et al., 2008 
    ! alpha1=12*(81D0*dp0**6D0+12D0*dp0**3D0*(3D0*B_Kohler/A_Kohler)**(3D0/2D0))**(1D0/2D0)
    ! alpha2=(108D0*dp0**3D0+8*(3D0*B_Kohler/A_Kohler)**(3D0/2D0)+alpha1)**(1D0/3D0)
    ! S_c= EXP( A_Kohler/(alpha2/6D0 + 2D0/3D0*(3D0*B_Kohler/A_Kohler)*1D0/alpha2+&
    ! 1D0/3D0*(3D0*B_Kohler/A_Kohler)**(1D0/2D0))+&
    ! B_Kohler/((alpha2/6D0 + 2D0/3D0*(3D0*B_Kohler/A_Kohler)*1D0/alpha2+&
    ! 1D0/3D0*(3D0*B_Kohler/A_Kohler)**(1D0/2D0))**3D0 - dp0**3D0))-1D0
  
    ! CCN=SUM(PACK(N_bins, S_c <= S_super)) ! Number of CCN (activated droplets at 0.2 % supersaturation)
    ! m_drops=LWC/CCN ! Average single droplet mass/volume of cloud droplets (kg or dm^3)
	! r_drops=(m_drops/qH2O*3D0/(4D0*pi))**(1D0/3D0) ! Radius of cloud droplets
    ! L_cloud_vol_frac=LWC/qH2O ! Liquid water volume fraction of total air parcel
    
	! ! saturation vapor pressure of SO2 and H2O2 above the droplet surfaces
    ! !Ks1=1.3D-2*EXP(1960D0*(1D0/T-1D0/298D0)) ! mol L^-1 water  
    ! !Ks2=6.6D-8*EXP(1500D0*(1D0/T-1D0/298D0)) ! mol  L^-1 water 	
    ! !H_SO2=1.3D0*EXP(2900D0*(1D0/T-1D0/298D0))  ! Henry's law coefficient SO2(g) mol  L^-1 water atm^-1
    ! !H_S=H_SO2*(1D0+Ks1/cH+Ks1*Ks2/(cH**2D0))   ! effective Henry's law constant SO2(g)
	! !H_H2O2=8.4D4*EXP(7600D0*(1D0/T-1D0/298D0)) ! Henry's law coefficient H2O2(g) mol  L^-1 water atm^-1 at 298 K.
	! ! H_DMS  = 0.56*EXP(4480D0*(1D0/T-1D0/298D0)); ! Henry's law coefficient DMS(g) mol  L^-1 water atm^-1
    ! ! H_DMSO = 1D7*EXP(2580D0*(1D0/T-1D0/298D0));  ! Henry's law coefficient DMSO(g) mol  L^-1 water atm^-1
    ! ! H_MSIA  = 1D8*EXP(1760D0*(1D0/T-1D0/298D0)); ! Henry's law coefficient MSIOA(g) mol  L^-1 water atm^-1
    ! ! H_O3 = 0.011*EXP(2400D0*(1D0/T-1D0/298D0));  ! Henry's law coefficient O3(g) mol  L^-1 water atm^-1
	! ! H_OH = 3D1*EXP(4500D0*(1D0/T-1D0/298D0)); ! Henry's law coefficient OH(g) mol  L^-1 water atm^-1
    	 
	! pOH=OH*1D6/Na*Rg*T/1.01325D5 ! Initial OH partial pressure atm
    ! pO3=O3*1D6/Na*Rg*T/1.01325D5 ! Initial OH partial pressure atm
    ! pH2O2=H2O2*1D6/Na*Rg*T/1.01325D5 ! H2O2 partial pressure atm
	! pSO2=SO2*1D6/Na*Rg*T/1.01325D5 ! SO2 partial pressure atm
	! pDMS=DMS*1D6/Na*Rg*T/1.01325D5 ! DMS partial pressure atm
	! pDMSO=DMSO*1D6/Na*Rg*T/1.01325D5 ! DMSO partial pressure atm
	! pMSIA=MSIA*1D6/Na*Rg*T/1.01325D5 ! MSIA partial pressure atm
    ! pHPMTF=HPMTF*1D6/Na*Rg*T/1.01325D5 ! HPMTF partial pressure atm
	
    ! M_DMSaq=cDMSaq/(LWC*Na*1D-6) ! mol/(L cloud water)
    ! M_DMSOaq=cDMSOaq/(LWC*Na*1D-6) ! mol/(L cloud water)
    ! M_MSIAaq=cMSIAaq/(LWC*Na*1D-6) ! mol/(L cloud water)
    ! M_OHaq=cOHaq/(LWC*Na*1D-6) ! mol/(L cloud water)
    ! M_O3aq=cO3aq/(LWC*Na*1D-6) ! mol/(L cloud water)
    ! M_HPMTFaq=cHPMTFaq/(LWC*Na*1D-6) ! mol/(L cloud water)
    ! M_SO2aq=cSO2aq/(LWC*Na*1D-6) ! mol/(L cloud water)
    ! M_H2O2aq=cH2O2aq/(LWC*Na*1D-6) ! mol/(L cloud water)

    ! conc_aq(1)=pOH
	! conc_aq(2)=M_OHaq
    ! conc_aq(3)=pO3
	! conc_aq(4)=M_O3aq
    ! conc_aq(5)=pDMS
	! conc_aq(6)=M_DMSaq
    ! conc_aq(7)=pDMSO
	! conc_aq(8)=M_DMSOaq
    ! conc_aq(9)=pMSIA
	! conc_aq(10)=M_MSIAaq
	! conc_aq(11)=pHPMTF
	! conc_aq(12)=M_HPMTFaq
	! conc_aq(13)=pH2O2
	! conc_aq(14)=M_H2O2aq
	! conc_aq(15)=pSO2
	! conc_aq(16)=M_SO2aq
    ! conc_aq(17)=0D0 ! MSA production mol/L
	! conc_aq(18)=0D0 ! S(VI) production mol/L
	

! CALL cloud_aq(0D0,dt,conc_aq,T,L_cloud_vol_frac,r_drops)

	
! OH=conc_aq(1)*1.01325D5*Na/(1D6*Rg*T)
! M_OHaq=conc_aq(2)
! O3=conc_aq(3)*1.01325D5*Na/(1D6*Rg*T)
! M_O3aq=conc_aq(4)
! DMS=conc_aq(5)*1.01325D5*Na/(1D6*Rg*T)
! M_DMSaq=conc_aq(6)
! DMSO=conc_aq(7)*1.01325D5*Na/(1D6*Rg*T)
! M_DMSOaq=conc_aq(8)
! MSIA=conc_aq(9)*1.01325D5*Na/(1D6*Rg*T)
! M_MSIAaq=conc_aq(10)
! HPMTF=conc_aq(11)*1.01325D5*Na/(1D6*Rg*T)
! M_HPMTFaq=conc_aq(12)
! H2O2=conc_aq(13)*1.01325D5*Na/(1D6*Rg*T)
! M_H2O2aq=conc_aq(14)
! SO2=conc_aq(15)*1.01325D5*Na/(1D6*Rg*T)
! M_SO2aq=conc_aq(16)

! I_nCH3SO3H=m_drops*conc_aq(17)/dt ! CH3SO3H aqueous phase formation rate  mol/s per droplet
! I_nSVI=m_drops*conc_aq(18)/dt ! S(VI) formation rate mol/s

! cDMSaq=M_DMSaq*LWC*Na*1D-6 ! molec/(cm^3 air)
! cDMSOaq=M_DMSOaq*LWC*Na*1D-6 ! molec/(cm^3 air)
! cMSIAaq=M_MSIAaq*LWC*Na*1D-6 ! molec/(cm^3 air)
! cOHaq=M_OHaq*LWC*Na*1D-6 ! molec/(cm^3 air)
! cO3aq=M_O3aq*LWC*Na*1D-6 ! molec/(cm^3 air)
! cHPMTFaq=M_HPMTFaq*LWC*Na*1D-6 ! molec/(cm^3 air)
! cSO2aq=M_SO2aq*LWC*Na*1D-6 ! molec/(cm^3 air)
! cH2O2aq=M_H2O2aq*LWC*Na*1D-6 ! molec/(cm^3 air)


 	
! END SUBROUTINE cloud_activation_processing
    
! SUBROUTINE particle_halog_chem(c_p,N_bins,HOCL,CL2,HOBr,BrCl,HOI,ICl,&
	! cHOCLaq,cCL2aq,cHOBraq,cBrClaq,cHOIaq,cIClaq,CS_HOCL,CS_HOBr,CS_HOI,&
	! CS_CL2,CS_BrCl,CS_ICL,pH,T,dt)
    ! REAL(dp), DIMENSION(Nlayers,NSPEC_P,nr_bins), INTENT(in) :: c_p
    ! REAL(dp), DIMENSION(nr_bins), INTENT(in) :: N_bins,pH
    ! REAL(dp), INTENT(in) :: T,dt
    ! REAL(dp), INTENT(inout) :: HOCL,CL2,HOBr,BrCl,HOI,ICl
	! REAL(dp), INTENT(inout) :: cHOCLaq,cCL2aq,cHOBraq,cBrClaq,cHOIaq,cIClaq
	! REAL(dp), INTENT(in)    :: CS_HOCL,CS_HOBr,CS_HOI,CS_CL2,CS_BrCl,CS_ICL
	
	! REAL(dp)                :: pHOCL,pCL2,pHOBr,pBrCl,pHOI,pICl
	! REAL(dp)                :: M_HOCLaq,M_CL2aq,M_HOBraq,M_BrClaq,M_HOIaq,M_IClaq
	! REAL(dp)                :: LWCp,CLp,M_CLaq,M_Haq,p_aq_frac
	! REAL(dp), DIMENSION(NSPEC_AQ_HALOGEN) :: conc_halog_aq
	
    
! LWCp=SUM(c_p(Nlayers,7,:))/Na*MH2O ! Inorganic bulk liquid water content	kg/cm^3 = L/cm^3
! CLp=SUM(c_p(Nlayers,3,:))/Na ! Bulk Cl content mol/cm^3
! M_CLaq=CLp/(LWCp) ! mol/(L cloud water)

! M_Haq=SUM(c_p(Nlayers,7,:)/Na*MH2O*10**(-pH))/LWCp ! Bulk H+ content mol/L water

! p_aq_frac=LWCp/(qH2O*1D-6) ! Liquid particle water volume fraction of total air parcel

! !write(*,*) p_aq_frac

! ! (kg/cm^3)/(kg/L)= L/cm^3 =>
    ! pHOCl=HOCL*1D6/Na*Rg*T/1.01325D5 ! Initial HOCl partial pressure atm
    ! pCL2=CL2*1D6/Na*Rg*T/1.01325D5 ! Initial Cl2 partial pressure atm
    ! pHOBr=HOBr*1D6/Na*Rg*T/1.01325D5 ! HOBr partial pressure atm
	! pBrCl=BrCl*1D6/Na*Rg*T/1.01325D5 ! BrCl partial pressure atm
	! pHOI=HOI*1D6/Na*Rg*T/1.01325D5 ! HOI partial pressure atm
	! pICl=ICl*1D6/Na*Rg*T/1.01325D5 ! ICl partial pressure atm
	
    ! M_HOCLaq=cHOCLaq/(LWCp*Na) ! mol/(L cloud water)
    ! M_CL2aq=cCL2aq/(LWCp*Na) ! mol/(L cloud water)
    ! M_HOBraq=cHOBraq/(LWCp*Na) ! mol/(L cloud water)
    ! M_BrClaq=cBrClaq/(LWCp*Na) ! mol/(L cloud water)
    ! M_HOIaq=cHOIaq/(LWCp*Na) ! mol/(L cloud water)
    ! M_ICLaq=cICLaq/(LWCp*Na) ! mol/(L cloud water)
    
	
    ! conc_halog_aq(1)=pHOCL
	! conc_halog_aq(2)=M_HOCLaq
    ! conc_halog_aq(3)=pCL2
	! conc_halog_aq(4)=M_CL2aq
    ! conc_halog_aq(5)=pHOBr
	! conc_halog_aq(6)=M_HOBraq
    ! conc_halog_aq(7)=pBrCl
	! conc_halog_aq(8)=M_BrClaq
    ! conc_halog_aq(9)=pHOI
	! conc_halog_aq(10)=M_HOIaq
	! conc_halog_aq(11)=pICL
	! conc_halog_aq(12)=M_ICLaq
	
    	
! CALL halogen_aq(0D0,dt,conc_halog_aq,T,M_Haq,M_CLaq,CS_HOCL,CS_HOBr,CS_HOI,CS_CL2,CS_BrCl,CS_ICL,p_aq_frac)


! HOCl=conc_halog_aq(1)*1.01325D5*Na/(1D6*Rg*T)
! Cl2=conc_halog_aq(3)*1.01325D5*Na/(1D6*Rg*T)
! HOBr=conc_halog_aq(5)*1.01325D5*Na/(1D6*Rg*T)
! BrCl=conc_halog_aq(7)*1.01325D5*Na/(1D6*Rg*T)
! HOI=conc_halog_aq(9)*1.01325D5*Na/(1D6*Rg*T)
! ICL=conc_halog_aq(11)*1.01325D5*Na/(1D6*Rg*T)

! M_HOCLaq=conc_halog_aq(2)
! M_CL2aq=conc_halog_aq(4)
! M_HOBraq=conc_halog_aq(6)
! M_BrClaq=conc_halog_aq(8)
! M_HOIaq=conc_halog_aq(10)
! M_ICLaq=conc_halog_aq(12)
	
! cHOCLaq=M_HOCLaq*(LWCp*Na) ! mol/(L cloud water)
! cCL2aq=M_CL2aq*(LWCp*Na) ! mol/(L cloud water)
! cHOBraq=M_HOBraq*(LWCp*Na) ! mol/(L cloud water)
! cBrClaq=M_BrClaq*(LWCp*Na) ! mol/(L cloud water)
! cHOIaq=M_HOIaq*(LWCp*Na) ! mol/(L cloud water)
! cICLaq=M_ICLaq*(LWCp*Na) ! mol/(L cloud water)
    	
! END SUBROUTINE particle_halog_chem 
    
    !--------------------------------------------------------------! 
    ! Calculate fraction of POA in each sizebin                    !
    !--------------------------------------------------------------!
    SUBROUTINE fraction_POA_marine(d_p,f_POA_marine) 
 
      REAL(dp), DIMENSION(nr_bins), INTENT(in) :: d_p
      REAL(dp), DIMENSION(nr_bins), INTENT(out)  :: f_POA_marine
      INTEGER  :: j
 
      ! fraction of OA using cubic fit of data-points from O'Dowd et al 2004
      f_POA_marine = -6.178D-10*(d_p*1D9)**3 + 1.599D-6*(d_p*1D9)**2 - 0.001681D0*d_p*1D9 + 1.004D0 
      DO j=1,nr_bins
 IF (f_POA_marine(j) > 0.99) THEN
   f_POA_marine(j) = 0.99
 ELSE IF (f_POA_marine(j) < 0.01) THEN
   f_POA_marine(j) = 0.01
 END IF 
      END DO
 
    END SUBROUTINE fraction_POA_marine

SUBROUTINE sea_salt_ice(U10m,temp,dp_dry,u,SHTF,E_snow_flux)
REAL(dp), INTENT(in) :: U10m,temp,u,SHTF
REAL(dp), DIMENSION(nr_bins), INTENT(in) :: dp_dry
REAL(dp), DIMENSION(nr_bins) :: vp_ssa, vp_snow,Dp_snow,f_snow,dM_snow_flux,dV_snow_flux
REAL(dp) :: A_age, Ut0, Ut, L_heat_sub, Qs_prim, Qs, temp_c,dens_NaCl,dens_ice,salt_snow,&
alpha,Beta,Tau,qbsalt,qb0      
REAL(dp), DIMENSION(nr_bins), INTENT(out) :: E_snow_flux

A_age=1D0 ! Empirical snow age factor
Ut0=6.975 ! m/s
L_heat_sub=2.838D6 ! J/kg Latent heat of sublimation Stigter et al. (2018)
Qs_prim=SHTF/L_heat_sub ! kg m^-2 s^-1 (Estimated snow sublimation rate)
temp_c=temp-273.15
dens_NaCl=2160D0 ! kg/m^3
dens_ice=900D0 ! kg/m^3
salt_snow=1D-5 ! Estimated snow salt content in mass %

IF (Qs_prim>0D0) THEN
vp_ssa=(pi*dp_dry**3D0)/6D0
vp_snow=(vp_ssa*dens_NaCl/salt_snow)/dens_ice 
Dp_snow=1D6*(vp_snow*6D0/pi)**(1D0/3D0) !um
alpha=3.2D0
Beta=131D0/alpha
Tau=GAMMA(alpha)
f_snow=(EXP(-Dp_snow/Beta)*Dp_snow**(alpha-1D0))/(Tau*Beta**alpha) ! Blowing snow spectrum 

Ut=Ut0+0.0033*(temp_c+27.27)**2D0

IF (U10m>Ut) THEN
qbsalt=(0.385*(1D0-Ut/U10m)**2.59D0)/u ! Salt layer blowing snow mixing rate kg/kg
ELSE
qbsalt=0D0
END IF

IF (U10m>Ut0) THEN
qb0=(0.385*(1D0-Ut0/U10m)**2.59D0)/u ! Salt layer blowing snow mixing rate ref kg/kg
ELSE
qb0=1D0
END IF

Qs=Qs_prim*A_age*qbsalt/qb0 ! Blowing snow bulk sublimation flux (kg m^-2 s^-1)

dM_snow_flux=Qs*f_snow !
dV_snow_flux=dM_snow_flux/dens_ice ! m^3 m^-2 s^-1
E_snow_flux=dV_snow_flux/vp_snow ! # m^-2 s^-1


ELSE
E_snow_flux=0D0
END IF

END SUBROUTINE sea_salt_ice  

    !--------------------------------------------------------------!
    ! Sea-spray particle emission                                  !
    !--------------------------------------------------------------!
SUBROUTINE sea_spray(wind10,temp,dlogdp,dp_dry,d_g,d,E_sea_salt,sea_spray_index)
      REAL(dp), INTENT(in) :: wind10,temp
      REAL(dp), DIMENSION(nr_bins), INTENT(in) :: dlogdp,dp_dry,d_g
	  REAL(dp), DIMENSION(nr_bins+1), INTENT(in) :: d
      REAL(dp), DIMENSION(nr_bins), INTENT(out) :: E_sea_salt
	  REAL(dp), DIMENSION(nr_bins) :: PV_sea_salt

      REAL(dp) :: A_whitecape,T_water,cc4,cc3,cc2,cc1,cc0,dd4,dd3,dd2,dd1,dd0, &
             ccc4,ccc3,ccc2,ccc1,ccc0,ddd4,ddd3,ddd2,ddd1,ddd0, &
             cccc4,cccc3,cccc2,cccc1,cccc0,dddd4,dddd3,dddd2,dddd1,dddd0
      REAL(dp), DIMENSION(nr_bins) :: Ak,Bk,dF_saltdlogDp
	  REAL(dp), DIMENSION(nr_bins) :: V_sea_salt,SST_corr1,SST_corr2,SST_corr3,SST_corr4,SST_corr
      REAL(dp), DIMENSION(5) :: dm,sigma,N_modes
      REAL(dp) :: N1,N2,N3,F_ent,A1,A2,A3,B1,B2,B3,C1,C2,C3,D1,D2,D3,temp1
	  INTEGER, INTENT(in) :: sea_spray_index
	  

      INTEGER :: i,modes=5


IF (sea_spray_index==1) THEN
      A_whitecape = 0.01*3.84D-4*wind10**3.41 ! estimated whitecap area fraction
      ! Estimate sea surface temperature based on temperature at 2 m
         T_water = temp

      DO i = 1,nr_bins
 IF (dp_dry(i) < 20D-9) THEN 
   cc4 = -2.576D35
   cc3 = 5.932D28
   cc2 = -2.867D21
   cc1 = -3.003D13
   cc0 = -2.881D6
   dd4 = 7.188D37
   dd3 = -1.616D31
   dd2 = 6.791D23
   dd1 = 1.829D16
   dd0 = 7.609D8   

   Ak(i) = (dp_dry(i)/20D-9)*(cc4*(20D-9)**4.+cc3*(20D-9)**3.+cc2*(20D-9)**2.+cc1*(20D-9)+cc0)
   Bk(i) = (dp_dry(i)/20D-9)*(dd4*(20D-9)**4.+dd3*(20D-9)**3.+dd2*(20D-9)**2.+dd1*(20D-9)+dd0)

 ELSE IF (dp_dry(i) < 145D-9) THEN
   cc4 = -2.576D35
   cc3 = 5.932D28
   cc2 = -2.867D21
   cc1 = -3.003D13
   cc0 = -2.881D6
    dd4 = 7.188D37
   dd3 = -1.616D31
   dd2 = 6.791D23
   dd1 = 1.829D16
   dd0 = 7.609D8

   Ak(i) = cc4*(dp_dry(i))**4.+cc3*(dp_dry(i))**3.+cc2*(dp_dry(i))**2.+cc1*(dp_dry(i))+cc0
   Bk(i) = dd4*(dp_dry(i))**4.+dd3*(dp_dry(i))**3.+dd2*(dp_dry(i))**2.+dd1*(dp_dry(i))+dd0  

 ELSE IF (dp_dry(i) < 419D-9) THEN
   ccc4 = -2.452D33
   ccc3 = 2.404D27
   ccc2 = -8.148D20
   ccc1 = 1.183D14
   ccc0 = -6.743D6
   ddd4 = 7.368D35
   ddd3 = -7.310D29
   ddd2 = 2.528D23
   ddd1 = -3.787D16
   ddd0 = 2.279D9

   Ak(i) = ccc4*(dp_dry(i))**4.+ccc3*(dp_dry(i))**3.+ccc2*(dp_dry(i))**2.+ccc1*(dp_dry(i))+ccc0
   Bk(i) = ddd4*(dp_dry(i))**4.+ddd3*(dp_dry(i))**3.+ddd2*(dp_dry(i))**2.+ddd1*(dp_dry(i))+ddd0

 ELSE IF (dp_dry(i) < 2800d-9) THEN
   cccc4 = 1.085D29
   cccc3 = -9.841D23
   cccc2 = 3.132D18
   cccc1 = -4.165D12
   cccc0 = 2.181D6
   dddd4 = -2.859D31
   dddd3 = 2.601D26
   dddd2 = -8.297D20
   dddd1 = 1.105D15
   dddd0 = -5.800D8

   Ak(i) = cccc4*(dp_dry(i))**4.+cccc3*(dp_dry(i))**3.+cccc2*(dp_dry(i))**2.+cccc1*(dp_dry(i))+cccc0
   Bk(i) = dddd4*(dp_dry(i))**4.+dddd3*(dp_dry(i))**3.+dddd2*(dp_dry(i))**2.+dddd1*(dp_dry(i))+dddd0

 ELSE 
   cccc4 = 1.085D29
   cccc3 = -9.841D23
   cccc2 = 3.132D18
   cccc1 = -4.165D12
   cccc0 = 2.181D6
   dddd4 = -2.859D31
   dddd3 = 2.601D26
   dddd2 = -8.297D20
   dddd1 = 1.105D15
   dddd0 = -5.800D8

   Ak(i) = (2800D-9/dp_dry(i))**20.*(cccc4*(2800D-9)**4.+cccc3*(2800D-9)**3.+cccc2*(2800D-9)**2.+cccc1*(2800D-9)+cccc0)
   Bk(i) = (2800D-9/dp_dry(i))**20.*(dddd4*(2800D-9)**4.+dddd3*(2800D-9)**3.+dddd2*(2800D-9)**2.+dddd1*(2800D-9)+dddd0)

 END IF 
   
 dF_saltdlogDp(i) = MAXVAL((/ A_whitecape*(Ak(i)*T_water+Bk(i)),0D0 /)) ! m^-2 s^-1
      END DO
      E_sea_salt = dF_saltdlogDp*dlogdp ! m^-2 s^-1 emission of sea salt particles
ELSEIF (sea_spray_index==2) THEN

temp1=temp
IF (temp1<271.15) THEN
temp1=271.15
END IF

! Alt. Salter et al., ACP 
F_ent=2D-8*wind10**3.41D0

A1=-5.2168D5;A2=0D0;A3=0D0
B1=3.31725D7;B2=7.374D5;B3=1.421D4
C1=-6.95275D8;C2=-2.4803D7;C3=1.4662D7
D1=1.0684D10;D2=7.7373D8;D3=1.7075D8

N1=F_ent*(A1*(temp1-273.15)**3D0+B1*(temp1-273.15)**2D0+C1*(temp1-273.15)+D1)
N2=F_ent*(A2*(temp1-273.15)**3D0+B2*(temp1-273.15)**2D0+C2*(temp1-273.15)+D2)
N3=F_ent*(A3*(temp1-273.15)**3D0+B3*(temp1-273.15)**2D0+C3*(temp1-273.15)+D3)

dm = (/ 95D-9, 6D-7, 1.5D-6, 2D-6, 2D-6 /) ! Mode diameter size distribution (m)
sigma = (/ 2.1, 1.72, 1.6, 1.5, 1.5 /)  ! Standard deviation for each mode
N_modes = (/ N1, N2, N3, 0D0, 0D0 /) ! Number concentration modes (#/m^3)

E_sea_salt=&
N1/(sqrt(2D0*pi)*log10(sigma(1)))*exp(-0.5D0*((log10(d_g)-log10(dm(1)))**2)/log10(sigma(1))**2)+&
N2/(sqrt(2D0*pi)*log10(sigma(2)))*exp(-0.5D0*((log10(d_g)-log10(dm(2)))**2)/log10(sigma(2))**2)+&
N3/(sqrt(2D0*pi)*log10(sigma(3)))*exp(-0.5D0*((log10(d_g)-log10(dm(3)))**2)/log10(sigma(3))**2)
E_sea_salt=E_sea_salt*dlogDp

!CALL dNdlogDp(d_g,dp_dry,dlogDp,modes,sigma,N_modes,dm,E_sea_salt,PV_sea_salt)


ELSEIF (sea_spray_index==3) THEN ! Sofiev et al. (2011)
temp1=temp
IF (temp1<271.15) THEN
temp1=271.15
END IF

! Sofiev et al. (2011)
E_sea_salt=(3.84D0*wind10**3.41D0)*&
(EXP(-9D-2/(d_g*1D6+3D-3))/(2D0+EXP(-5D0/(d_g*1D6))))*&
((1D0+5D-2*(d_g*1D6)**1.05D0)/((d_g*1D6)**3))*&
1D1**(1.05D0*EXP(-((0.27D0-LOG10(d_g*1D6))/1.1D0)**2))

! E_sea_salt2=(3.84*wind10.^3.41).*...
! (exp(-0.09./(dp_dry*1D6+0.003))./(2D0+exp(-5D0./(dp_dry*1D6)))).*...
! ((1D0+0.05*(dp_dry*1D6).^1.05)./((dp_dry*1D6).^3)).*...
! 1D1.^(1.05*exp(-((0.27-log10(dp_dry*1D6))/1.1).^2));



A1=0.092D0
A2=0.15D0
A3=0.48D0
B1=-0.96D0
B2=-0.88D0
B3=-0.36D0
SST_corr1=A1*(dp_dry*1D6)**B1
SST_corr2=A2*(dp_dry*1D6)**B2
SST_corr3=A3*(dp_dry*1D6)**B3
SST_corr4=1D0
IF (temp1<=271.15) THEN
SST_corr=SST_corr1
ELSEIF (temp1<=278.15) THEN
SST_corr=SST_corr1*(278.15-temp1)/7D0 + SST_corr2*(temp1-271.15)/7D0
ELSEIF (temp1<=288.15) THEN
SST_corr=SST_corr2*(288.15-temp1)/1D1 + SST_corr3*(temp1-278.15)/1D1
ELSEIF (temp1<=298.15) THEN
SST_corr=SST_corr3*(298.15-temp1)/1D1 + SST_corr4*(temp1-288.15)/1D1
ELSE
SST_corr=SST_corr4
END IF

E_sea_salt=E_sea_salt*SST_corr*(d(2:nr_bins+1)*1D6-d(1:nr_bins)*1D6);

END IF

END SUBROUTINE sea_spray                   




SUBROUTINE coagulation_clusterin(N_bins,V_bins,c_p,c_p_backg,d_p,dp_dry,dt,dens_p,MX,qX,chem_1,chem_2,chem_3,chem_4,clustering_vapors,three_comp_sys, &
  log_array, ncases)

  REAL(dp), DIMENSION(nr_bins), INTENT(inout) :: N_bins
  REAL(dp), DIMENSION(NSPEC_P,nr_bins), INTENT(inout) :: c_p
  REAL(dp), DIMENSION(nr_bins), INTENT(inout) :: d_p, dp_dry
  REAL(dp), DIMENSION(nr_bins), INTENT(inout) :: dens_p
  REAL(dp), DIMENSION(NSPEC_P), INTENT(in) :: MX,qX
  REAL(dp), INTENT(in) :: dt
  REAL(dp), DIMENSION(nr_bins), INTENT(out) :: V_bins
  REAL(dp), DIMENSION(NSPEC_P,nr_bins), INTENT(in) :: c_p_backg 

  REAL(dp) :: a,Vp_coag,r1,r2,Coag_source, Coag_sinktot
  REAL(dp), DIMENSION(nr_bins) :: vp_dry,vp_wet
  
  REAL(dp), DIMENSION(nr_bins+1) :: Vp
  REAL(dp), DIMENSION(NSPEC_P,nr_bins) :: c_p_single
 
  
  INTEGER :: i,j,m,n, ic,ff

  !! for cluster 
  type(clustering_mod), intent(in), optional:: chem_1, chem_2,chem_3, chem_4
  integer,intent(in) :: clustering_vapors
  INTEGER, allocatable :: nclust(:), list_nclust(:)
  INTEGER:: clustering_systems                            ! Number of clusters
  REAL(dp), allocatable :: dN_coag_clust(:,:), c_p_single_clust(:,:), Vp_clust(:) ! Concentrations, compositions and volumes of coagulated clusters
  REAL(dp), DIMENSION(nr_bins+1) :: dN                     ! Changes in aerosol number and composition; note: dN, not dNdt as the changes are already integrated by
  REAL(dp), DIMENSION(NSPEC_P,nr_bins+1) :: dcp            ! the cluster dynamics routine
  REAL(dp) :: N_bins_old(nr_bins), N_fix, N_fix_tmp,dp_max
  REAL(DP) ::  c_p_coag_cluster(NSPEC_P)
  logical, intent(in):: three_comp_sys
  logical, intent(in):: log_array(4)
  integer, intent(in) :: ncases 
  

  dp_max = dp_dry(nr_bins)*dp_dry(nr_bins)/dp_dry(nr_bins-1)
  Vp(1:nr_bins) = (4D0*pi*((dp_dry/2D0))**3D0)/3D0  ! Dry singel particle volume in each size bin m^3.
  Vp(nr_bins+1) = (pi*dp_max**3.)/6.

  !!!!!!!!! CLUSTER COAGULATION
  !!!! no self coagulation
    dN = 0.
    dcp=0.
    ff=0

    ! c_p_single=c_p!sum(c_p_single,dim=1)
    DO j=1,nr_bins
      c_p_single(:,j) = c_p(:,j)/(N_bins(j)*1D-6) ! concentration of each condensable compound in particle phase in a single particle in each bin
    END DO
    
    ! allocate(list_nclust(clustering_systems))
    allocate(nclust(clustering_systems))
    
    ! list_nclust=(/chem_1%nclust_syst, chem_2%nclust_syst,chem_3%nclust_syst,chem_4%nclust_syst/)

    ! do i=1,size(log_array)
    !   if (log_array(i)) then
    !     ff=ff+1
    !     nclust(ff)=list_nclust(i)
    !   end if
    ! end do

    nclust=(/chem_1%nclust_syst, chem_2%nclust_syst,chem_3%nclust_syst,chem_4%nclust_syst/)
    ! c_p_single=c_p!sum(c_p_single,dim=1)
  
    Do ic=1,2 !1, clustering_systems
    ! Coagulation source:
      if (ic==1) then
        ! write(*,*) 'Here L 262, chem 1'
        allocate(dN_coag_clust(nr_bins,nclust(ic)), c_p_single_clust(NSPEC_P,nclust(ic)), Vp_clust(nclust(ic)))
        vp_clust=chem_1%v_clust
        dN_coag_clust=chem_1%conc_coag_clust
        c_p_single_clust=chem_1%c_p_clust
      elseif (ic==2) then
        allocate(dN_coag_clust(nr_bins,nclust(ic)), c_p_single_clust(NSPEC_P,nclust(ic)), Vp_clust(nclust(ic)))
        ! write(*,*) 'Here L 262, chem 2'
        vp_clust=chem_2%v_clust
        dN_coag_clust=chem_2%conc_coag_clust
        c_p_single_clust=chem_2%c_p_clust
      elseif (ic==3) then
        allocate(dN_coag_clust(nr_bins,nclust(ic)), c_p_single_clust(NSPEC_P,nclust(ic)), Vp_clust(nclust(ic)))
        ! write(*,*) 'Here L 262, chem 3'
        vp_clust=chem_3%v_clust
        dN_coag_clust=chem_3%conc_coag_clust
        c_p_single_clust=chem_3%c_p_clust
      elseif (ic==4) then
        allocate(dN_coag_clust(nr_bins,nclust(ic)), c_p_single_clust(NSPEC_P,nclust(ic)), Vp_clust(nclust(ic)))
        ! write(*,*) 'Here L 262, chem 4'
        vp_clust=chem_4%v_clust
        dN_coag_clust=chem_4%conc_coag_clust
        c_p_single_clust=chem_4%c_p_clust
        
      end if

    
    ! Do ic=1, size(nclust)!clustering_systems
    ! ! Coagulation source:
    !   if (ic==1) then
    !     ! write(*,*) 'Here L 262, chem 1'
    !     ! allocate(dN_coag_clust(nr_bins,nclust(ic)), c_p_single_clust(NSPEC_P,nclust(ic)), Vp_clust(nclust(ic)))
    !     ! vp_clust=chem_1%v_clust
    !     ! dN_coag_clust=chem_1%conc_coag_clust
    !     ! c_p_single_clust=chem_1%c_p_clust
    !     ! write(*,*) sum(Vp_clust),SUM(dN_coag_clust),SUM(c_p_single_clust)
    !   ! elseif (ic==2) then
    !     ! write(*,*) 'Here L 262, chem 2'
    !     ! allocate(dN_coag_clust(nr_bins,nclust(ic)), c_p_single_clust(NSPEC_P,nclust(ic)), Vp_clust(nclust(ic)))
    !     ! vp_clust=chem_2%v_clust
    !     ! dN_coag_clust=chem_2%conc_coag_clust
    !     ! c_p_single_clust=chem_2%c_p_clust
    !     ! write(*,*) 'l927'
    !     ! write(*,*) sum(Vp_clust),SUM(dN_coag_clust),SUM(c_p_single_clust)
    !     ! elseif (ic==3) then
    !       ! write(*,*) 'skip HIO3-HIO2 system'
    !     allocate(dN_coag_clust(nr_bins,nclust(ic)), c_p_single_clust(NSPEC_P,nclust(ic)), Vp_clust(nclust(ic)))
    !     write(*,*) 'Here L 262, chem 3'
    !     vp_clust=chem_3%v_clust
    !     dN_coag_clust=chem_3%conc_coag_clust
    !     ! c_p_single_clust=chem_3%c_p_clust
    !     ! write(*,*) sum(Vp_clust),SUM(dN_coag_clust),SUM(c_p_single_clust)
    !     ! elseif (ic==4) then
    !       ! write(*,*) 'Here L 262, chem 4'
    !       !   allocate(dN_coag_clust(nr_bins,nclust(ic)), c_p_single_clust(NSPEC_P,nclust(ic)), Vp_clust(nclust(ic)))
    !       !   vp_clust=chem_4%v_clust
    !       !   dN_coag_clust=chem_4%conc_coag_clust
    !       !   c_p_single_clust=chem_4%c_p_clust  
    !       !   write(*,*) sum(Vp_clust),SUM(dN_coag_clust),SUM(c_p_single_clust)
    !   end if
    
      DO j = 1,nr_bins
        DO m = 1,nclust(ic)!nclust
          Vp_coag = vp_clust(m) +Vp(j) !Vp_clust(m)
          c_p_coag_cluster = c_p_single_clust(:,m)+c_p_single(:,j)!c_p_single_clust(:,m)+c_p_single(:,j)
          DO  i = j,nr_bins
            IF (Vp_coag >= Vp(i) .and. Vp_coag < Vp(i+1)) THEN
              r1 = (Vp(i+1)-Vp_coag)/(Vp(i+1)-Vp(i))
              r2 = 1.-r1
              Coag_source = dN_coag_clust(j,m)!dN_coag_clust(j,m) ! (#/m^3)
              dN(i) = dN(i)+r1*Coag_source
              dN(i+1) = dN(i+1)+r2*Coag_source
              dcp(:,i) = dcp(:,i)+Vp(i)/Vp_coag*r1*Coag_source*c_p_coag_cluster ! (molec/(m^3))
              dcp(:,i+1) = dcp(:,i+1)+Vp(i+1)/Vp_coag*r2*Coag_source*c_p_coag_cluster
            END IF
          END DO
        END DO
      END DO

      DO j = 1,nr_bins
        Coag_sinktot = sum(dN_coag_clust(j,:))!dN_coag_clust(j,:)) ! (#/m^3)
        dN(j) = dN(j)-Coag_sinktot
        dcp(:,j) = dcp(:,j)-Coag_sinktot*c_p_single(:,j) ! (molec/m^3)
      END DO

      N_bins = N_bins+dN(1:nr_bins) ! New particle concentration in each size bin
      c_p(:,:) = c_p(:,:)+dcp(:,1:nr_bins)*1D-6 ! molec/ cm^3

      !! reset for every ic
      dN=0D0
      dcp=0D0

      deallocate(dN_coag_clust, vp_clust,c_p_single_clust)

      DO j = 1,nr_bins
        IF (N_bins(j)<1D-6) THEN
          N_bins(j)=1D-3		
          c_p(:,j)=c_p_backg(:,j)*1D-3
          ! c_p(:,j)=c_p(:,j)*Vp(j)/SUM(c_p(index_dry,j)/Na*MX(index_dry)/qX(index_dry))*N_bins(j)*1D-6
        end if
      end do

    End do

    ! write(*,*) 'here'

  DO j = 1,nr_bins
    vp_dry(j)=SUM(c_p(index_dry,j)/Na*MX(index_dry)/qX(index_dry)/(N_bins(j)*1D-6)) ! m^3   
    vp_wet(j)=SUM(c_p(:,j)/Na*MX/qX/(N_bins(j)*1D-6)) ! m^3
    V_bins(j)=N_bins(j)*vp_wet(j)
    dens_p(j)=SUM(c_p(:,j)*MX*1D6)/Na/V_bins(j) ! Total particle density
  END DO
  dp_dry=(vp_dry*6D0/pi)**(1D0/3D0) ! Dry particle diameters
  d_p=(vp_wet*6D0/pi)**(1D0/3D0) ! Dry particle diameters

END SUBROUTINE coagulation_clusterin


    !------------------------------------------------------!
    ! Coagulation                                          !
    !------------------------------------------------------!


!!!! original coagulation

SUBROUTINE coagulation(N_bins,V_bins,c_p,d_p,dp_dry,dt,T,p,dens_p,MX,qX)
   
  REAL(dp), DIMENSION(nr_bins), INTENT(inout) :: N_bins
  REAL(dp), DIMENSION(NSPEC_P,nr_bins), INTENT(inout) :: c_p
  REAL(dp), DIMENSION(nr_bins), INTENT(inout) :: d_p, dp_dry
  REAL(dp), DIMENSION(nr_bins), INTENT(inout) :: dens_p
  REAL(dp), DIMENSION(NSPEC_P), INTENT(in) :: MX,qX
  REAL(dp), INTENT(in) :: dt,T,p
  REAL(dp), DIMENSION(nr_bins), INTENT(out) :: V_bins

  REAL(dp) :: dyn_visc,dp_max,a,Vp_coag,r1,r2,Coag_source, Coag_sinktot, l_gas
  REAL(dp), DIMENSION(nr_bins) :: C,D,m_p,speed_p,free_path_p,dist, beta,coag_sink,&
                                  vp_dry,vp_wet
  REAL(dp), DIMENSION(nr_bins,nr_bins) :: K
  REAL(dp), DIMENSION(nr_bins+1) :: Vp,dNdt
  REAL(dp), DIMENSION(NSPEC_P,nr_bins) :: c_p_single
  REAL(dp), DIMENSION(NSPEC_P,nr_bins+1) :: dcpdt
  REAL(dp), DIMENSION(NSPEC_P) :: c_p_coag
  
  INTEGER :: i,j,m,n

  dyn_visc = 1.8D-5*(T/298.)**0.85  ! dynamic viscosity of air
  l_gas=2D0*dyn_visc/(p*SQRT(8D0*Mair/(pi*Rg*T))) ! Gas mean free path in air (m)

  DO j=1,nr_bins
    c_p_single(:,j) = c_p(:,j)/(N_bins(j)*1D-6) ! concentration of each condensable compound in particle phase in a single particle in each bin
  END DO

  C = 1D0+(2D0*l_gas/(d_p))*(1.257+0.4*exp(-1.1/(2D0*l_gas/d_p))) ! Cunninghams correction factor (seinfeld and Pandis eq 9.34
  D = C*kb*T/(3D0*pi*dyn_visc*d_p)                  ! Diffusivitys for the different particle sizes m^2/s
  m_p = (dens_p*pi*(d_p)**3D0)/6D0                    ! mass of particles  (kg)
  speed_p = SQRT(8D0*kb*T/(pi*m_p))                 ! speed of particles (m/s)
  free_path_p = 8D0*D/(pi*speed_p)                   ! particle mean free path
  dist = (1./(3.*d_p*free_path_p))*((d_p+free_path_p)**3. &
-(d_p**2.+free_path_p**2.)**(3./2.))-d_p ! mean distance from the center of a sphere reached by particles leaving the sphere's surface

  DO n = 1,nr_bins
    beta = ((d_p+d_p(n))/(d_p+d_p(n)+2.*(dist**2.+dist(n)**2.)**0.5)+ &
8.*(D+D(n))/(((speed_p**2.+speed_p(n)**2.)**0.5)*(d_p+d_p(n))))**(-1.) ! Fuchs correction factor from Seinfeld and Pandis, 2006, p. 600
    K(n,:) = 2.*pi*beta*(d_p*D(n)+d_p*D+d_p(n)*D+d_p(n)*D(n))              ! coagulation rates between two particles of all size combinations  (m^3/s)    
  END DO

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  dp_max = dp_dry(nr_bins)*dp_dry(nr_bins)/dp_dry(nr_bins-1)
  Vp(1:nr_bins) = (4D0*pi*((dp_dry/2D0))**3D0)/3D0  ! Dry singel particle volume in each size bin m^3.
  Vp(nr_bins+1) = (pi*dp_max**3.)/6.

  dNdt = 0.
  dcpdt=0.
  ! Coagulation source:
  DO j = 1,nr_bins
    DO m = 1,j
      IF (m == j) THEN 
        a = 0.5 ! self-coagulation
      ELSE
        a = 1.
      END IF
      Vp_coag = Vp(m)+Vp(j)   ! Volume of new particles formed by coagulation
      c_p_coag = c_p_single(:,m)+c_p_single(:,j) ! compounds in formed single particle (molecues/#) 
      DO  i = j,nr_bins
        IF (Vp_coag >= Vp(i) .and. Vp_coag < Vp(i+1)) THEN
        ! If i=j some of the particles will stay in the same size bin as before. 
        ! This is anyway treated as a source into the old size bin because it is 
        ! assumed that all particles are lost due to coagulation in that same size 
        ! bin below
        ! Splitting parameters:
          r1 = (Vp(i+1)-Vp_coag)/(Vp(i+1)-Vp(i)) ! Fraction of particles in size bin i
          r2 = 1.-r1 ! Fraction of particles in next size bin (i+1)
          Coag_source = a*K(m,j)*N_bins(j)*N_bins(m) ! (m³/s)
          dNdt(i) = dNdt(i)+r1*Coag_source ! (#/m^3 s)
          dNdt(i+1) = dNdt(i+1)+r2*Coag_source ! (#/m^3 s)
          dcpdt(:,i) = dcpdt(:,i)+Vp(i)/Vp_coag*r1*Coag_source*c_p_coag ! (molec/(m^3 s))
          dcpdt(:,i+1) = dcpdt(:,i+1)+Vp(i+1)/Vp_coag*r2*Coag_source*c_p_coag ! (molec/(m^3 s))
        END IF
END DO
END DO
  END DO

  ! Coagulation sink included:
  DO j = 1,nr_bins
    Coag_sink = K(j,:)*N_bins*N_bins(j) ! (m^3/s)
    Coag_sinktot = sum(Coag_sink)
    dNdt(j) = dNdt(j)-Coag_sinktot
    dcpdt(:,j) = dcpdt(:,j)-Coag_sinktot*c_p_single(:,j) ! (molec/s)
  END DO


  N_bins = N_bins+dNdt(1:nr_bins)*dt ! New particle concentration in each size bin
  c_p = c_p+dcpdt(:,1:nr_bins)*dt*1D-6 ! molec/ cm^3

  DO j = 1,nr_bins
    vp_dry(j)=SUM(c_p(index_dry,j)/Na*MX(index_dry)/qX(index_dry)/(N_bins(j)*1D-6)) ! m^3   
    vp_wet(j)=SUM(c_p(:,j)/Na*MX/qX/(N_bins(j)*1D-6)) ! m^3
    V_bins(j)=N_bins(j)*vp_wet(j)
    dens_p(j)=SUM(c_p(:,j)*MX*1D6)/Na/V_bins(j) ! Total particle density
  END DO
  dp_dry=(vp_dry*6D0/pi)**(1D0/3D0) ! Dry particle diameters
  d_p=(vp_wet*6D0/pi)**(1D0/3D0) ! Dry particle diameters
END SUBROUTINE coagulation

    !--------------------------------------------------------------!
    ! Condensation                                                 !
    !--------------------------------------------------------------!
	
	!!!!!!!!!!!!!!!!!! Condensation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE condensation(N_bins,V_bins,p,T,RH,d_p,dp_dry,aX,MX,qX,dX,dens,corg,cH2SO4,cHNO3,&
cHCl,cCH3SO3H,cHIO3,cNH3,cDMA,cHIO2,c_p,W,Kprim_HNO3,Kprim_HCl,Kprim_CH3SO3H,Kprim_HIO3,Hprim_NH3,Kprim_NH3,&
fHSO4,fSO4,fNO3,fCl,fCH3SO3,fHIO3,mHCO3,mCO3,mOH,psat_org,yorg,Dorg,CS_H2SO4,CS_air,c_p_backg,dt,&
 Nconc_evap1, comp_evap, clust_evap, flag_clusterin)    

REAL(dp), DIMENSION(nr_bins), INTENT(in) :: W,Kprim_HNO3,Kprim_HCl,Kprim_CH3SO3H,Kprim_HIO3,Hprim_NH3,&
Kprim_NH3,fHSO4,fSO4,fNO3,fCl,fCH3SO3,fHIO3,mHCO3,mCO3,mOH
REAL(dp), DIMENSION(NSPEC_P,nr_bins), INTENT(in) :: c_p_backg
REAL(dp), DIMENSION(nr_bins), INTENT(inout) :: d_p, dp_dry, dens, N_bins
REAL(dp), DIMENSION(nr_bins), INTENT(out) :: V_bins
REAL(dp), DIMENSION(NSPEC_P,nr_bins), INTENT(inout) :: c_p
REAL(dp), DIMENSION(NSPEC_P,nr_bins) :: c_p_fixed, c_p_old
REAL(dp), INTENT(in) :: T,p,RH
REAL(dp), DIMENSION(NCOND,nr_bins), INTENT(in) :: yorg
REAL(dp), INTENT(out) :: CS_H2SO4, CS_air
REAL(dp), INTENT(inout) :: cH2SO4, cHNO3, cNH3, cHCl, cCH3SO3H, cHIO3, cDMA,cHIO2
REAL(dp), DIMENSION(NCOND), INTENT(inout) :: corg
REAL(dp), DIMENSION(NCOND), INTENT(in)    :: psat_org
REAL(dp), DIMENSION(NSPEC_P), INTENT(in) :: MX,qX,aX,dX
REAL(dp), DIMENSION(NSPEC_P) :: surf_tens
REAL(dp), DIMENSION(NSPEC_P,nr_bins) :: S_Kelvin
REAL(dp) :: aH2SO4,aHNO3,aHCl,aCH3SO3H,aHIO3,aDMA,aHIO2,CH2SO4old,CHNO3old,CHClold,CCH3SO3Hold,CHIO3old,CDMAold,CNH3old,CHIO2old,DHNO3,DHCl,DCH3SO3H,DHIO3,DHIO2,DDMA,Dair,DH2SO40,&
Keq1,Keq2,DH2SO4,mH2SO4,mHNO3,mHCl,mCH3SO3H,mHIO3,mDMA,mHIO2,m_air,speedH2SO4,speedHNO3,speedHCl,speedCH3SO3H,speedHIO3,speedHIO2,speedDMA,speedair,CH2SO4sat,&
CCH3SO3Hsat,CHIO3sat,CDMAsat,CH2SO4m,CHNO3m,CHClm,CCH3SO3Hm,CHIO3m,CHIO2m,CDMAm,CNH3m,CH2SO4tot,cNH3tot,CHNO3tot,CHCltot,CCH3SO3Htot,CHIO3tot,CDMAtot,&
Corgm,Corgtot, errorCNH3, fn,fprimn, dp_max, r1, r2, l_gas, dyn_visc, d_H2SO4, d_air, dens_air

REAL(dp), DIMENSION(NCOND,nr_bins) :: corgp,corgpold,korg,corgp_eq

REAL(dp), DIMENSION(NCOND), INTENT(in) :: Dorg
REAL(dp), DIMENSION(NCOND)  :: Corgold,speedorg,aorg,morg,m_org,dm_orgdt,m_orgold

REAL(dp), DIMENSION(nr_bins) :: cH2SO4pold,cHNO3pold,cHClpold,cCH3SO3Hpold,cHIO3pold,cHIO2pold,cDMApold,cNH3pold,KnH2SO4,&
KnHNO3,KnHCl,KnCH3SO3H,KnHIO3,KnHIO2,KnDMA,Knair,f_corH2SO4,f_corHNO3,f_corHCl,f_corCH3SO3H,f_corHIO3,f_corHIO2,f_corDMA,f_cororg,f_corair,DH2SO4eff,&
DHNO3eff,DHCleff,DCH3SO3Heff,DHIO3eff,DHIO2eff,DDMAeff,Daireff,kH2SO4,&
kHNO3,kHCl,kCH3SO3H,kHIO3,kHIO2,kDMA,kair,S_KelvinH2SO4,S_KelvinHNO3,S_KelvinHCl,S_KelvinCH3SO3H,S_KelvinHIO3,S_KelvinDMA,cH2SO4p,cHNO3p,cHClp,&
cCH3SO3Hp,cHIO3p,cDMAp,cHIO2p, term1, term2,cNap, cCO3p, cHCO3p, cOHp, cNH3p, charge, cNH4ion, cNH3aq,vp_dry,vp_wet, N_bins_fixed, &
Cunningh, Diff_p,m_p,speed_p, gasmeanfpH2SO4, gasmeanfpHNO3, gasmeanfpHCl, gasmeanfpCH3SO3H, gasmeanfpHIO3,  gasmeanfpHIO2, gasmeanfpDMA,gasmeanfporg, &
gasmeanfpair, Knorg, Dorgeff,Corgsat,xorg, N_bins_old,x_SVI

REAL(dp), DIMENSION(nr_bins+1) :: vp_fixed
INTEGER :: j,i,jj,a
REAL(dp), INTENT(in) :: dt

!!!! carlton addition for clusterin
logical, INTENT(IN):: clust_evap
REAL(DP), INTENT(OUT):: Nconc_evap1, comp_evap(:)

REAL(dp), DIMENSION(NSPEC_P,nr_bins) :: const1
logical, intent(in) :: flag_clusterin

!! initialize them
! type(clustering_mod), intent(out):: chem_1, chem_2

! chem_1%Nconc_evap=0D0; chem_2%Nconc_evap=0D0; 
comp_evap=0D0; Nconc_evap1=0D0

if (flag_clusterin) then

  N_bins_old=N_bins
  c_p_old=c_p
  dyn_visc=1.8D-5*(T/298.)**0.85                ! dynamic viscosity of air

  surf_tens(1:8)=(76.1-0.155*(T-273.15))*1D-3 ! (kg s^-2; % inorganics
  surf_tens(9:NSPEC_P)=0.05 ! organics
  surf_tens(9)=0.053 ! MSA
  surf_tens(10)=0.053 ! Assumed HIO3
  surf_tens(11)=0.053 ! Assumed DMA
  surf_tens(12)=0.053 ! Assumed HIO2

  DO j=1,NSPEC_P
  S_Kelvin(j,:)=1D0+2D0*surf_tens(j)*MX(j)/(d_p/2D0*dens*Rg*T) ! equilibrium saturation ratio of condensing gas (Kelvin effect)
  END DO

  l_gas=2D0*dyn_visc/(p*SQRT(8D0*Mair/(pi*Rg*T))) ! Gas mean free path in air (m)
  Cunningh=1D0+(2D0*l_gas/d_p)*(1.257+0.4*EXP(-1.1/(2D0*l_gas/d_p)))  ! Cunninghams correction factor
  Diff_p=Cunningh*kb*T/(3D0*pi*dyn_visc*d_p)                          ! Diffusivitys for the different particle sizes m^2/s
  m_p=(dens*pi*d_p**3D0)/6D0  ! mass of particles 
  speed_p=SQRT(8D0*kb*T/(pi*m_p)) ! speed of particles

  ! Condensation inorganics:
  cH2SO4pold=c_p(1,:)/Na*1D6 ! mol/m^3  
  cHNO3pold=c_p(2,:)/Na*1D6  ! mol/m^3  
  cHClpold=c_p(3,:)/Na*1D6   ! mol/m^3  
  cCH3SO3Hpold=c_p(9,:)/Na*1D6   ! mol/m^3  
  cHIO3pold=c_p(10,:)/Na*1D6   ! mol/m^3  
  cNH3pold=c_p(4,:)/Na*1D6   ! mol/m^3  
  cNap=c_p(5,:)/Na*1D6       ! mol/m^3 Na+
  cDMApold=c_p(11,:)/Na*1D6       ! mol/m^3 DMA
  cHIO2pold=c_p(12,:)/Na*1D6       ! mol/m^3 DMA

  corgpold(1:NCOND,:)=c_p(13:NSPEC_P,:)/Na*1D6 ! mol/m^3 cond org comp
  corgp=corgpold ! mol/m^3 cond org comp

  aH2SO4=aX(1);aHNO3=aX(2);aHCl=aX(3); aCH3SO3H=aX(9); aHIO3=aX(10); aDMA=aX(11);  aHIO2=aX(12); aorg=aX(13:NSPEC_P)

  CH2SO4old=cH2SO4/Na*1D6 ! vapor mole concentration (mol m^-3)
  CHNO3old=cHNO3/Na*1D6   ! vapor mole concentration (mol m^-3)
  CHClold=cHCl/Na*1D6     ! vapor mole concentration (mol m^-3)
  CNH3old=cNH3/Na*1D6     ! vapor mole concentration (mol m^-3)
  CCH3SO3Hold=cCH3SO3H/Na*1D6 ! vapor mole concentration (mol m^-3)     
  CHIO3old=cHIO3/Na*1D6 ! vapor mole concentration (mol m^-3)     
  CDMAold=cDMA/Na*1D6 ! vapor mole concentration (mol m^-3)
  CHIO2old=cHIO2/Na*1D6 ! vapor mole concentration (mol m^-3)

  Corgold=corg/Na*1D6      

  dens_air=Mair*p/(Rg*T) ! Air density
  DHNO3=(0.22/1.87)*1D-4 ! gas diffusivity HNO3 m^2/s
  DHCl=(0.22/1.42)*1D-4  ! gas diffusivity HCl m^2/s
  DCH3SO3H=1D-5       ! gas diffusivity CH3SO3H (MSA) m^2/s
  DHIO3=1D-5       ! gas diffusivity HIO3  m^2/s
  DHIO2=1D-5       ! gas diffusivity HIO3  m^2/s
  DDMA=1.5D-5       ! gas diffusivity DMA  m^2/s

  Dair=0.22D-4           ! Air molecule gas diffusivity m^2/s

  ! RH dependent Diffusion, diameter and mass of H2SO4 molecules:
  DH2SO40=0.094D-4       ! gas diffusivity H2SO4 m^2/s RH=0%
  Keq1=0.13D0
  Keq2=0.016D0
  DH2SO4 = (DH2SO40+0.85D0*DH2SO40*Keq1*RH*100.+0.76D0*DH2SO40*Keq1*Keq2*(RH*100.)**2D0)&
  /(1D0+Keq1*RH*100.+Keq1*Keq2*(RH*100.)**2D0) ! Diffusivity H2SO4 at ambient RH

  d_H2SO4=(((MX(1)/qX(1)/Na)*6D0/pi)**(1D0/3D0)+&
  Keq1*RH*100.*(((MX(1)+18D-3)/qX(1)/Na)*6D0/pi)**(1D0/3D0)+&
  Keq1*Keq2*(RH*100.)**2D0*(((MX(1)+36D-3)/qX(1)/Na)*6D0/pi)**(1D0/3D0))/&
  (1D0+Keq1*RH*100.+Keq1*Keq2*(RH*100.)**2D0) ! RH dependent H2SO4 diameter
  d_air=3D-10 ! Estimated molecular diameter of an average air molecule
  mH2SO4=(MX(1)/Na+Keq1*RH*100.*(MX(1)+18D-3)/Na+Keq1*Keq2*(RH*100.)**2D0*(MX(1)+36D-3)/Na)/&
  (1D0+Keq1*RH*100.+Keq1*Keq2*(RH*100.)**2D0) ! RH H2SO4 molecular mass

  mHNO3=MX(2)/Na; mHCl=MX(3)/Na ; mCH3SO3H=MX(9)/Na; mHIO3=MX(10)/Na; mDMA=MX(11)/Na; mHIO2=MX(12)/Na
  morg=MX(13:NSPEC_P)/Na
  m_air=Mair/Na

  speedH2SO4=SQRT(8D0*kb*T/(pi*mH2SO4)) ! speed of H2SO4 molecule
  speedHNO3=SQRT(8D0*kb*T/(pi*mHNO3))   ! speed of HNO3 molecules
  speedHCl=SQRT(8D0*kb*T/(pi*mHCl))     ! speed of HCl molecules
  speedCH3SO3H=SQRT(8D0*kb*T/(pi*mCH3SO3H)) ! speed of CH3SO3H molecules
  speedHIO3=SQRT(8D0*kb*T/(pi*mHIO3)) ! speed of HIO3 molecules
  speedHIO2=SQRT(8D0*kb*T/(pi*mHIO2)) ! speed of HIO3 molecules
  speedDMA=SQRT(8D0*kb*T/(pi*mDMA)) ! speed of DMA molecules

  speedorg=SQRT(8D0*kb*T/(pi*morg))     ! speed of organic molecules
  !speedair=SQRT(8D0*kb*T/(pi*m_air))     ! speed of air molecule
  speedair=4.625D2 ! Speed of air ions m/s at 294 K, Fuchs, 1963

  ! Gas mean free path, Knudsen number (Lehtinen and Kulmala, 2003):
  gasmeanfpH2SO4=3D0*(DH2SO4+Diff_p)/SQRT(speedH2SO4**2D0+speed_p**2D0) 
  gasmeanfpHNO3=3D0*(DHNO3+Diff_p)/SQRT(speedHNO3**2D0+speed_p**2D0)    
  gasmeanfpHCl=3D0*(DHCl+Diff_p)/SQRT(speedHCl**2D0+speed_p**2D0)       
  gasmeanfpCH3SO3H=3D0*(DCH3SO3H+Diff_p)/SQRT(speedCH3SO3H**2D0+speed_p**2D0)       
  gasmeanfpHIO3=3D0*(DHIO3+Diff_p)/SQRT(speedHIO3**2D0+speed_p**2D0)       
  gasmeanfpHIO2=3D0*(DHIO2+Diff_p)/SQRT(speedHIO2**2D0+speed_p**2D0)       
  gasmeanfpDMA=3D0*(DDMA+Diff_p)/SQRT(speedDMA**2D0+speed_p**2D0)   

  !gasmeanfpair=3D0*(Dair+Diff_p)/SQRT(speedair**2D0+speed_p**2D0)        
  gasmeanfpair=2.21D-8 !Mean free path of air ions (m), Hoppel and Frick, 1986

  KnH2SO4=2D0*gasmeanfpH2SO4/(d_p+d_H2SO4) ! Knudsen number H2SO4
  KnHNO3=2D0*gasmeanfpHNO3/(d_p+dX(2))     ! Knudsen number HNO3
  KnHCl=2D0*gasmeanfpHCl/(d_p+dX(3))       ! Knudsen number HCl
  KnCH3SO3H=2D0*gasmeanfpCH3SO3H/(d_p+dX(9)) ! Knudsen number MSA (CH3SO3H)
  KnHIO3=2D0*gasmeanfpHIO3/(d_p+dX(10))    ! Knudsen number HIO3
  KnHIO2=2D0*gasmeanfpHIO2/(d_p+dX(10))    ! Knudsen number HIO3
  KnDMA=2D0*gasmeanfpDMA/(d_p+dX(11)) ! Knudsen number DMA

  Knair=2D0*gasmeanfpair/(d_p+d_air)       ! Knudsen number H2SO4

  f_corH2SO4=(0.75*aH2SO4*(1D0+KnH2SO4))&
  /(KnH2SO4**2D0+KnH2SO4+0.283*KnH2SO4*aH2SO4+0.75*aH2SO4) ! Fuchs-Sutugin correction factor for transit and kinetic regime a=accomodation coefficient

  f_corHNO3=(0.75*aHNO3*(1D0+KnHNO3))&
  /(KnHNO3**2D0+KnHNO3+0.283*KnHNO3*aHNO3+0.75*aHNO3)  ! Fuchs-Sutugin correction factor for transit and kinetic regime

  f_corHCl=(0.75*aHCl*(1D0+KnHCl))/&
  (KnHCl**2D0+KnHCl+0.283*KnHCl*aHCl+0.75*aHCl) ! Fuchs-Sutugin correction factor for transit and kinetic regime

  f_corCH3SO3H=(0.75*aCH3SO3H*(1D0+KnCH3SO3H))/&
  (KnCH3SO3H**2D0+KnCH3SO3H+0.283*KnCH3SO3H*aCH3SO3H+0.75*aCH3SO3H) ! Fuchs-Sutugin correction factor for transit and kinetic regime

  f_corHIO3=(0.75*aHIO3*(1D0+KnHIO3))/&
  (KnHIO3**2D0+KnHIO3+0.283*KnHIO3*aHIO3+0.75*aHIO3) ! Fuchs-Sutugin correction factor for transit and kinetic regime
  
  f_corHIO2=(0.75*aHIO2*(1D0+KnHIO2))/&
  (KnHIO2**2D0+KnHIO2+0.283*KnHIO2*aHIO2+0.75*aHIO2) ! Fuchs-Sutugin correction factor for transit and kinetic regime

  f_corDMA=(0.75*aDMA*(1D0+KnDMA))/&
  (KnDMA**2D0+KnDMA+0.283*KnDMA*aDMA+0.75*aDMA) ! Fuchs-Sutugin correction factor for transit and kinetic regime

  f_corair=(0.75*1D0*(1D0+Knair))&
        /(Knair**2D0+Knair+0.283*Knair*1D0+0.75*1D0) ! Fuchs-Sutugin correction factor for transit and kinetic regime

  Daireff=(Dair+Diff_p)*f_corair                         ! m^2/s
  DH2SO4eff=(DH2SO4+Diff_p)*f_corH2SO4                   ! m^2/s
  DHNO3eff=(DHNO3+Diff_p)*f_corHNO3                      ! m^2/s
  DHCleff=(DHCl+Diff_p)*f_corHCl                         ! m^2/s
  DCH3SO3Heff=(DCH3SO3H+Diff_p)*f_corCH3SO3H             ! m^2/s
  DHIO3eff=(DHIO3+Diff_p)*f_corHIO3                      ! m^2/s
  DHIO2eff=(DHIO2+Diff_p)*f_corHIO2                      ! m^2/s
  DDMAeff=(DDMA+Diff_p)*f_corDMA                      ! m^2/s

  CH2SO4sat=0D0 ! H2SO4 saturation mole concentration (mol m^-3) p/(R*T)

  kH2SO4=N_bins*2D0*pi*(d_p+d_H2SO4)*DH2SO4eff    ! mass transfer coefficient s^-1
  kHNO3=N_bins*2D0*pi*(d_p+dX(2))*DHNO3eff        ! mass transfer coefficient s^-1
  kHCl=N_bins*2D0*pi*(d_p+dX(3))*DHCleff          ! mass transfer coefficient s^-1
  kCH3SO3H=N_bins*2D0*pi*(d_p+dX(9))*DCH3SO3Heff  ! mass transfer coefficient s^-1
  kHIO3=N_bins*2D0*pi*(d_p+dX(9))*DHIO3eff        ! mass transfer coefficient s^-1
  kHIO2=N_bins*2D0*pi*(d_p+dX(9))*DHIO2eff        ! mass transfer coefficient s^-1
  kDMA=N_bins*2D0*pi*(d_p+dX(11))*DDMAeff        ! mass transfer coefficient s^-1

  kair=N_bins*2D0*pi*(d_p+d_air)*Daireff          ! mass transfer coefficient s^-1


  DO j=1,NCOND
  gasmeanfporg=3D0*(Dorg(j)+Diff_p)/SQRT(speedorg(j)**2D0+speed_p**2D0) 
  Knorg=2D0*gasmeanfporg/(d_p+dX(8+j))       			 ! Knudsen number organic comp
  f_cororg=(0.75*aorg(j)*(1D0+Knorg))/&
  (Knorg**2D0+Knorg+0.283*Knorg*aorg(j)+0.75*aorg(j))  ! Fuchs-Sutugin correction factor for transit and kinetic regime
  Dorgeff=(Dorg(j)+Diff_p)*f_cororg                    ! m^2/s
  korg(j,:)=N_bins*2D0*pi*(d_p+dX(8+j))*Dorgeff        ! mass transfer coefficient s^-1
  END DO



  ! Condensation of organic compounds:
  corgpold=corgp;  

  DO j=1,nr_bins
  corgp_eq(:,j)=Corgold*(SUM(corgpold(:,j))+c_p(8,j)/Na*1D6)/(S_Kelvin(13:NSPEC_P,j)*psat_org/(Rg*T)*yorg(1:NCOND,j)) ! Approximate equilibrium concentration (mol/m^3) of each compound in each size bin   
  ! corgp_eq(:,j)=Corgold*(SUM(corgpold(:,j))+SUM(c_p(7:8,j)))/(S_Kelvin(12:NSPEC_P,j)*psat_org/(Rg*T)*yorg(1:NCOND,j)) ! Approximate equilibrium concentration (mol/m^3) of each compound in each size bin   
  END DO

  DO j=1,NCOND-1
  !IF (psat_org(j)<1D-2 .AND. corg(j)>1D4) THEN ! Only consider condensation of compounds with substantial concentration and relatively low vapour pressure
  IF ((corg(j)+SUM(c_p(12+j,:)))>1D6) THEN ! Only consider condensation/evaporation of compounds with substantial concentration gas+particles
  !xorg=corgpold(j,:)/(SUM(corgpold(1:NCOND,:),DIM=1)+SUM(c_p(7:8,:),DIM=1)/Na*1D6+1D-100) ! mole fraction of each organic compound in the organic + water particle phase
  xorg=corgpold(j,:)/(SUM(corgpold(1:NCOND,:),DIM=1)+c_p(13,:)/Na*1D6+1D-100) ! mole fraction of each organic compound in the organic + water particle phase
  Corgsat=yorg(j,:)*xorg*psat_org(j)/(Rg*T) 	! Saturation concentration (mol/m^3)  
  Corgm=(Corgold(j)+dt*SUM(korg(j,:)*S_Kelvin(12+j,:)*Corgsat))/(1D0+dt*SUM(korg(j,:)))
  Corgtot=Corgold(j)+SUM(corgpold(j,:)) ! total vapor + particle mole concentration (mol m^-3 air) 
  corgp(j,:)=corgpold(j,:)+dt*korg(j,:)*(MIN(Corgm,Corgtot)-S_Kelvin(12+j,:)*Corgsat)
  WHERE (corgp(j,:)<0D0) corgp(j,:)=0D0
  WHERE (corgp(j,:)>corgp_eq(j,:) .AND. Corgold(j)<psat_org(j)/(Rg*T)*S_Kelvin(12+j,:)*yorg(j,:)) corgp(j,:)=corgp_eq(j,:) ! Don't allow particles to grow more than to the saturation concentration limit   
  Corgm=Corgtot-SUM(corgp(j,:)) ! mol m^-3
  corg(j)=Corgm*1D-6*Na  ! Update gas-phase concentrations (molecules cm^-3)
  END IF
  END DO

  ! H2SO4:
  S_KelvinH2SO4=S_Kelvin(1,:)
  CH2SO4m=(CH2SO4old+dt*SUM(kH2SO4*S_KelvinH2SO4*CH2SO4sat))/(1D0+dt*SUM(kH2SO4)) ! condensation
  CH2SO4tot=CH2SO4old+SUM(cH2SO4pold) ! total vapor + particle mole concentration (mol m^-3 air) 
  cH2SO4p=cH2SO4pold+dt*kH2SO4*(MIN(CH2SO4m,CH2SO4tot)-S_KelvinH2SO4*CH2SO4sat)
  where (cH2SO4p<0D0) cH2SO4p=0D0
  CH2SO4m=CH2SO4tot-SUM(cH2SO4p) ! mol m^-3
  cH2SO4=CH2SO4m*1D-6*Na ! molecules cm^-3 (gas-phase H2SO4 concentration updated)
      
  S_KelvinHNO3=S_Kelvin(2,:);

  CS_H2SO4=SUM(kH2SO4);
  CS_air=SUM(kair)                              ! air ion condensation sink s^-1

  IF (thermodyn_index==1) THEN

  ! UPDATE HNO3(g):

  ! Solid salt NH4NO3 particles:
  !where (CHNO3s>0D0) term1=kHNO3*S_KelvinHNO3*CHNO3s*dt
  ! Liquid NH4NO3 particles:
  !where (CHNO3s<=0D0) term1=cHNO3pold*(1D0-EXP(-dt*S_KelvinHNO3*kHNO3/Kprim_HNO3))
  term1=cHNO3pold*(1D0-EXP(-dt*S_KelvinHNO3*kHNO3/Kprim_HNO3))

  ! Solid salt NH4NO3 particles:
  !where (CHNO3s>0D0) term2=kHNO3*dt
  ! Liquid NH4NO3 particles:
  !where (CHNO3s<=0D0) term2=Kprim_HNO3/S_KelvinHNO3*(1D0-EXP(-dt*S_KelvinHNO3*kHNO3/Kprim_HNO3))
  term2=Kprim_HNO3/S_KelvinHNO3*(1D0-EXP(-dt*S_KelvinHNO3*kHNO3/Kprim_HNO3))

  CHNO3m=(CHNO3old+SUM(term1))/(1D0+SUM(term2))
      
  ! UPDATE NO3(p):
  term1=0D0
  !where (CHNO3s>0D0) term1=cHNO3pold
  CHNO3tot=CHNO3old+SUM(term1) ! total vapor + particle mole concentration of solid NH4NO3 particles (mol m^-3 air) 

  !where (CHNO3s>0D0) cHNO3p=cHNO3pold+dt*kHNO3*(MIN(CHNO3m,CHNO3tot)-S_KelvinHNO3*CHNO3s)
  !where (CHNO3s<=0D0) cHNO3p=Kprim_HNO3*CHNO3m/S_KelvinHNO3+(cHNO3pold-Kprim_HNO3*CHNO3m/ &
  cHNO3p=Kprim_HNO3*CHNO3m/S_KelvinHNO3+(cHNO3pold-Kprim_HNO3*CHNO3m/ &
  S_KelvinHNO3)*EXP(-dt*S_KelvinHNO3*kHNO3/Kprim_HNO3)

  ! For safety, don't allow [NO3(p)]<0
  DO j=1,nr_bins
  IF (cHNO3p(j)<0D0) THEN
  cHNO3p(j)=0D0
  END IF
  END DO

  CHNO3tot=CHNO3old+SUM(cHNO3pold) ! total vapor + particle mole concentration (mol m^-3 air) 
  CHNO3m=CHNO3tot-SUM(cHNO3p) ! mol m^-3
  cHNO3=CHNO3m*1D-6*Na ! molecules cm^-3

        ! UPDATE HCl(g):
          term1=0D0; term2=0D0

          S_KelvinHCl=S_Kelvin(3,:);

          ! Solid salt NH4Cl particles:
          ! WHERE (CHCls>0D0) term1=kHCl*S_KelvinHCl*CHCls*dt
          ! ! Liquid NH4Cl particles:
          ! WHERE (CHCls<=0D0) term1=cHClpold*(1D0-EXP(-dt*S_KelvinHCl*kHCl/Kprim_HCl))
          term1=cHClpold*(1D0-EXP(-dt*S_KelvinHCl*kHCl/Kprim_HCl))
          ! ! Solid salt NH4Cl particles:
          ! WHERE (CHCls>0D0) term2=kHCl*dt
          ! ! Liquid NH4Cl particles:
          ! WHERE (CHCls<=0D0) term2=Kprim_HCl/S_KelvinHCl*(1D0-EXP(-dt*S_KelvinHCl*kHCl/Kprim_HCl))
          term2=Kprim_HCl/S_KelvinHCl*(1D0-EXP(-dt*S_KelvinHCl*kHCl/Kprim_HCl))
          CHClm=(CHClold+SUM(term1))/(1D0+SUM(term2))
      
          ! UPDATE Cl(p):
          !term1=0D0
          !WHERE (CHCls>0D0) term1=cHClpold
          !CHCltot=CHClold+SUM(term1) ! total vapor + particle mole concentration of solid NH4Cl particles (mol m^-3 air) 

          ! WHERE (CHCls>0D0) cHClp=cHClpold+dt*kHCl*(MIN(CHClm,CHCltot)-S_KelvinHCl*CHCls)
          !  WHERE (CHCls<=0D0) cHClp=Kprim_HCl*CHClm/S_KelvinHCl+(cHClpold-Kprim_HCl*CHClm/ &
          ! S_KelvinHCl)*EXP(-dt*S_KelvinHCl*kHCl/Kprim_HCl)
          cHClp=Kprim_HCl*CHClm/S_KelvinHCl+(cHClpold-Kprim_HCl*CHClm/ &
          S_KelvinHCl)*EXP(-dt*S_KelvinHCl*kHCl/Kprim_HCl)

          ! For safety, don't allow [Cl(p)]<0
          DO j=1,nr_bins
          IF (cHClp(j)<0D0) THEN
              cHClp(j)=0D0
          END IF
          END DO

          CHCltot=CHClold+SUM(cHClpold) ! total vapor + particle mole concentration (mol m^-3 air) 
          CHClm=CHCltot-SUM(cHClp) ! mol m^-3
          cHCl=CHClm*1D-6*Na ! molecules cm^-3
          ! For safety, if negative gas-phase concentrations use use previous time step conc:
          IF (cHCl<0D0) THEN
            cHCl=CHClold*1D-6*Na ! molecules cm^-3
            cHClp=cHClpold
          END IF
      

  ! UPDATE MSA(g) (CH3SO3H(g):
  term1=0D0; term2=0D0
  S_KelvinCH3SO3H=S_Kelvin(9,:);

  ! Liquid CH3SO3- particles:
  term1=cCH3SO3Hpold*(1D0-EXP(-dt*S_KelvinCH3SO3H*kCH3SO3H/Kprim_CH3SO3H))

  ! Liquid CH3SO3- particles:
  term2=Kprim_CH3SO3H/S_KelvinCH3SO3H*(1D0-EXP(-dt*S_KelvinCH3SO3H*kCH3SO3H/Kprim_CH3SO3H))

  CCH3SO3Hm=(CCH3SO3Hold+SUM(term1))/(1D0+SUM(term2))
      
  ! UPDATE CH3SO3(p):
  term1=0D0
  cCH3SO3Hp=Kprim_CH3SO3H*CCH3SO3Hm/S_KelvinCH3SO3H+(cCH3SO3Hpold-Kprim_CH3SO3H*CCH3SO3Hm/ &
  S_KelvinCH3SO3H)*EXP(-dt*S_KelvinCH3SO3H*kCH3SO3H/Kprim_CH3SO3H)

  ! For safety, don't allow [CH3SO3(p)]<0
  DO j=1,nr_bins
  IF (cCH3SO3Hp(j)<0D0) THEN
  cCH3SO3Hp(j)=0D0
  END IF
  END DO

  CCH3SO3Htot=CCH3SO3Hold+SUM(cCH3SO3Hpold) ! total vapor + particle mole concentration (mol m^-3 air) 
  CCH3SO3Hm=CCH3SO3Htot-SUM(cCH3SO3Hp) ! mol m^-3
  cCH3SO3H=CCH3SO3Hm*1D-6*Na ! molecules cm^-3

  !! Alternative solution assuming effectively non-volatile MSA:
  !CCH3SO3Hsat=0D0
  !S_KelvinCH3SO3H=S_Kelvin(9,:);
  !CCH3SO3Hm=(CCH3SO3Hold+dt*SUM(kCH3SO3H*S_KelvinCH3SO3H*CCH3SO3Hsat))/(1D0+dt*SUM(kCH3SO3H)) ! condensation
  !CCH3SO3Htot=CCH3SO3Hold+SUM(cCH3SO3Hpold) ! total vapor + particle mole concentration (mol m^-3 air) 
  !cCH3SO3Hp=cCH3SO3Hpold+dt*kCH3SO3H*(MIN(CCH3SO3Hm,CCH3SO3Htot)-S_KelvinCH3SO3H*CCH3SO3Hsat)
  !where (cCH3SO3Hp<0D0) cCH3SO3Hp=0D0
  !CCH3SO3Hm=CCH3SO3Htot-SUM(cCH3SO3Hp) ! mol m^-3
  !cCH3SO3H=CCH3SO3Hm*1D-6*Na ! molecules cm^-3 (gas-phase H2SO4 concentration updated)


  ! UPDATE HIO3(g):
  term1=0D0; term2=0D0
  S_KelvinHIO3=S_Kelvin(10,:);

  ! Liquid IO3- particles:
  term1=cHIO3pold*(1D0-EXP(-dt*S_KelvinHIO3*kHIO3/Kprim_HIO3))

  ! Liquid IO3- particles:
  term2=Kprim_HIO3/S_KelvinHIO3*(1D0-EXP(-dt*S_KelvinHIO3*kHIO3/Kprim_HIO3))

  CHIO3m=(CHIO3old+SUM(term1))/(1D0+SUM(term2))
      
  ! UPDATE HIO3(p):
  term1=0D0
  cHIO3p=Kprim_HIO3*CHIO3m/S_KelvinHIO3+(cHIO3pold-Kprim_HIO3*CHIO3m/ &
  S_KelvinHIO3)*EXP(-dt*S_KelvinHIO3*kHIO3/Kprim_HIO3)

  ! For safety, don't allow [IO3(p)]<0
  DO j=1,nr_bins
  IF (cHIO3p(j)<0D0) THEN
  cHIO3p(j)=0D0
  END IF
  END DO

  CHIO3tot=CHIO3old+SUM(cHIO3pold) ! total vapor + particle mole concentration (mol m^-3 air) 
  CHIO3m=CHIO3tot-SUM(cHIO3p) ! mol m^-3
  cHIO3=CHIO3m*1D-6*Na ! molecules cm^-3

  ! DMA:
  S_KelvinDMA=S_Kelvin(11,:)
  cDMAm=(CDMAold+dt*SUM(kDMA*S_KelvinDMA*CDMAsat))/(1D0+dt*SUM(kDMA)) ! condensation
  CDMAtot=CDMAold+SUM(cDMApold) ! total vapor + particle mole concentration (mol m^-3 air) 
  cDMAp=cDMApold+dt*kDMA*(MIN(CDMAm,CDMAtot)-S_KelvinDMA*CDMAsat)
  where (cDMAp>cH2SO4p+cHNO3p+cHClp+cCH3SO3Hp-cNap) cDMAp=cH2SO4p+cHNO3p+cHClp+cCH3SO3Hp-cNap
  where (cDMAp<0D0) cDMAp=0D0
  CDMAm=CDMAtot-SUM(cDMAp) ! mol m^-3
  cDMA=CDMAm*1D-6*Na ! molecules cm^-3 (gas-phase DMA concentration updated)
      
      
  ! Calculate the new gas phase NH3 concentration after equilibration with
  ! the particle phase according to the Newton-Raphson iteration in Jacobson, 2005 (eq. 17.126)
  cCO3p=mCO3*W*N_bins ! mol/m^3 CO32-
  cHCO3p=mHCO3*W*N_bins ! mol/m^3 HCO3-
  cOHp=mOH*W*N_bins ! mol/m^3 OH-

  ! charge=-fHSO4*cH2SO4p-2D0*fSO4*cH2SO4p-fCl*cHClp-fNO3*cHNO3p-fCH3SO3*cCH3SO3Hp-fHIO3*cHIO3p+cNap+cDMAp ! charge imbalance after condensation of HNO3, HCl, H2SO4, MSA (CH3SO3H), HIO3 and DMA
  charge=-fHSO4*cH2SO4p-2D0*fSO4*cH2SO4p-fCl*cHClp-fNO3*cHNO3p-fCH3SO3*cCH3SO3Hp+cNap+cDMAp ! charge imbalance after condensation of HNO3, HCl, H2SO4, MSA (CH3SO3H), HIO3 and DMA

  where (charge>0D0) charge=0D0 ! correct for positive charge imbalance which is only possible if the acids evaporates

  CNH3tot=CNH3old+sum(cNH3pold) ! total vapor + particle mole concentration (mol m^-3 air) 
  CNH3m=0D0 ! initial guess of gas phase NH3
  errorCNH3=1D0

  DO WHILE (errorCNH3>1D-6)
  fn=CNH3m+SUM(CNH3m*Hprim_NH3-charge*CNH3m*Hprim_NH3*Kprim_NH3/ &
  (CNH3m*Hprim_NH3*Kprim_NH3+1D0))-CNH3tot
  fprimn=1D0+SUM(Hprim_NH3-charge*Hprim_NH3*Kprim_NH3/(CNH3m*Hprim_NH3*Kprim_NH3+1D0)+ &
  charge*CNH3m*(Hprim_NH3*Kprim_NH3)**2D0/((CNH3m*Hprim_NH3*Kprim_NH3+1D0)**2D0))
  CNH3m=CNH3m-fn/fprimn
  errorCNH3=ABS(fn/fprimn/CNH3m) ! error mol/m^3
  END DO
  cNH3=CNH3m*1D-6*Na ! molecules cm^-3

  cNH4ion=-charge*CNH3m*Hprim_NH3*Kprim_NH3/(CNH3m*Hprim_NH3*Kprim_NH3+1D0) ! ammount of dissolved NH4+ mol/m^3 air in each size bin
  cNH3aq=CNH3m*Hprim_NH3 ! amount of dissolved NH3(aq) mol/m^3 air in each size bin
  cNH3p=cNH4ion+cNH3aq ! amount of NH3(aq) and NH4+ mol/m^3 air in each size bin

  c_p(2,:)=cHNO3p*1D-6*Na
  c_p(3,:)=cHClp*1D-6*Na
  c_p(4,:)=cNH3p*1D-6*Na
  c_p(9,:)=cCH3SO3Hp*1D-6*Na
  c_p(10,:)=cHIO3p*1D-6*Na
  c_p(11,:)=cDMAp*1D-6*Na
  END IF


  c_p(1,:)=cH2SO4p*1D-6*Na
  c_p(13:NSPEC_P,:)=corgp(1:NCOND,:)*1D-6*Na
end if !! if flag_clusterin


   ! Output particle properties:
      DO j = 1,nr_bins
        vp_dry(j)=SUM(c_p(index_dry,j)/Na*MX(index_dry)/qX(index_dry)/(N_bins(j)*1D-6)) ! m^3     
      END DO
      vp_fixed(1:nr_bins)=dp_dry**3D0*pi/6D0 ! dry single particle particle volume
      N_bins_fixed=0D0
      c_p_fixed=0D0
      dp_max = dp_dry(nr_bins)*dp_dry(nr_bins)/dp_dry(nr_bins-1)
      vp_fixed(nr_bins+1)=(pi*dp_max**3D0)/6D0

      IF (fullstat_index==1) THEN ! Full-stationary 
        DO j = 1,nr_bins
          a = MINLOC(vp_fixed-vp_dry(j),1,mask = (vp_fixed-vp_dry(j)) >= 0D0)
          IF (a > j) THEN ! growth
            r1 = (vp_fixed(a)-vp_dry(j))/(vp_fixed(a)-vp_fixed(a-1)) ! Fraction of particles in size bin a-1
            r2 = 1D0-r1 ! Fraction of particles in next size bin (a)
            IF (a > nr_bins+1) THEN
            ELSE IF (a == nr_bins+1) THEN
              N_bins_fixed(a-1) = N_bins_fixed(a-1)+r1*N_bins(j)
              c_p_fixed(:,a-1) = c_p_fixed(:,a-1)+vp_fixed(a-1)/vp_dry(j)*r1*c_p(:,j)
            ELSE
              N_bins_fixed(a-1) = N_bins_fixed(a-1)+r1*N_bins(j)
              N_bins_fixed(a) = N_bins_fixed(a)+r2*N_bins(j)
              c_p_fixed(:,a-1) = c_p_fixed(:,a-1)+vp_fixed(a-1)/vp_dry(j)*r1*c_p(:,j)
              c_p_fixed(:,a) = c_p_fixed(:,a)+vp_fixed(a)/vp_dry(j)*r2*c_p(:,j)
            END IF
          ELSE ! evaporation
            IF (a == 0) THEN           
            ELSEIF (a == 1) THEN
              r1 = (vp_dry(j))/(vp_fixed(1)) ! Fraction of particles in size bin a
              N_bins_fixed(1) = N_bins_fixed(1)+r1*N_bins(j)
              c_p_fixed(:,1) = c_p_fixed(:,1)+r1*c_p(:,j)
              ! c_p_fixed(:,1) = c_p_fixed(:,1)+vp_fixed(1)/vp_dry(j)*r1*c_p(:,j)
              ! write(*,*) 'Carlton debug L1460, a, j and r1', '', a, '', j, '', r1,'',vp_fixed(1)/vp_dry(j)*r1, flag_clusterin
                if(clust_evap) then
                  r2 = 1D0-r1
                  Nconc_evap1 = Nconc_evap1 + r2*N_bins(j)
                  comp_evap = comp_evap+ r2*c_p(:,j) !molec/cm3
                end if
           
            ELSE
              r1 = (vp_dry(j)-vp_fixed(a-1))/(vp_fixed(a)-vp_fixed(a-1)) ! Fraction of particles in size bin i
              r2 = 1D0-r1 ! Fraction of particles in previous size bin (i-1)
              N_bins_fixed(a) = N_bins_fixed(a)+r1*N_bins(j)
              N_bins_fixed(a-1) = N_bins_fixed(a-1)+r2*N_bins(j)
              c_p_fixed(:,a) = c_p_fixed(:,a)+vp_fixed(a)/vp_dry(j)*r1*c_p(:,j)
              c_p_fixed(:,a-1) = c_p_fixed(:,a-1)+vp_fixed(a-1)/vp_dry(j)*r2*c_p(:,j)
              ! if (a==2 .or. a==3) then
              ! !   if (clust_evap) then
              ! ! !     ! Nconc_evap2 = Nconc_evap2 + r2*N_bins(j)
              !       comp_evap(11) = c_p_fixed(11,a) !comp_evap+ r2*c_p(:,j) !molec/cm3
              ! !   end if     
              ! end if
              ! end if
            END IF
          END IF
        END DO
        c_p=c_p_fixed
        N_bins=N_bins_fixed
      END IF
       DO j = 1,nr_bins
        IF (N_bins(j)<1D-6) THEN ! This can occur with the full-stationary method if some particles grow to more than the next particle size bin within one condensation time step
          IF (N_bins(j)<0D0) THEN
            WRITE(*,*) 'Neg. conc. in bin ',j,' after cond.:',N_bins(j),' m^-3'
          END IF
          
          N_bins(j)=1D-3		
          c_p(:,j)=c_p_backg(:,j)*1D-3;       
          c_p(:,j)=c_p(:,j)*vp_fixed(j)/SUM(c_p(index_dry,j)/Na*MX(index_dry)/qX(index_dry))*N_bins(j)*1D-6   
        END IF
        DO i = 1,NSPEC_P
          IF (c_p(i,j)<0D0) THEN
              WRITE(*,*) 'Neg. conc. for species ',i,' in bin ',j,' after cond.:',c_p(i,j),' cm^-3'
            c_p(i,j) = c_p_backg(i,j)*1D-3
          END IF
        END DO
        
        vp_dry(j)=SUM(c_p(index_dry,j)/Na*MX(index_dry)/qX(index_dry)/(N_bins(j)*1D-6)) ! m^3     
        vp_wet(j)=SUM(c_p(:,j)/Na*MX/qX/(N_bins(j)*1D-6)) ! m^3
        V_bins(j)=N_bins(j)*vp_wet(j)
        dens(j)=SUM(c_p(:,j)*MX*1D6)/Na/V_bins(j) ! Total particle density
      END DO
      dp_dry=(vp_dry*6D0/pi)**(1D0/3D0) ! Dry particle diameters
      d_p=(vp_wet*6D0/pi)**(1D0/3D0) ! Dry particle diameters 

END SUBROUTINE condensation

    !-----------------------------------------------------------------------!
    ! Dry deposition gases (SO2,O3,NO2,NO,HNO3,H2O2,HCHO,CH3OOH,HONO,H2SO4) !
    !-----------------------------------------------------------------------!

   SUBROUTINE dry_dep_gases(temp,press,RH,u,PBLH,SHTF,DSWF,landuse_index,&
    vd_gas,rain, windspeed,season_index,Henrys_coeff,Dorg)
    
      REAL(dp), INTENT(in) :: temp,press,RH,u,PBLH,SHTF,DSWF,rain,windspeed ! meteorological parameters:
      ! temp [K],pressure [Pa],RH,u-resp v-comp of mom flux [N/m^2],boundary layer height [m], sensible heat net flux at surface [W/m^2] (GDAS) and downward short wave rad fux (W/m²),rain(mm/h)
      ! and time along trajectory in hour and trajectory time step
      INTEGER, INTENT(in)  :: landuse_index,season_index
      REAL(dp), DIMENSION(NCOND), INTENT(in) :: Henrys_coeff, Dorg
      REAL(dp), DIMENSION(NCOND+12), INTENT(out) :: vd_gas ! deposition loss-rate coefficient for each particels size bin [m/s]
      REAL(dp), DIMENSION(NCOND+12)              :: rb,rc
      REAL(dp) :: ra
      ! empirical coefficients
      REAL(dp), DIMENSION(8,5) :: z0,rj,rlu,rac,rgsS,rgsO,rclS,rclO
      ! Variables used to calculate the air density
      REAL(dp) :: R_d,R_v,L,ratio,p_0,T_0,p_sat,p_v_H2O,ratio_gasconst,mix_ratio,T_v,dens_air
      ! Variables needed to calculate the Richardson nr
      REAL(dp) :: zr,t_s,t_s_pot,ka,g,C_p_dry,C_p,T_star,L_Ob,Rfr,Rf0
      ! Variables used to calculate aerodynamic resistance
      REAL(dp) :: y0,yr,Pr,beta,gam,denominator
      ! Variables used to calculate the qausi-laminar resistance
      REAL(dp) :: v
      REAL(dp), DIMENSION(NCOND+12) :: Diff,Sc,D_ratio
      REAL(dp) :: z_rough
      ! VAriables needed to calculate the surface resistance
      REAL(dp) :: rp,rdc,kg,Sc_CO2
      REAL(dp), DIMENSION(NCOND+12) :: rst,rlui,rcli,rgsi,H_eff,f0,H_T0,H_exp,H,H_dimless,Sc_water,Diff_water,Sc_ratio,kl,alpha,bet
      REAL(dp) :: Dg0,Keq1,Keq2,Dg1,Dg2
      !REAL(4) :: DH2SO4, DH2O
      REAL(dp) :: DH2SO4, DH2O
      INTEGER  :: i,j

      ! Rougness length for diffreent land-use catagories and seasons, also table 19.2 [m]
      z0(1,:) = (/ 0.8,0.9,0.9,0.9,0.8 /) ! evergreen, needleleaf trees
      z0(2,:) = (/ 1.05,1.05,0.95,0.55,0.75 /) ! deciduous broadleaf trees
      z0(3,:) = (/ 0.1,0.1,0.05,0.02,0.05 /) ! grass
      z0(4,:) = (/ 0.04,0.04,0.04,0.04,0.04 /) ! desert
      z0(5,:) = (/ 0.1,0.1,0.1,0.1,0.1 /) ! shrubs and interrupted woodlands
      z0(6,:) = (/ 1D-3,1D-3,1D-3,1D-3,1D-3 /) ! sea, from table 8.1 in Jacobson
      z0(7,:) = (/ 0.26,0.26,0.26,0.26,0.26 /) ! small urban area, from table 8.1 in Jacobson
      z0(8,:) = (/ 1D-3,1D-3,1D-3,1D-3,1D-3 /) ! Snow and ice
      
      ! Resistance components used when calculating the surface resistance for gases, table 19.3 Seinfeld and Pandis (9999 indicate that there is no air-surface exchange via that resistance pathway
      ! The minimum, bulk canopy stomatal resistance for water vapor:
      rj(1,:) = (/ 130.,250.,250.,400.,250. /) ! Evergreen, needleleaf
      rj(2,:) = (/ 70D0,1D10,1D10,1D10,140D0 /) ! Deciduous broadleaf
      rj(3,:) = (/ 80D0,1D10,1D10,1D10,160D0 /) ! grass (nonforested wetland in table 19.3)
      rj(4,:) = (/ 1D10,1D10,1D10,1D10,1D10 /) ! desert (barren land, mostly desert)
      rj(5,:) = (/ 150D0,1D10,1D10,1D10,300D0 /) ! shrubs (rocky open aeras with low-growing shrubs)
      rj(6,:) = (/ 1D10,1D10,1D10,1D10,1D10 /) ! sea (water, both salt and fresh)
      rj(7,:) = (/ 1D10,1D10,1D10,1D10,1D10 /) ! urban
      rj(8,:) = (/ 1D10,1D10,1D10,1D10,1D10 /) ! snow and ice
      ! The resistance of the outer surfaces in the upper canopy
      rlu(1,:) = (/ 2000.,4000.,4000.,6000.,2000. /) ! Evergreen
      rlu(2,:) = (/ 2000D0,9000D0,9000D0,1D10,4000D0 /) ! Deciduous
      rlu(3,:) = (/ 2500D0,9000D0,9000D0,1D10,4000D0 /) ! Wetland
      rlu(4,:) = (/ 1D10,1D10,1D10,1D10,1D10 /) ! barren
      rlu(5,:) = (/ 4000.,9000.,9000.,9000.,800. /) ! shrubs
      rlu(6,:) = (/ 1D10,1D10,1D10,1D10,1D10 /) ! water
      rlu(7,:) = (/ 1D10,1D10,1D10,1D10,1D10 /) ! urban
      rlu(8,:) = (/ 1D10,1D10,1D10,1D10,1D10 /) ! snow and ice
      ! transfer resistance on the ground (that depends only on canopy height)
      rac(1,:) = (/ 2000.,2000.,2000.,2000.,2000. /) ! evergreen
      rac(2,:) = (/ 2000.,1500.,1000.,1000.,1200. /) ! deciduous
      rac(3,:) = (/ 300.,200.,100.,50.,200. /) ! wetland
      rac(4,:) = (/ 1D-3,1D-3,1D-3,1D-3,1D-3 /) ! barren
      rac(5,:) = (/ 200.,140.,120.,50.,120. /) ! shrubs
      rac(6,:) = (/ 1D-3,1D-3,1D-3,1D-3,1D-3 /) ! water
      rac(7,:) = (/ 100.,100.,100.,100.,100. /) ! urban
      rac(8,:) = (/ 100.,100.,100.,100.,100. /) ! snow and ice
      ! resistance for uptake by soil, leaf litter, and so on at the ground SO2
      rgsS(1,:) = (/ 500.,500.,500.,100.,500. /) ! evergreen
      rgsS(2,:) = (/ 500.,500.,500.,100.,500. /) ! deciduous
      rgsS(3,:) = (/ 1D-3,1D-3,1D-3,1D2,1D-3 /) ! wetland
      rgsS(4,:) = (/ 1000.,1000.,1000.,1000.,1000. /) ! barren
      rgsS(5,:) = (/ 400.,400.,400.,50.,400. /) ! shrubs
      rgsS(6,:) = (/ 1D-3,1D-3,1D-3,1D-3,1D-3 /) ! water
      rgsS(7,:) = (/ 400.,400.,400.,100.,500. /) ! urban
      rgsS(8,:) = (/ 100.,100.,100.,100.,100. /) ! ! snow and ice
      ! restistance for uptake by soil, leaf litter, and so on at the ground, O3
      rgsO(1,:) = (/ 200.,200.,200.,3500.,200. /) ! evergreen
      rgsO(2,:) = (/ 200.,200.,200.,3500.,200. /) ! deciduous
      rgsO(3,:) = (/ 1000.,800.,1000.,3500.,1000. /) ! wetland
      rgsO(4,:) = (/ 400.,400.,400.,400.,400. /) ! barren
      rgsO(5,:) = (/ 200.,200.,200.,3500.,200. /) ! shrubs
      !rgsO(6,:) = (/ 2000.,2000.,2000.,2000.,2000. /) ! water
	  rgsO(6,:) = (/ 1D4,1D4,1D4,1D4,1D4 /) ! water (If model is run in polar regions, see e.g. Pound et al. https://doi.org/10.5194/acp-2019-1043)
      rgsO(7,:) = (/ 300.,300.,300.,600.,300. /) ! urban
      rgsO(8,:) = (/ 1D6,1D6,1D6,1D6,1D6 /) ! snow and ice (assumed)
      ! resistance for uptake by leaves,twigs, and other exposed surfaces, SO2
      rclS(1,:) = (/ 2000.,2000.,3000.,200.,2000. /) ! evergreen
      rclS(2,:) = (/ 2000.,9000.,9000.,9000.,4000. /) ! deciduous
      rclS(3,:) = (/ 2500.,9000.,9000.,9000.,4000. /) ! wetland
      rclS(4,:) = (/ 1D10,1D10,1D10,1D10,1D10 /) ! barren
      rclS(5,:) = (/ 4000.,9000.,9000.,9000.,8000. /) ! shrubs
      rclS(6,:) = (/ 1D10,1D10,1D10,1D10,1D10 /) ! water
      rclS(7,:) = (/ 1D10,1D10,1D10,1D10,1D10 /) ! urban
      rclS(8,:) = (/ 1D10,1D10,1D10,1D10,1D10 /) ! snow and ice
      ! resistance for uptake by leaves, twigs, and other exposed surfaces, O3
      rclO(1,:) = (/ 1000.,1000.,1000.,1500.,1500. /) ! evergreen    
      rclO(2,:) = (/ 1000.,400.,400.,400.,500. /) ! deciduous
      rclO(3,:) = (/ 1000.,400.,800.,800.,600./) ! wetland
      rclO(4,:) = (/ 1D10,1D10,1D10,1D10,1D10 /) ! barren
      rclO(5,:) = (/ 1000.,400.,600.,800.,800. /) ! shrubs
      rclO(6,:) = (/ 1D10,1D10,1D10,1D10,1D10 /) ! water
      rclO(7,:) = (/ 1D10,1D10,1D10,1D10,1D10 /) ! urban
      rclO(8,:) = (/ 1D10,1D10,1D10,1D10,1D10 /) ! snow and ice
   
      ! Calculation of aerodynamic resistance      
   
      ! Calculation of air density using the Clasius-Clapeyron eq to calculate the saturation vapor pressure 
      R_d = 287.053 ! specific gas constant for dry air (J/(kgK))
      R_v = 461.495 ! specific gas constant for water vapor (J/(kgK))
      L = 2.5D6 ! latent heat of vaporization (J/kg)
      ratio = L/R_v ! (K)
      p_0 = 0.611 ! (kPa)
      T_0 = 273.15 ! (K)
      p_sat = p_0*exp(ratio*(1D0/T_0-1D0/temp)) ! saturation vapor pressure of water (kPa)
      p_v_H2O = ((RH/1D2)*p_sat)*1D3 ! vapor pressure of water (Pa)
      ratio_gasconst = 0.622 ! ratio of gas constant for dry air to that of water vapor (g_v/g_dry)
      mix_ratio = ratio_gasconst*p_v_H2O/(press-p_v_H2O) ! mixing ratio
      T_v = temp*(1D0+0.61*mix_ratio) ! Virtual temperature (K)
      dens_air = Mair*press/(Rg*temp)  ! Air density
      ka = 0.4 ! Von Karman constant
  
      ! Calculation of the friction velocity using vertical turbulent flux of
      ! horizontal momentum (as computed in HYSPLIT)
  !    flux_horiz = dens_air*((-UMOF/dens_air)**2.+(-VMOF/dens_air)**2.)**0.5 ! Vertical turbulent flux of horizontal momentum
  !    u = (flux_horiz/dens_air)**0.5;            ! Friction velocity (m/s)
        
      ! Calculation of Rfr = z/L (L is the Obukhov length) computed from the friction values to be consistent with model derived flux fields
      ! Diffusion coefficients of selected gases
      DH2O = 0.234D-4 ! Diffusion coefficient of water vapor in air (m²/s), table 16.2 Seinfeld and Pandis
      Dg0=0.094D-4 ! gas diffusivity H2SO4 m^2/s RH=0%
      Keq1=0.13
      Keq2=0.16
      Dg1=0.85*Dg0
      Dg2=0.76*Dg0
      DH2SO4=(Dg0+Dg1*Keq1*RH+Dg2*Keq1*Keq2*RH**2.)/(1+Keq1*RH+Keq1*Keq2*RH**2.) ! Diffusivity H2SO4


      D_ratio(1:12) = (/ 1.89D0,1.63D0,1.6D0,1.29D0,1.87D0,1.37D0,1.29D0,1.6D0,1.62D0,DH2O/DH2SO4,1.42D0,0.97D0/) ! ratio between diff of water vapor and specific gas: so2,o3,no2,no,hno3,h2o2,hcho,ch3ooh,hono,h2so4,hcl,nh3, HOM monomer, HOM dimer from table 19.4 Seinfeld and Pandis, h2so4 from parameterized value of h2so4 diff in condensation subroutine, hcl assumed same as HNO3
   
      D_ratio(13:NCOND+12)=DH2O/Dorg

      Diff = DH2O/D_ratio ! Diffusivity of the specific gases (see line above)

      z_rough = z0(landuse_index,season_index)
      !z_rough = Diff/(ka*u)
      zr = 10. ! reference height where dry dep occur. ?? 0.1.*PBLH; % estimated height of the surface layer (bottom 10 % of PBL) (m)
!     IF (zr>1D2) THEN ! if using zr = 0.1*PBLH
!       zr=1D2
!     END IF
  
      t_s = T_v ! temperature 2 m agl taken as temperature at z_s (K)
      t_s_pot = T_v*(1D5/press)**0.286
      g = 9.81 ! gravitational acceleration (m/s^2)
      C_p_dry = 1004.67 ! specific heat (J/kgK)
      C_p = C_p_dry*(1D0+0.84*mix_ratio) ! specific heat, moist air (J/kgK)
      T_star = -SHTF/(C_p*dens_air*u) ! friction temperature (K)
      L_Ob = u**2.*t_s_pot/(ka*g*T_star)

      Rfr = zr/L_Ob  ! z/L; if > 0 stable, if < 0 unstable, if = 0 neutral
      !DO j = 1,NCOND+12
       Rf0 = z_rough/L_Ob  ! Richards number ground

      ! Calculation of aerodynamic resistance (Jacobson  8.4.2.3 and p.667 eq 20.12) 
       Pr = 0.95 ! turbulent Prandtl number (when ka = 0.4 (Hogstrom, 1988))
       beta = 7.8 ! when ka = 0.4 (Hogstrom, 1988)
       gam = 11.6 ! when ka = 0.4 (Hogstrom, 1988)
       denominator = ka*u
       IF (Rfr > 1D-3) THEN
        ra = (Pr*log(zr/z_rough)+beta/L_Ob*(zr-z_rough))/denominator ! stable
       ELSE IF (Rfr<1D-3 .AND. Rfr>-1D-3) THEN
        ra = Pr*log(zr/z_rough)/denominator ! neutral
       ELSE
        ra = (Pr*(log(((1.-gam*zr/L_Ob)**0.5-1.)/((1.-gam*zr/L_Ob)**0.5+1.))-log(((1.-gam*z_rough/L_Ob)**0.5-1.)/((1.-gam*z_rough/L_Ob)**0.5+1.))))/(ka*u) ! unstable
       END IF
      !END DO 

      ! Calculation of Quasi-laminar resistance (rb)                 

!     Vorg = Morg/qHC/Na ! molecule volume (m^3 molec)
!     d_org = Vorg**(1./3.) ! Molecular diameter (m) Bird, R. B.; Stewart, W. E.; Lightfoot, E. N. Transport Phenomena; John Wiley & Sons: New York, 1960 p 22
!     Dorg=5./(16.*Na*d_org**2.*dens_air)*((Rg*temp*Mair/(2.*pi)*((Morg+Mair)/Morg)))**0.5 ! Diffusivity  m^2 s^-1
      
      v = 1.34D-5 ! (m²/s) kinematic viscosity for air at 273.15 K
      Sc = v/Diff ! dimensionless Schmidt number for each gas considered
!      IF (landuse_index == 6) THEN ! quasi-laminar resistance over ocean (Hicks and Liss, 1976)
!        rb = 1/(ka*u)*log(z_rough*ka*u/Diff)
!      ELSE   
!        rb = 5.*Sc**(2./3.)/u ! quasi-laminar resistance for each gas considered
!      END IF   
      rb = 5.*Sc**(2./3.)/u ! quasi-laminar resistance for each gas considered
      ! Calculation of surface resistance (rc)                         

    ! IF (landuse_index == 6 .OR. landuse_index == 8) THEN ! surface resistance for water (19.5.1 Seinfeld and Pandis)
      ! kg = 0.0013*windspeed ! Hicks and Liss 1976, gas phase mass transfer coeff
      ! H_T0(1:12) = (/ 1.23D0,1.1D-2,1D-2,1.9D-3,2.1D5,1D5,2.5D0,310D0,49D0,2.1D5,1.1D0,62D0 /) ! Henry's law coeff. at T=298K from table 7.2 Seinfeld and Pandis (mol L⁻¹ atm⁻¹)
      ! H_exp(1:12) = (/ -3020.,2300.,-2500.,-1480.,-8650.,-6800.,-6500.,-5600.,-4800.,-8650.,-2020.,-3400. /) ! Henry's law exponent used to calculate the temp dep henry's law coeff (table 19.4 Seinfeld and Pandis)
      ! H = H_T0*exp(H_exp*(1./298.-1./temp)) ! temp dependent henry's law coeff, table 19.4 (mol L⁻¹ atm⁻¹)
      ! H_dimless = H*Rg*0.009869*temp ! Dimensionless Henry's law coeff, see App A in Jacobson for unit conversion
      ! Sc_CO2 = 2073.1-125.62*19.85+3.6276*19.85**2.-0.043219*19.85**3. ! The Schmidth nr of CO2 in seawater at 293 K (19.85 degrees Celcius)
      ! Diff_water=1D-9
	  ! Diff_water(1:12) = (/ 1D-9,1D-9,1D-9,1D-9,1D-9,1D-9,1D-9,1D-9,1D-9,1D-9,1D-9,1D-9 /)! gas diffusivity in water (typical values of any gas in water, fundamentals of environmental engineering)
      ! Sc_water = 1D-6/Diff_water ! Sc nr
      ! Sc_ratio = Sc_CO2/Sc_water ! Sc nr ratio at the temp of interet
      ! ! Estimations of liquid-phase mass transfer coeff (Liss and Merlivat, 1986)
      ! IF (windspeed <= 3.6) THEN
        ! kl = (0.17*windspeed*Sc_ratio**(2./3.))/3600.*1D-2 ! (m/s)
      ! ELSE IF (windspeed > 3.6 .and. windspeed <= 13.) THEN
        ! kl = (0.612*Sc_ratio**(2./3.)+(2.85*windspeed-10.26)*Sc_ratio**0.5)/3600.*1D-2
      ! ELSE
        ! kl = (0.612*Sc_ratio**(2./3.)+(5.9*windspeed-49.9)*Sc_ratio**0.5)/3600.*1D-2
      ! END IF
      ! rc(1:12) = 1./kg+1./(kl(1:12)*H_dimless(1:12)) ! (s/m)

! ! Parameterization from Zhang et al 2003 where they set rc_ground(O3)=2000s/m and rc_ground(SO2)=20s/m and scale the other gases with coefficients from Zhang et al 2002
    ! ! alpha = 1D0
    ! ! bet = 1D0
    ! ! alpha(1:12) = (/ 1.,0.,0.,0.,10.,1.,0.8,0.8,2.,1.,10.,1./) ! table 1 in zhang et al 2002
    ! ! bet(1:12) = (/ 0.,1.,0.8,0.8,10.,1.,0.2,0.2,2.,1.,10.,0./) ! table 1 in Zhang et al 2002
    ! ! rc(1:12) = 1./(alpha(1:12)/20.+bet(1:12)/2000.) ! s/m      


    ! ELSE
     ! HENRYWIN calculations gives H_eff>1D12 M/atm for HOM, all  H_eff>1D10  gives 0 resistance.
       H_eff(1:12) = (/ 1D5,1D-6,1D-2,2D-3,1D14,1D5,6D3,220D0,1D5,1D14,1D14,2D4/) ! effective henry's lay const, table 19.4, H2SO4 and hcl assumed same as nitric acid (both strong acids)
       H_eff(13:NCOND+12)=Henrys_coeff
       f0=0D0
       f0(1:12) = (/ 0.,1.,0.1,0.,0.,1.,0.,0.3,0.1,0.,0.,0./) ! noramlized reactivity, table 19.4
       IF (temp >= 273.15 .and. temp <= 40+273.15) THEN
         rp = rj(landuse_index,season_index)*(1.+(200./(DSWF+0.1))**2.*(400./((temp-273.15)*(40.-(temp-273.15))))) ! bulk canopy stomatal resistance
       ELSE
         rp = 1D10 ! outside temp range less than 0 cel or more than 40 cel stomata is assumed to be closed
       END IF    
       rst = rp*D_ratio+1./(3.3D-4*H_eff+100.*f0) ! combined min stomatal and mesophyll resistance  
       rlui = rlu(landuse_index,season_index)*(1./(1D-5*H_eff+f0)) ! resistance of the outer surfaces in the upper canopy for each specific gas
       !IF (rain >= 0.1) THEN 
       ! rlui(1) = (1./5000.+1./(3*rlui(1)))**(-1) ! Upper canopy resistance for so2 when surface is covered by rainwater
       ! rlui(2) = (1./1000.+1./(3*rlui(1)))**(-1) ! Upper canopy resistance for o3 when surface is covered by rainwater
       ! rlui(3:10) = (1/(3*rlui(3:10))+1D-7*H_eff(3:10)+f0(3:10)/rlui(2))**(-1) ! Upper canopy resistance for rest of species when surface is covered by rainwater 
       !END IF     
       rdc = 100.*(1.+1000./(DSWF+10.)) ! resistance to transfer by buoyant convection
       rcli = (1D-5*H_eff/rclS(landuse_index,season_index)+f0/rclO(landuse_index,season_index))**(-1.)
       IF (rgsS(landuse_index,season_index) < 1D-10) THEN
        rgsi = 1D-50
       ELSE   
        rgsi = (1D-5*H_eff/rgsS(landuse_index,season_index)+f0/rgsO(landuse_index,season_index))**(-1.) ! the resistance of the exposed surfaces on the groud (soil,leaf litter, ground)
       END IF
	   
       !rc = rst!s1./rst+1./rlui+1./(rdc+rcli)+1./(rac(landuse_index,season_index)+rgsi))**(-1.) ! (s/m)
       rc = (1./rst+1./rlui+1./(rdc+rcli)+1./(rac(landuse_index,season_index)+rgsi))**(-1.) 
 !      rc(2)=rst(2) ! Only stomatal uptake of O3 in order to match observed O3 flux
!     END IF   

      ! dry deposition velocity 
      vd_gas = 1D0/(ra+rb+rc) ! [m/s]

      
      WHERE (vd_gas >= 0.1) vd_gas = 0.1

    END SUBROUTINE dry_dep_gases

    subroutine full_stationary_rebinning(N_bins,V_bins,d_p,dp_dry,MX,qX,dens, &
      c_p,c_p_backg,dt,Nconc_evap1, comp_evap)

      implicit none

      REAL(DP), dimension(nr_bins) :: vp_dry,vp_wet, N_bins_fixed
      REAL(dp), DIMENSION(NSPEC_P,nr_bins), INTENT(in) :: c_p_backg
      REAL(dp), DIMENSION(nr_bins), INTENT(inout) :: d_p, dp_dry, dens, N_bins
      REAL(dp), DIMENSION(nr_bins), INTENT(out) :: V_bins
      REAL(dp), DIMENSION(NSPEC_P,nr_bins), INTENT(inout) :: c_p
      REAL(dp), DIMENSION(NSPEC_P,nr_bins) :: c_p_fixed
      REAL(dp), DIMENSION(NSPEC_P), INTENT(in) :: MX,qX
      REAL(dp), DIMENSION(NSPEC_P) :: surf_tens
      REAL(dp), DIMENSION(NSPEC_P,nr_bins) :: S_Kelvin
      REAL(dp) :: dp_max, r1, r2

      REAL(dp), DIMENSION(nr_bins+1) :: vp_fixed
      INTEGER :: j,a
      REAL(dp), INTENT(in) :: dt


      !!!! carlton addition for clusterin
    
      REAL(DP), INTENT(OUT):: Nconc_evap1, comp_evap(:)

      comp_evap=0D0; Nconc_evap1=0D0;
      
      ! write(*,*) Nconc_evap1

      DO j = 1,nr_bins
        vp_dry(j)=SUM(c_p(index_dry,j)/Na*MX(index_dry)/qX(index_dry)/(N_bins(j)*1D-6)) ! m^3     
      END DO
      
      vp_fixed(1:nr_bins)=dp_dry**3D0*pi/6D0 ! dry single particle particle volume
      N_bins_fixed=0D0
      c_p_fixed=0D0
      dp_max = dp_dry(nr_bins)*dp_dry(nr_bins)/dp_dry(nr_bins-1)
      vp_fixed(nr_bins+1)=(pi*dp_max**3D0)/6D0

      DO j = 1,nr_bins
        a = MINLOC(vp_fixed-vp_dry(j),1,mask = (vp_fixed-vp_dry(j)) >= 0D0)
        IF (a > j) THEN ! growth
          r1 = (vp_fixed(a)-vp_dry(j))/(vp_fixed(a)-vp_fixed(a-1)) ! Fraction of particles in size bin a-1
          r2 = 1D0-r1 ! Fraction of particles in next size bin (a)
          IF (a > nr_bins+1) THEN
          ELSE IF (a == nr_bins+1) THEN
            N_bins_fixed(a-1) = N_bins_fixed(a-1)+r1*N_bins(j)
            c_p_fixed(:,a-1) = c_p_fixed(:,a-1)+vp_fixed(a-1)/vp_dry(j)*r1*c_p(:,j)
          ELSE
            N_bins_fixed(a-1) = N_bins_fixed(a-1)+r1*N_bins(j)
            N_bins_fixed(a) = N_bins_fixed(a)+r2*N_bins(j)
            c_p_fixed(:,a-1) = c_p_fixed(:,a-1)+vp_fixed(a-1)/vp_dry(j)*r1*c_p(:,j)
            c_p_fixed(:,a) = c_p_fixed(:,a)+vp_fixed(a)/vp_dry(j)*r2*c_p(:,j)
          END IF
        ELSE ! evaporation
          IF (a == 0) THEN           
          ELSEIF (a == 1) THEN
            r1 = (vp_dry(j))/(vp_fixed(1)) ! Fraction of particles in size bin a
            N_bins_fixed(1) = N_bins_fixed(1)+r1*N_bins(j)
            c_p_fixed(:,1) = c_p_fixed(:,1)+r1*c_p(:,j)
            ! c_p_fixed(:,1) = c_p_fixed(:,1)+vp_fixed(1)/vp_dry(j)*r1*c_p(:,j)
            r2 = 1D0-r1
            Nconc_evap1 = Nconc_evap1 + r2*N_bins(j)
            comp_evap = comp_evap+ r2*c_p(:,j) !molec/cm3
            ! write(*,*) 'Nconc_evap' , Nconc_evap1
              
          ELSE
            r1 = (vp_dry(j)-vp_fixed(a-1))/(vp_fixed(a)-vp_fixed(a-1)) ! Fraction of particles in size bin i
            r2 = 1D0-r1 ! Fraction of particles in previous size bin (i-1)
            N_bins_fixed(a) = N_bins_fixed(a)+r1*N_bins(j)
            N_bins_fixed(a-1) = N_bins_fixed(a-1)+r2*N_bins(j)
            c_p_fixed(:,a) = c_p_fixed(:,a)+vp_fixed(a)/vp_dry(j)*r1*c_p(:,j)
            c_p_fixed(:,a-1) = c_p_fixed(:,a-1)+vp_fixed(a-1)/vp_dry(j)*r2*c_p(:,j)
          END IF
        END IF
      END DO

      c_p=c_p_fixed
      N_bins=N_bins_fixed
 
       DO j = 1,nr_bins
      
        IF (N_bins(j)<1D-6) THEN ! This can occur with the full-stationary method if some particles grow to more than the next particle size bin within one condensation time step
          N_bins(j)=1D-3		
          c_p(:,j)=c_p_backg(:,j)*1D-3;       
          ! c_p(:,j)=c_p(:,j)*vp_fixed(j)/SUM(c_p(index_dry,j)/Na*MX(index_dry)/qX(index_dry))*N_bins(j)*1D-6   
        END IF
       
        vp_dry(j)=SUM(c_p(index_dry,j)/Na*MX(index_dry)/qX(index_dry)/(N_bins(j)*1D-6)) ! m^3     
        vp_wet(j)=SUM(c_p(:,j)/Na*MX/qX/(N_bins(j)*1D-6)) ! m^3
        V_bins(j)=N_bins(j)*vp_wet(j)
        dens(j)=SUM(c_p(:,j)*MX*1D6)/Na/V_bins(j) ! Total particle density

      END DO

      dp_dry=(vp_dry*6D0/pi)**(1D0/3D0) ! Dry particle diameters
      d_p=(vp_wet*6D0/pi)**(1D0/3D0) ! Dry particle diameters 

    end subroutine full_stationary_rebinning
	

END MODULE dynamicsDMS_atm
