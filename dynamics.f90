MODULE dynamics
    USE constants
    USE second_Precision, ONLY : dp    ! KPP Numerical type
!    USE second_Parameters, !ONLY : NSPEC, ind_H2SO4, ind_HNO3

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: dry_dep_particles, dNdlogDp, fraction_POA_marine, &
 condensation, coagulation, wet_deposition, nucleation,sea_spray, &
 dry_dep_gases, wet_deposition_gas !, cloud_activation_processing

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
        dNdlogDp_modes(ii,:) = ( N_modes(ii) / 2.506628 / LOG10(s(ii)))*&       
        EXP(-((LOG10(d_g) - log10(dm(ii)))**2 ) / (2*LOG10(s(ii))**2 ))
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

      REAL(dp), DIMENSION(Nlayers,NSPEC_P,nr_bins), INTENT(inout) :: c_p
      REAL(dp), DIMENSION(nr_bins), INTENT(inout) :: N_bins,V_bins,dens_p,d_p,dp_dry
      REAL(dp), INTENT(in) :: cH2SO4,VOC_nucl,dt,temp,RH,q_ion,cHOM,CS_air
      REAL(dp), DIMENSION(NSPEC_P), INTENT(in) :: MX,qX
      REAL(dp), DIMENSION(Nlayers,NSPEC_P), INTENT(in)  :: c_p_nucl
      REAL(dp), DIMENSION(NSPEC_P,nr_bins)  :: c_p_tot
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
        c_p(:,:,5) = c_p(:,:,5)+c_p_nucl*Jnucl_tot*dt
        c_p_tot=SUM(c_p,DIM=1) ! molec/cm^3 of each species in each size bin
        vp_dry=SUM(c_p_tot(index_dry,5)/Na*MX(index_dry)/qX(index_dry)/(N_bins(5)*1D-6)) ! m^3
        vp_wet=SUM(c_p_tot(:,5)/Na*MX/qX/(N_bins(5)*1D-6)) ! m^3
        V_bins(5)=N_bins(5)*vp_wet
        dens_p(5)=SUM(c_p_tot(:,5)*MX*1D6)/Na/V_bins(5) ! Total particle density
        dp_dry(5)=(vp_dry*6D0/pi)**(1D0/3D0) ! Dry particle diameters
        d_p(5)=(vp_wet*6D0/pi)**(1D0/3D0) ! Particle diameters
 
    END SUBROUTINE nucleation

    !-----------------------------------------------------------!
    ! Dry deposition for particles                              !
    !-----------------------------------------------------------!
    SUBROUTINE dry_dep_particles(landuse_index,d_p,dens_p,temp,press,RH,u,PBLH,SHTF,windspeed,tr,v_dep,season_index)

      INTEGER, INTENT(in) :: landuse_index
      REAL(dp), DIMENSION(nr_bins), INTENT(in) :: d_p, dens_p ! arithmetic mean diameter and dry particle densities
      REAL(dp), INTENT(in) :: temp,press,RH,u,PBLH,SHTF,windspeed ! meteorological parameters:
        ! temp [K],pressure [Pa],RH,u-resp v-comp of mom flux [N/m^2],boundary layer height [m], sensible heat net flux at surface [W/m^2], windspeed at 10 m ([m/s] (GDAS)
        ! and time along trajectory in hour and trajectory time step
      INTEGER, INTENT(in)  :: tr,season_index
      REAL(dp), DIMENSION(nr_bins), INTENT(out) :: v_dep ! deposition loss-rate coefficient for each particels size bin [m/h]
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
      REAL(dp) :: v,mu,k_b,q,Cc,D,Sc,vs,St,R1,rb,th_speed,Kn
        
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
        vs = d_p(i)**2.*(dens_p(i)-dens_air)*Cc/(9.*mu) ! sedimentation velosity, eq 20.4 (m/s)
 
        IF (landuse_index == 6 .OR. landuse_index == 8) THEN ! Slinn and Slinn 1980
!   rb = u**2./(kar*windspeed)*(Sc**(-0.5)+10.**(-3./(vs*u**2./(v*g))))
         rb = 1./(u*(Sc**(-0.5)+10.**(-3./(vs*u**2./(g*v))))) ! E.-Y Nho-Kim et al 2004, Atmos Envi
          v_dep(i) = 1D0/(ra+rb+ra*rb*vs)+vs ! over ocean
        ELSE IF (landuse_index == 7) THEN
          v_dep(i) = 1./ra+vs ! no dry dep via rb urban areas
        ELSE 
          St = vs*u/(g*r_coll(landuse_index,season_index)) ! Stokes number vegetation
          R1 = EXP(-(St)**(0.5)) ! fraction of particles, once in contact, that sticks to the surface
          rb = 1D0/(3.*u*R1*(Sc**(-j_landuse(landuse_index))+ & 
          (St/(a_landuse(landuse_index)+St))**2+0.5*(d_p(i)/r_coll(landuse_index,season_index))**2))
          v_dep(i) =1D0/(ra+rb+ra*rb*vs)+vs ! m/s
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
    SUBROUTINE wet_deposition_gas(rain,Scaveng_gas,RH)
    REAL(dp), INTENT(in)        :: rain ! rainfall intensity (mm/h)
    REAL(dp), INTENT(in)        :: RH ! %
    REAL(dp), DIMENSION(6), INTENT(out) :: Scaveng_gas
    REAL(dp), DIMENSION(6) :: Win,Wsub
    REAL(dp) :: Precipitation, hs
    Precipitation=rain/3.6D3 ! Presipitation rate in kg m^-2 s^-1
    hs=1000 ! characteristic scavenging depth (assumed to be 1000 m)
    !     SO2  HNO3 HONO NH3 H2O2 HCHO
    Win=(/0.3, 1.4, 1.4, 1.4, 1.4, 0.1/)*1D6
    Wsub=(/0.15, 0.5, 0.5, 0.5, 0.5, 0.03/)*1D6
    
    IF (RH>=98) THEN
    Scaveng_gas=Win*Precipitation/(hs*qH2O) ! s^-1
    ELSE
    Scaveng_gas=Wsub*Precipitation/(hs*qH2O) ! s^-1
    END IF
    
    END SUBROUTINE wet_deposition_gas
    
    ! IN cloud SO2 and H2O2 dissolution and S(VI) formation:
    
    SUBROUTINE cloud_activation_processing(c_p,N_bins,SO2,H2O2,T,S_c,S_super,I_nSVI,LWC)
    REAL(dp), DIMENSION(Nlayers,NSPEC_P,nr_bins), INTENT(in) :: c_p
    REAL(dp), DIMENSION(nr_bins), INTENT(in) :: N_bins
    REAL(dp), DIMENSION(nr_bins), INTENT(out) :: S_c
    REAL(dp), INTENT(out) :: S_super, I_nSVI
    REAL(dp), INTENT(in) :: T,LWC
    REAL(dp), DIMENSION(nr_bins) :: ns, B_Kohler, vp0,dp0, alpha1, alpha2
    REAL(dp), INTENT(inout) :: SO2, H2O2
    REAL(dp) :: CCN, dp_drops, A_Kohler, surf_tens, m_drops, cH
    REAL(dp) :: Ks1,Ks2,H_SO2,H_H2O2,H_S,pSO2,pH2O2,cSO2aq,cH2O2aq,kSIV,K_SIV,xHSO3
   
    
    !LWC = 5D-4 ! Assumed low level cloud liquid water content (kg/m^3)
    S_super = 2D-3 ! Assumed maximum supersaturation in the clouds (fraction)
    cH = 1D-4 ! Assumed cloud droplet acidity (mol H+ dm^-3, pH = 4)
    
    surf_tens=(76.1-0.155*(T-273.15))*1D-3 ! Surface tension of water (kg s^-2) 
    
    ns=(SUM(SUM(c_p(:,1:5,:), DIM=1),DIM=1)+SUM(SUM(c_p(:,9:NSPEC_P,:),DIM=1),DIM=1))/(N_bins*1D-6)/Na ! Moles of soluble material in each single particle inorganic ions + all organic compounds but not soot
    
    vp0=SUM(c_p(:,6,:)/Na*MEC/qEC)/(N_bins*1D-6) ! Equivalent volume of insoluble core (considered to be soot only)
    
    dp0=(vp0*6D0/pi)**(1D0/3D0) ! Equivalent diameter of insoluble core (considered to be soot only)
     
    A_Kohler=4*MH2O*surf_tens/(Rg*T*qH2O)
    B_Kohler=6*ns*MH2O/(pi*qH2O)
    
    ! Calculate the critical supersaturation for droplet activation according to Kokkola et al., 2008 
    alpha1=12*(81D0*dp0**6D0+12D0*dp0**3D0*(3D0*B_Kohler/A_Kohler)**(3D0/2D0))**(1D0/2D0)
    alpha2=(108D0*dp0**3D0+8*(3D0*B_Kohler/A_Kohler)**(3D0/2D0)+alpha1)**(1D0/3D0)
    S_c= EXP( A_Kohler/(alpha2/6D0 + 2D0/3D0*(3D0*B_Kohler/A_Kohler)*1D0/alpha2+&
    1D0/3D0*(3D0*B_Kohler/A_Kohler)**(1D0/2D0))+&
    B_Kohler/((alpha2/6D0 + 2D0/3D0*(3D0*B_Kohler/A_Kohler)*1D0/alpha2+&
    1D0/3D0*(3D0*B_Kohler/A_Kohler)**(1D0/2D0))**3D0 - dp0**3D0))-1D0
  
    CCN=SUM(PACK(N_bins, S_c <= S_super)) ! Number of CCN (activated droplets at 0.2 % supersaturation)
    m_drops=LWC/CCN ! Average single droplet mass/volume of cloud droplets (kg or dm^3)
    
    ! saturation vapor pressure of SO2 and H2O2 above the droplet surfaces
    Ks1=1.3D-2*EXP(1960D0*(1D0/T-1D0/298D0)) ! mol L^-1 water  
    Ks2=6.6D-8*EXP(1500D0*(1D0/T-1D0/298D0)) ! mol  L^-1 water 	
    H_SO2=1.23D0*EXP(-6.25D0/Rg*(1D0/298D0+1D0/T))  ! Henry's law coefficient SO2(g) mol  L^-1 water atm^-1
    H_H2O2=1D5*EXP(-14.5D0/Rg*(1D0/298D0+1D0/T)) ! Henry's law coefficient H2O2(g) mol  L^-1 water atm^-1 at 298 K.
    H_S=H_SO2*(1D0+Ks1/cH+Ks1*Ks2/(cH**2D0)) ! effective Henry's law constant SO2(g)
	
	!Phase transfer:
    !H_DMS  = 0.56*EXP(4480D0*(1D0/T-1D0/298D0)); ! Henry's law coefficient DMS(g) mol  L^-1 water atm^-1
    !H_DMSO = 1D7*EXP(2580D0*(1D0/T-1D0/298D0));  ! Henry's law coefficient DMSO(g) mol  L^-1 water atm^-1
    !H_MSIA  = 1D8*EXP(1760D0*(1D0/T-1D0/298D0)); ! Henry's law coefficient MSIOA(g) mol  L^-1 water atm^-1
    !H_O3 = 0.011*EXP(2400D0*(1D0/T-1D0/298D0));  ! Henry's law coefficient O3(g) mol  L^-1 water atm^-1
    !H_OH = 3D1*EXP(4500D0*(1D0/TEMP-1D0/298D0)); ! Henry's law coefficient OH(g) mol  L^-1 water atm^-1
     
    pSO2=(1D0/101325D0)*1D6*SO2*Rg*T/Na    ! Vapour pressure of SO2 in atm
    pH2O2=(1D0/101325D0)*1D6*H2O2*Rg*T/Na  ! Vapour pressure of H2O2 in atm
    !pDMS=(1D0/101325D0)*1D6*DMS*Rg*T/Na   ! Vapour pressure of DMS in atm
	!pDMSO=(1D0/101325D0)*1D6*DMSO*Rg*T/Na ! Vapour pressure of DMSO in atm
	!pMSIA=(1D0/101325D0)*1D6*MSIA*Rg*T/Na ! Vapour pressure of DMSO in atm
	!pOH=(1D0/101325D0)*1D6*OH*Rg*T/Na    ! Vapour pressure of OH in atm
	!pO3=(1D0/101325D0)*1D6*O3*Rg*T/Na ! Vapour pressure of O3 in atm
	
	! Aqueous phase equilibrium concentrations:
    cSO2aq=pSO2*H_S ! moles L^-1
    cH2O2aq=pH2O2*H_H2O2 ! moles L^-1
	!cDMSaq=pDMS*H_DMS ! moles L^-1
    !cDMSOaq=pDMSO*H_DMSO ! moles L^-1
    !cMSIAaq=pMSIA*H_MSIA ! moles L^-1
	!cOHaq=pOH*H_OH ! moles L^-1
	!cO3aq=pO3*H_O3 ! moles L^-1
   
   ! Table S13 in Hoffman et al., PNAS 2016: 
   ! DMS, DMSO, MSIA and MSA oxidation in the aquesous phase is driven by O3 and OH oxidation!!!!
   !k_DMS_O3 = 8.61D8*EXP(-2600*(1D0/TEMP-1D0/298D0))  ! DMS + O3 -> DMSO L^3 mol^-1 s^-1 
   !k_DMSO_OH = 6.65D9*EXP(-1270*(1D0/TEMP-1D0/298D0)) ! DMSO + OH -> MSIA L^3 mol^-1 s^-1
   !k_MSIA_O3 = 3.5D7 ! MSIA + O3 -> MSA  ! L^3 mol^-1 s^-1
   !k_MSIA_OH = 6D9 ! MSIA + OH -> MSA  L^3 mol^-1 s^-1
   

    !dDMSdt=-k_DMS_O3*cO3aq*cDMSaq*m_drops ! DMS(aq) loss rate mol s^-1 in each single droplet
    !dDMSOdt=k_DMS_O3*cO3aq*cDMSaq*m_drops-k_DMSO_OH*cOHaq*cDMSOaq*m_drops ! DMSO(aq) formation rate mol s^-1 in each single droplet
	!dMSIAdt=k_DMSO_OH*cOHaq*cDMSOaq*m_drops ! MSIA(aq) formation rate mol s^-1 in each single droplet
	!dMSAdt=k_MSIA_OH*cOHaq*cMSIAaq*m_drops +k_MSIA_O3*cO3aq*cMSIAaq*m_drops ! MSA(aq) formation rate mol s^-1 in each single droplet
	!dOHdt=-k_MSIA_OH*cOHaq*cMSIAaq*m_drops-k_DMSO_OH*cOHaq*cDMSOaq*m_drops ! OH(aq) loss rate mol s^-1 in each single droplet
	!dO3dt=-k_MSIA_O3*cO3aq*cMSIAaq*m_drops-k_DMS_O3*cO3aq*cDMSaq*m_drops ! O3(aq) loss rate mol s^-1 in each single droplet
	
    kSIV=7.45D7*EXP(-15.96D0*(298D0/T-1D0)) ! H2O2+S(IV)-> S(VI), Jacobson, 2005 table B.8
    K_SIV=13D0  ! mol  L^-1 water atm^-1
    xHSO3=(1D0+cH/Ks1+Ks2/cH)**(-1D0) ! molar fraction HSO3- of total S(IV)    I_nSVI=dp*0;
    I_nSVI=m_drops*cH2O2aq*xHSO3*cSO2aq*cH*kSIV/(1D0+K_SIV*cH) ! S(VI) formation rate mol/s
    
    ! Update the particle and gas-phase chemical composition after cloud processing and heterogeneous S(VI) formation
    !where (S_c<=S_super) c_p(Nlayers,1,:)=c_p(Nlayers,1,:)+I_nSVI*Na*N_bins*1D-6*dt
   
    !SO2=MAXVAL((/SO2-SUM(I_nSVI*Na*PACK(N_bins, S_c <= S_super)*1D-6)*dt,1D-10/))
    !H2O2=MAXVAL((/H2O2-SUM(I_nSVI*Na*PACK(N_bins, S_c <= S_super)*1D-6)*dt,1D-10/))
    
	
    END SUBROUTINE cloud_activation_processing
    
    
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

    !--------------------------------------------------------------!
    ! Sea-spray particle emission                                  !
    !--------------------------------------------------------------!
    SUBROUTINE sea_spray(wind10,temp,dlogdp,dp_dry,E_sea_salt)
      REAL(dp), INTENT(in) :: wind10,temp
      REAL(dp), DIMENSION(nr_bins), INTENT(in) :: dlogdp,dp_dry
      REAL(dp), DIMENSION(nr_bins), INTENT(out) :: E_sea_salt

      REAL(dp) :: A_whitecape,T_water,cc4,cc3,cc2,cc1,cc0,dd4,dd3,dd2,dd1,dd0, &
             ccc4,ccc3,ccc2,ccc1,ccc0,ddd4,ddd3,ddd2,ddd1,ddd0, &
             cccc4,cccc3,cccc2,cccc1,cccc0,dddd4,dddd3,dddd2,dddd1,dddd0
      REAL(dp), DIMENSION(nr_bins) :: Ak,Bk,dF_saltdlogDp
      INTEGER :: i

      A_whitecape = 0.01*3.84D-4*wind10**3.41 ! estimated whitecap area fraction
      ! Estimate sea surface temperature based on temperature at 2 m
      IF (temp > 283.15) THEN
 T_water = 283.15
      ELSE IF (temp > 273.15) THEN
 T_water = 273.15 + 5.
      ELSE
 T_water = 273.15
      END IF

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

    END SUBROUTINE sea_spray              

    !------------------------------------------------------!
    ! Coagulation                                          !
    !------------------------------------------------------!
    SUBROUTINE coagulation(N_bins,V_bins,c_p,d_p,dp_dry,dt,T,p,dens_p,MX,qX)
   
      REAL(dp), DIMENSION(nr_bins), INTENT(inout) :: N_bins
      REAL(dp), DIMENSION(Nlayers,NSPEC_P,nr_bins), INTENT(inout) :: c_p
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
      REAL(dp), DIMENSION(Nlayers,NSPEC_P,nr_bins) :: c_p_single
      REAL(dp), DIMENSION(Nlayers,NSPEC_P,nr_bins+1) :: dcpdt
      REAL(dp), DIMENSION(Nlayers,NSPEC_P) :: c_p_coag
      REAL(dp), DIMENSION(NSPEC_P,nr_bins) :: c_p_tot
 
      INTEGER :: i,j,m,n

      dyn_visc = 1.8D-5*(T/298.)**0.85  ! dynamic viscosity of air
      l_gas=2D0*dyn_visc/(p*SQRT(8D0*Mair/(pi*Rg*T))) ! Gas mean free path in air (m)

      DO j=1,nr_bins
        c_p_single(:,:,j) = c_p(:,:,j)/(N_bins(j)*1D-6) ! concentration of each condensable compound in particle phase in a single particle in each bin
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
          c_p_coag = c_p_single(:,:,m)+c_p_single(:,:,j) ! compounds in formed single particle (molecues/#) 
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
              dcpdt(:,:,i) = dcpdt(:,:,i)+Vp(i)/Vp_coag*r1*Coag_source*c_p_coag ! (molec/(m^3 s))
              dcpdt(:,:,i+1) = dcpdt(:,:,i+1)+Vp(i+1)/Vp_coag*r2*Coag_source*c_p_coag ! (molec/(m^3 s))
            END IF
   END DO
 END DO
      END DO

      ! Coagulation sink included:
      DO j = 1,nr_bins
        Coag_sink = K(j,:)*N_bins*N_bins(j) ! (m^3/s)
        Coag_sinktot = sum(Coag_sink)
        dNdt(j) = dNdt(j)-Coag_sinktot
        dcpdt(:,:,j) = dcpdt(:,:,j)-Coag_sinktot*c_p_single(:,:,j) ! (molec/s)
      END DO


      N_bins = N_bins+dNdt(1:nr_bins)*dt ! New particle concentration in each size bin
      c_p = c_p+dcpdt(:,:,1:nr_bins)*dt*1D-6 ! molec/ cm^3

      c_p_tot(1:NSPEC_P,1:nr_bins)=sum(c_p,DIM=1)

      DO j = 1,nr_bins
        vp_dry(j)=SUM(c_p_tot(index_dry,j)/Na*MX(index_dry)/qX(index_dry)/(N_bins(j)*1D-6)) ! m^3   
        vp_wet(j)=SUM(c_p_tot(:,j)/Na*MX/qX/(N_bins(j)*1D-6)) ! m^3
        V_bins(j)=N_bins(j)*vp_wet(j)
        dens_p(j)=SUM(c_p_tot(:,j)*MX*1D6)/Na/V_bins(j) ! Total particle density
      END DO
      dp_dry=(vp_dry*6D0/pi)**(1D0/3D0) ! Dry particle diameters
      d_p=(vp_wet*6D0/pi)**(1D0/3D0) ! Dry particle diameters
    END SUBROUTINE coagulation

    !--------------------------------------------------------------!
    ! Condensation                                                 !
    !--------------------------------------------------------------!
    SUBROUTINE condensation(N_bins,V_bins,p,T,RH,d_p,dp_dry,aX,MX,qX,molec_radius,dens,corg,cH2SO4,cHNO3,&
     cHCl,cNH3,cSO2,cH2O2,c_p,W,Kprim_HNO3,Kprim_HCl,Hprim_NH3,Kprim_NH3,fHSO4,fSO4,fNO3,fCl,&
     mHCO3,mCO3,mOH,mCOO,CHNO3s,CHCls,psat_org,Dorg,dt,tr,yorg,c_p_backg,indexALDEHY,indexROOH,CS_H2SO4,&
     dimerC,dimerO,dimerN,dimerH,NrC,NrO,NrN,NrH,CS_air,LWC)    

      REAL(dp), DIMENSION(nr_bins), INTENT(in) :: W,Kprim_HNO3,Kprim_HCl,Hprim_NH3,&
       Kprim_NH3,fHSO4,fSO4,fNO3,fCl,mHCO3,mCO3,mOH,mCOO,CHNO3s,CHCls
      REAL(dp), DIMENSION(NCOND+1,nr_bins), INTENT(in) :: yorg
      REAL(dp), DIMENSION(Nlayers,NSPEC_P,nr_bins), INTENT(in) :: c_p_backg
      REAL(dp), DIMENSION(nr_bins), INTENT(inout) :: d_p, dp_dry, dens, N_bins, &
      dimerC, dimerO, dimerN, dimerH
      REAL(dp), DIMENSION(nr_bins), INTENT(out) :: V_bins
      REAL(dp), DIMENSION(Nlayers,NSPEC_P,nr_bins), INTENT(inout) :: c_p
      REAL(dp), DIMENSION(Nlayers,NSPEC_P,nr_bins) :: c_p_fixed,c_p_old
      REAL(dp), INTENT(in) :: T,p,RH,LWC
      REAL(dp), INTENT(inout) :: cHNO3, cNH3, cHCl,cSO2,cH2O2
      REAL(dp), INTENT(in) :: cH2SO4
      REAL(dp), INTENT(out)   :: CS_H2SO4, CS_air
      REAL(dp), DIMENSION(NCOND), INTENT(inout) :: corg
      REAL(dp), DIMENSION(NCOND), INTENT(in)    :: psat_org,indexALDEHY, indexROOH, &
      NrC,NrO,NrN,NrH
      REAL(dp), DIMENSION(NSPEC_P), INTENT(in) :: MX,qX,aX,molec_radius
      REAL(dp), DIMENSION(NSPEC_P) :: surf_tens, dX
      REAL(dp), DIMENSION(NSPEC_P,nr_bins) :: c_p_surf,S_Kelvin, c_p_tot
      REAL(dp) :: aH2SO4,aHNO3,aHCl,CH2SO4old,CHNO3old,CHClold,CNH3old,DHNO3,DHCl,DH2SO40,Dair,&
       Keq1,Keq2,DH2SO4,mH2SO4,mHNO3,mHCl,m_air,speedH2SO4,speedHNO3,speedHCl,speedair,CH2SO4sat,CH2SO4m,&
       CHNO3m,CHClm,CNH3m,CH2SO4tot,cNH3tot,CHNO3tot,CHCltot,Corgm,Corgtot, errorCNH3, fn,&
       fprimn, dp_max, r1, r2, l_gas, dyn_visc, d_H2SO4, d_air,dt, k_PHAolig, k_HOMfrag, mROOH, mALDEHY,&
       S_super,I_nSVI
      REAL(dp), DIMENSION(NCOND,nr_bins) :: corgp,corgpold,korg,corgp_eq
      REAL(dp), DIMENSION(NCOND) :: Corgold,speedorg,aorg,morg,m_org,dm_orgdt,m_orgold,nrOOH,&
      dm_orgdt_frag, dm_orgdt_frag2
      REAL(dp), DIMENSION(NCOND), INTENT(in) :: Dorg

      REAL(dp), DIMENSION(nr_bins) :: cH2SO4pold,cHNO3pold,cHClpold,cNH3pold,KnH2SO4,&
       KnHNO3,KnHCl,Knair,f_corH2SO4,f_corHNO3,f_corHCl,f_cororg,f_corair,DH2SO4eff,DHNO3eff,DHCleff,Daireff,kH2SO4,&
       kHNO3,kHCl,kair,S_KelvinH2SO4,S_KelvinHNO3,S_KelvinHCl,cH2SO4p,cHNO3p,cHClp, term1, term2, &
       cNap, cCO3p, cHCO3p, cOHp, cNH3p, charge, cNH4ion, cNH3aq,vp_dry,vp_wet, N_bins_fixed, &
       Cunningh, Diff_p,m_p,speed_p, gasmeanfpH2SO4, gasmeanfpHNO3, gasmeanfpHCl,gasmeanfporg, &
       gasmeanfpair,Knorg, Dorgeff,Corgsat,xorg,N_bins_old,S_c

      REAL(dp), DIMENSION(nr_bins+1) :: vp_fixed
      INTEGER :: j,jj,a
      INTEGER, INTENT(in) :: tr
 
      N_bins_old=N_bins
      c_p_old=c_p

      dX=molec_radius*2 ! Molecule diameter (m)

      dyn_visc=1.8D-5*(T/298.)**0.85                ! dynamic viscosity of air

      c_p_surf=c_p(Nlayers,:,:) ! molec cm^-3 in particle surface layers

      surf_tens(1:8)=(76.1-0.155*(T-273.15))*1D-3 ! (kg s^-2; % inorganics
      surf_tens(9:NSPEC_P)=0.05 ! organics

      DO j=1,NSPEC_P
        S_Kelvin(j,:)=1D0+2D0*surf_tens(j)*MX(j)/(d_p/2D0*dens*Rg*T) ! equilibrium saturation ratio of condensing gas (Kelvin effect)
      END DO

      l_gas=2D0*dyn_visc/(p*SQRT(8D0*Mair/(pi*Rg*T))) ! Gas mean free path in air (m)
      Cunningh=1D0+(2D0*l_gas/d_p)*(1.257+0.4*EXP(-1.1/(2D0*l_gas/d_p)))  ! Cunninghams correction factor
      Diff_p=Cunningh*kb*T/(3D0*pi*dyn_visc*d_p)                          ! Diffusivitys for the different particle sizes m^2/s
      m_p=(dens*pi*d_p**3D0)/6D0  ! mass of particles 
      speed_p=SQRT(8D0*kb*T/(pi*m_p)) ! speed of particles

      ! Condensation inorganics:
      cH2SO4pold=c_p_surf(1,:)/Na*1D6 ! mol/m^3  
      cHNO3pold=c_p_surf(2,:)/Na*1D6  ! mol/m^3  
      cHClpold=c_p_surf(3,:)/Na*1D6   ! mol/m^3  
      cNH3pold=c_p_surf(4,:)/Na*1D6   ! mol/m^3  
      cNap=c_p_surf(5,:)/Na*1D6       ! mol/m^3 Na+

      DO j=1,NCOND
        corgpold(j,:)=c_p_surf(8+j,:)/Na*1D6 ! mol/m^3 cond org comp
      END DO
      corgp=corgpold

      aH2SO4=aX(1);aHNO3=aX(2);aHCl=aX(3); aorg=aX(9:NSPEC_P) ! Mass accommodation coefficients of condensable compounds

      CH2SO4old=cH2SO4/Na*1D6 ! vapor mole concentration (mol m^-3)
      CHNO3old=cHNO3/Na*1D6   ! vapor mole concentration (mol m^-3)
      CHClold=cHCl/Na*1D6     ! vapor mole concentration (mol m^-3)
      CNH3old=cNH3/Na*1D6     ! vapor mole concentration (mol m^-3)
      Corgold=corg/Na*1D6      

      DHNO3=(0.22/1.87)*1D-4 ! gas diffusivity HNO3 m^2/s
      DHCl=(0.22/1.42)*1D-4  ! gas diffusivity HCl m^2/s
      Dair=0.22D-4           ! Air molecule gas diffusivity m^2/s

      ! RH dependent Diffusion, diameter and mass of H2SO4 molecules:
      DH2SO40=0.094D-4       ! gas diffusivity H2SO4 m^2/s RH=0%
      Keq1=0.13D0
      Keq2=0.016D0
      DH2SO4 = (DH2SO40+0.85D0*DH2SO40*Keq1*RH+0.76D0*DH2SO40*Keq1*Keq2*(RH)**2D0)&
      /(1D0+Keq1*RH+Keq1*Keq2*(RH)**2D0) ! Diffusivity H2SO4 at ambient RH

      d_H2SO4=(((MX(1)/qX(1)/Na)*6D0/pi)**(1D0/3D0)+&
       Keq1*RH*(((MX(1)+18D-3)/qX(1)/Na)*6D0/pi)**(1D0/3D0)+&
       Keq1*Keq2*(RH)**2D0*(((MX(1)+36D-3)/qX(1)/Na)*6D0/pi)**(1D0/3D0))/&
       (1D0+Keq1*RH+Keq1*Keq2*(RH)**2D0) ! RH dependent H2SO4 diameter
       d_air=3D-10 ! Estimated molecular diameter of an average air molecule

      mH2SO4=(MX(1)/Na+Keq1*RH*(MX(1)+18D-3)/Na+Keq1*Keq2*(RH)**2D0*(MX(1)+36D-3)/Na)/&
      (1D0+Keq1*RH+Keq1*Keq2*(RH)**2D0) ! RH H2SO4 molecular mass

      mHNO3=MX(2)/Na; mHCl=MX(3)/Na
      morg=MX(9:NSPEC_P)/Na
      m_air=Mair/Na

      speedH2SO4=SQRT(8D0*kb*T/(pi*mH2SO4)) ! speed of H2SO4 molecule
      speedHNO3=SQRT(8D0*kb*T/(pi*mHNO3))   ! speed of HNO3 molecules
      speedHCl=SQRT(8D0*kb*T/(pi*mHCl))     ! speed of HCl molecules
      speedorg=SQRT(8D0*kb*T/(pi*morg))     ! speed of organic molecules
      speedair=SQRT(8D0*kb*T/(pi*m_air))     ! speed of air molecule
      

      ! Gas mean free path, Knudsen number (Lehtinen and Kulmala, 2003):
      gasmeanfpH2SO4=3D0*(DH2SO4+Diff_p)/SQRT(speedH2SO4**2D0+speed_p**2D0) 
      gasmeanfpHNO3=3D0*(DHNO3+Diff_p)/SQRT(speedHNO3**2D0+speed_p**2D0)    
      gasmeanfpHCl=3D0*(DHCl+Diff_p)/SQRT(speedHCl**2D0+speed_p**2D0)
      gasmeanfpair=3D0*(Dair+Diff_p)/SQRT(speedair**2D0+speed_p**2D0)        

      KnH2SO4=2D0*gasmeanfpH2SO4/(d_p+d_H2SO4) ! Knudsen number H2SO4
      KnHNO3=2D0*gasmeanfpHNO3/(d_p+dX(2))     ! Knudsen number HNO3
      KnHCl=2D0*gasmeanfpHCl/(d_p+dX(3))       ! Knudsen number HCl
      Knair=2D0*gasmeanfpair/(d_p+d_air)       ! Knudsen number H2SO4
      

      f_corH2SO4=(0.75*aH2SO4*(1D0+KnH2SO4))&
      /(KnH2SO4**2D0+KnH2SO4+0.283*KnH2SO4*aH2SO4+0.75*aH2SO4) ! Fuchs-Sutugin correction factor for transit and kinetic regime a=accomodation coefficient

      f_corHNO3=(0.75*aHNO3*(1D0+KnHNO3))&
      /(KnHNO3**2D0+KnHNO3+0.283*KnHNO3*aHNO3+0.75*aHNO3)  ! Fuchs-Sutugin correction factor for transit and kinetic regime

      f_corHCl=(0.75*aHCl*(1D0+KnHCl))/&
      (KnHCl**2D0+KnHCl+0.283*KnHCl*aHCl+0.75*aHCl) ! Fuchs-Sutugin correction factor for transit and kinetic regime
      
      f_corair=(0.75*1D0*(1D0+Knair))&
      /(Knair**2D0+Knair+0.283*Knair*1D0+0.75*1D0) ! Fuchs-Sutugin correction factor for transit and kinetic regime


      Daireff=(Dair+Diff_p)*f_corair                         ! m^2/s
      DH2SO4eff=(DH2SO4+Diff_p)*f_corH2SO4                   ! m^2/s
      DHNO3eff=(DHNO3+Diff_p)*f_corHNO3                      ! m^2/s
      DHCleff=(DHCl+Diff_p)*f_corHCl                         ! m^2/s

      CH2SO4sat=0D0 ! H2SO4 saturation mole concentration (mol m^-3) p/(R*T)

      kH2SO4=N_bins*2D0*pi*(d_p+d_H2SO4)*DH2SO4eff  ! mass transfer coefficient s^-1
      kHNO3=N_bins*2D0*pi*(d_p+dX(2))*DHNO3eff     ! mass transfer coefficient s^-1
      kHCl=N_bins*2D0*pi*(d_p+dX(3))*DHCleff        ! mass transfer coefficient s^-1
      kair=N_bins*2D0*pi*(d_p+d_air)*Daireff        ! mass transfer coefficient s^-1

      CS_H2SO4=SUM(kH2SO4)                          ! H2SO4 condensation sink s^-1
      CS_air=SUM(kair)                              ! air ion condensation sink s^-1


      DO j=1,NCOND
        gasmeanfporg=3D0*(Dorg(j)+Diff_p)/SQRT(speedorg(j)**2D0+speed_p**2D0) 
        Knorg=2D0*gasmeanfporg/(d_p+dX(8+j))       ! Knudsen number organic comp
        f_cororg=(0.75*aorg(j)*(1D0+Knorg))/&
        (Knorg**2D0+Knorg+0.283*Knorg*aorg(j)+0.75*aorg(j))     ! Fuchs-Sutugin correction factor for transit and kinetic regime
        Dorgeff=(Dorg(j)+Diff_p)*f_cororg                    ! m^2/s
        korg(j,:)=N_bins*2D0*pi*(d_p+dX(8+j))*Dorgeff        ! mass transfer coefficient s^-1
      END DO
      
! Condensation and possible heterogeneous reactions of organic compounds:

k_PHAolig=2D-4!200D0 !12D0 !200D0 ! Peroxyhemiacetal formation rate (kg/mole s^-1)
k_HOMfrag=1D-4 !6D-5 ! 2.6E-4 s^-1 (Krapf et al., 2016)
! The molality of ROOH and ALDEHY can become ~5 each and thus dm_orgdt=k_PHAolig*10

 dm_orgdt=0D0
 !x_SVI=c_p_surf(1,:)/(SUM(c_p_surf(9:NSPEC_P-1,:),DIM=1)+c_p_surf(1,:)) ! Mole fraction S(VI) in each particle surface layer
 !where (x_SVI<1D-3) x_SVI=1D-3 ! Uncatalyzed limit of peroxyhemiacetyle formation
 !where (x_SVI>1D0) x_SVI=1D0 ! Uncatalyzed limit of peroxyhemiacetyle formation

IF (index_oligomerization==1) THEN
!  x_SVI=1D0
  ! Peroxyhemiacetal formation involving HOM with multiple -OOH groups:  
  nrOOH = 0D0 
 nrOOH(635:742)=(/0,0,1,1,2,2,3,3,4,4,5,0,0,1,1,2,2,3,3,4,4,5,0,0,1,1,2,2,3,3,4,4,5,6,0,0,1,1,2,2,&
                 3,3,4,4,5,0,0,0,1,1,2,2,3,3,4,4,5,5,0,0,1,1,2,2,3,3,4,4,5,1,1,1,1,2,2,1,1,1,1,2,&
                 2,1,1,1,1,2,3,3,1,1,1,1,2,2,1,1,1,1,2,2,2,3,1,1,1,1,1,2/) !number of hydroperoxide / peroxy acid functional groups in the molecules
                       
  nrOOH(NCOND-5)=3  ! OLIGOMER1
  nrOOH(NCOND-4)=2; ! OLIGOMER2
  nrOOH(NCOND-3)=1; ! OLIGOMER3
  nrOOH(NCOND-2)=0; ! OLIGOMER4
  
    
  DO j=1,nr_bins
  !mROOH=SUM(corgpold(:,j)*indexROOH)/(SUM(corgpold(:,j)*MX(9:NSPEC_P))+1D-100) ! Molality of ROOH in each particle  
  !mALDEHY=SUM(corgpold(:,j)*indexALDEHY)/(SUM(corgpold(:,j)*MX(9:NSPEC_P))+1D-100) ! Molality of R=O in each particle
  !mHOM=SUM(corgpold(635:703,j))/(SUM(corgpold(:,j)*MX(9:NSPEC_P))+1D-100) ! Molality of HOM in each particle
 
  
  m_orgold=corgpold(:,j)/(SUM(corgpold(:,j)*MX(9:NSPEC_P))+1D-100) ! Molality of each organic compound in the organic + water particle phase
  m_org=m_orgold

  dm_orgdt_frag=0D0 ! Fragmentation of HOM monomers
  dm_orgdt_frag2=0D0 ! Fragmentation of HOM dimers
  where (nrOOH>1D-3) dm_orgdt_frag=m_orgold*k_HOMfrag ! First order fragmentation of HOM monomers to form products with 2, 4, 6, and 8 carbon atoms
  dm_orgdt_frag2(680:703)=dm_orgdt_frag(680:703) ! Fragmentation of HOM dimer reaction products forming C10 products
  dm_orgdt_frag2(729:742)=dm_orgdt_frag(729:742) ! Fragmentation of HOM dimer reaction products forming C10 products
  dm_orgdt_frag(680:703)=0D0
  dm_orgdt_frag(729:742)=0D0 
  dm_orgdt_frag2(NCOND-5:NCOND-2)=dm_orgdt_frag(NCOND-5:NCOND-2) ! Fragmentation of C20-Peroxyhemiacetals forming C10 products
  dm_orgdt_frag(NCOND-5:NCOND-2)=0D0 ! No fragmentation of heterogeneous oligomer reaction products
  
  mROOH=SUM(corgpold(:,j)*nrOOH)/(SUM(corgpold(:,j)*MX(9:NSPEC_P))+1D-100) ! Molality of -O-OH in each particle
  
  
  dm_orgdt=k_PHAolig*m_orgold*mROOH ! org
  dm_orgdt(NCOND-1:NCOND)=0D0

  ! dm_orgdt=k_PHAolig*m_orgold*indexALDEHY*mROOH+k_PHAolig*m_orgold*indexROOH*mALDEHY
  
  !dm_orgdt=k_PHAolig*m_orgold*mHOM ! org
  
  !dm_orgdt(663)=0D0
  
  !DO jj=1,NCOND-26
  !m_org(jj,j)=corgpold(jj,j)/(SUM(corgpold(:,j)*MX(9:NSPEC_P))+1D-100) ! Molality of each organic compound in the organic + water particle phase
  !m_orgold(jj,j)=m_org(jj,j)
  !IF (jj<NCOND-1) THEN
  !dm_orgdt(jj,j)=x_SVI(j)*k_PHAolig*m_org(jj,j)*indexALDEHY(jj)*mROOH+x_SVI(j)*k_PHAolig*m_org(jj,j)*indexROOH(jj)*mALDEHY
  !END IF
  !END DO
  
  where (m_orgold(1:NCOND-2)<1.001D0*dm_orgdt(1:NCOND-2)*dt) dm_orgdt(1:NCOND-2)=m_orgold(1:NCOND-2)/(dt*1.001D0)
  
  ! Update the moles of C, O, N and H atoms in the particles:
  dimerC(j)=dimerC(j)+SUM(NrC(1:NCOND-6)*(dm_orgdt(1:NCOND-6)+dm_orgdt_frag(1:NCOND-6)+dm_orgdt_frag2(1:NCOND-6)))*dt
  dimerO(j)=dimerO(j)+SUM(NrO(1:NCOND-6)*(dm_orgdt(1:NCOND-6)+dm_orgdt_frag(1:NCOND-6)+dm_orgdt_frag2(1:NCOND-6)))*dt
  dimerN(j)=dimerN(j)+SUM(NrN(1:NCOND-6)*(dm_orgdt(1:NCOND-6)+dm_orgdt_frag(1:NCOND-6)+dm_orgdt_frag2(1:NCOND-6)))*dt
  dimerH(j)=dimerH(j)+SUM(NrH(1:NCOND-6)*(dm_orgdt(1:NCOND-6)+dm_orgdt_frag(1:NCOND-6)+dm_orgdt_frag2(1:NCOND-6)))*dt
  
  ! Consider changes in the concentration of MCM compounds + HOM due to heterogeneous oligomerization:
  !m_org(1:NCOND-2)=m_orgold(1:NCOND-2)-dm_orgdt(1:NCOND-2)*dt ! Only consider dimer SOA for HOM x other compounds or HOM x HOM 
  !m_org(NCOND-1)=m_orgold(NCOND-1)+SUM(dm_orgdt(1:NCOND-2)*MX(9:NSPEC_P-2)/MX(NSPEC_P-1))*dt
  
  !!m_org(1:NCOND-2,j)=m_org(1:NCOND-2,j)+1D-9*dm_orgdt(1:NCOND-2,j)/(SUM(dm_orgdt(1:NCOND-2,j),DIM=1)+1D-100)*m_org(NCOND-1,j)*MX(NSPEC_P-1)/MX(9:NSPEC_P-2)*dt_cond ! Degradation of dimers back to monomers
  !!m_org(NCOND-1,j)=m_org(NCOND-1,j)-1D-9*m_org(NCOND-1,j)*dt_cond ! Degradation of dimers back to monomers
  
! Peroxyhemiacetal formation involving HOM with multiple -OOH groups:
 DO jj=1,NCOND-2 
 IF (nrOOH(jj)==4) THEN
 m_org(NCOND-5)=m_org(NCOND-5)+dm_orgdt(jj)*MX(jj+8)/MX(NSPEC_P-5)*dt ! Oligomer SOA of compounds still contaning 3 -O-O-H groups
 m_org(jj) = m_org(jj)-dm_orgdt(jj)*dt
 ELSE IF (nrOOH(jj)==3) THEN
 m_org(NCOND-4)=m_org(NCOND-4)+dm_orgdt(jj)*MX(jj+8)/MX(NSPEC_P-4)*dt ! Oligomer SOA of compounds still contaning 2 -O-O-H groups
 m_org(jj) = m_org(jj)-dm_orgdt(jj)*dt
 ELSE IF (nrOOH(jj)==2) THEN
 m_org(NCOND-3)=m_org(NCOND-3)+dm_orgdt(jj)*MX(jj+8)/MX(NSPEC_P-3)*dt ! Oligomer SOA of compounds still contaning 1 -O-O-H groups
 m_org(jj) = m_org(jj)-dm_orgdt(jj)*dt
 ELSE IF (nrOOH(jj)<=1) THEN !(All other compounds reacting end up in the last oligomer SOA species class with no remaning -OOH groups)
 m_org(NCOND-2)=m_org(NCOND-2)+dm_orgdt(jj)*MX(jj+8)/MX(NSPEC_P-2)*dt ! Oligomer SOA of compounds contaning 0 -O-O-H groups
 m_org(jj) = m_org(jj)-dm_orgdt(jj)*dt
 END IF 
 END DO
 
 ! Consider fragmentation of HOM species:
 DO jj=635,742
 m_org(jj) = m_org(jj)-(dm_orgdt_frag(jj)+dm_orgdt_frag2(jj))*dt
 m_org(NCOND-10) = m_org(NCOND-10)+MX(jj+8)/MX(NSPEC_P-10)*0.1017*dm_orgdt_frag(jj)*dt ! C2 species
 m_org(NCOND-9) = m_org(NCOND-9)+MX(jj+8)/MX(NSPEC_P-9)*0.1966*dm_orgdt_frag(jj)*dt ! C4 species
 m_org(NCOND-8) = m_org(NCOND-8)+MX(jj+8)/MX(NSPEC_P-8)*0.3017*dm_orgdt_frag(jj)*dt ! C6 species
 m_org(NCOND-7) = m_org(NCOND-7)+MX(jj+8)/MX(NSPEC_P-7)*0.4000*dm_orgdt_frag(jj)*dt ! C8 species
 m_org(NCOND-6) = m_org(NCOND-6)+MX(jj+8)/MX(NSPEC_P-6)*dm_orgdt_frag2(jj)*dt ! C10 species from C20 fragmentation 
 END DO
 
  ! Consider fragmentation of C20-Peroxyhemiacetals species forming C10 fragments:
 DO jj=NCOND-5,NCOND-2
 m_org(jj) = m_org(jj)-dm_orgdt_frag2(jj)*dt
 m_org(NCOND-6) = m_org(NCOND-6)+MX(jj+8)/MX(NSPEC_P-6)*dm_orgdt_frag2(jj)*dt ! C10 species from C20 fragmentation 
 END DO


 corgp(:,j)=m_org*(SUM(corgpold(:,j)*MX(9:NSPEC_P))+1D-100) ! Convert it back to mole m^-3 air

 END DO
 END IF
    
corgpold=corgp;  

      DO j=1,nr_bins
        corgp_eq(:,j)=Corgold*SUM(corgpold(:,j))/(S_Kelvin(9:NSPEC_P,j)*psat_org/(Rg*T)) ! Approx equilibrium concentration (mol/m^3) of each compound in each size bin
      END DO

      ! Condensation of organic compounds
      !write(*,*) psat_org(NCOND)
      !write(*,*) corg(NCOND)
      !write(*,*) SUM(korg(NCOND,:))

      DO j=1,NCOND-1 ! Don't consider condensation of HOA
       IF (psat_org(j)<1D-2 .AND. corg(j)>1D4 .OR. j==(NCOND-1)) THEN
        xorg=corgpold(j,:)/(SUM(corgpold,DIM=1)+SUM(c_p_surf(7:8,:),DIM=1)/Na*1D6+1D-100) ! mole fraction of each organic compound in the organic + water particle phase
        Corgsat=xorg*psat_org(j)/(Rg*T)  ! Saturation concentration (mol/m^3)  
        Corgm=(Corgold(j)+dt*SUM(korg(j,:)*S_Kelvin(8+j,:)*Corgsat))/(1D0+dt*SUM(korg(j,:)))
        Corgtot=Corgold(j)+SUM(corgpold(j,:)) ! total vapor + particle mole concentration (mol m^-3 air) 
        corgp(j,:)=corgpold(j,:)+dt*korg(j,:)*(MIN(Corgm,Corgtot)-S_Kelvin(8+j,:)*Corgsat)
        WHERE (corgp(j,:)<0D0) corgp(j,:)=0D0
!        WHERE (corgp(j,:)>corgp_eq(j,:)) corgp(j,:)=corgp_eq(j,:) ! Don't allow particles to grow more than to the saturation concentration limit
        WHERE (corgp(j,:)>corgp_eq(j,:) .AND. Corgold(j)<psat_org(j)/(Rg*T)*S_Kelvin(8+j,:)) corgp(j,:)=corgp_eq(j,:) ! Don't allow particles to grow more than to the saturation concentration limit     
        Corgm=Corgtot-SUM(corgp(j,:)) ! mol m^-3
        IF (Corgm < 0D0) Corgm = 0D0
        corg(j)=Corgm*1D-6*Na  ! molecules cm^-3
       END IF
      END DO

      ! H2SO4:
      S_KelvinH2SO4=S_Kelvin(1,:)
      CH2SO4m=(CH2SO4old+dt*SUM(kH2SO4*S_KelvinH2SO4*CH2SO4sat))/(1D0+dt*SUM(kH2SO4)) ! condensation
      CH2SO4tot=CH2SO4old+SUM(cH2SO4pold) ! total vapor + particle mole concentration (mol m^-3 air) 
      cH2SO4p=cH2SO4pold+dt*kH2SO4*(MIN(CH2SO4m,CH2SO4tot)-S_KelvinH2SO4*CH2SO4sat)
      WHERE (cH2SO4p<0D0) cH2SO4p=0D0
      CH2SO4m=CH2SO4tot-SUM(cH2SO4p) ! mol m^-3
      !cH2SO4=CH2SO4m*1D-6*Na ! molecules cm^-3
    
    ! Consider SO2(aq) + H2O2(aq) -> S(VI) in cloud droplets:
    IF (LWC>=1D-6) THEN ! Assume that clouds are present in the grid cell if RH>=98 %
    CALL cloud_activation_processing(c_p,N_bins,cSO2,cH2O2,T,S_c,S_super,I_nSVI,LWC)
    
    ! Update the particle and gas-phase chemical composition after cloud processing and heterogeneous S(VI) formation
    where (S_c<=S_super) cH2SO4p=cH2SO4p+I_nSVI*N_bins*dt
   
    cSO2=MAXVAL((/cSO2-SUM(I_nSVI*Na*PACK(N_bins, S_c <= S_super)*1D-6)*dt,1D-10/))
    cH2O2=MAXVAL((/cH2O2-SUM(I_nSVI*Na*PACK(N_bins, S_c <= S_super)*1D-6)*dt,1D-10/))
    END IF
    

       S_KelvinHNO3=S_Kelvin(2,:);

      IF (thermodyn_index==1) THEN

        ! UPDATE HNO3(g):

        ! Solid salt NH4NO3 particles:
        WHERE (CHNO3s>0D0) term1=kHNO3*S_KelvinHNO3*CHNO3s*dt
        ! Liquid NH4NO3 particles:
        WHERE (CHNO3s<=0D0) term1=cHNO3pold*(1D0-EXP(-dt*S_KelvinHNO3*kHNO3/Kprim_HNO3))

        ! Solid salt NH4NO3 particles:
        WHERE (CHNO3s>0D0) term2=kHNO3*dt
        ! Liquid NH4NO3 particles:
        WHERE (CHNO3s<=0D0) term2=Kprim_HNO3/S_KelvinHNO3*(1D0-EXP(-dt*S_KelvinHNO3*kHNO3/Kprim_HNO3))

        CHNO3m=(CHNO3old+SUM(term1))/(1D0+SUM(term2))
     
        ! UPDATE NO3(p):
        term1=0D0
        WHERE (CHNO3s>0D0) term1=cHNO3pold
        CHNO3tot=CHNO3old+SUM(term1) ! total vapor + particle mole concentration of solid NH4NO3 particles (mol m^-3 air) 

        WHERE (CHNO3s>0D0) cHNO3p=cHNO3pold+dt*kHNO3*(MIN(CHNO3m,CHNO3tot)-S_KelvinHNO3*CHNO3s)
        WHERE (CHNO3s<=0D0) cHNO3p=Kprim_HNO3*CHNO3m/S_KelvinHNO3+(cHNO3pold-Kprim_HNO3*CHNO3m/ &
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
!         cHNO3p=cHNO3pold

        ! UPDATE HCl(g):
        term1=0D0; term2=0D0

        S_KelvinHCl=S_Kelvin(3,:);

        ! Solid salt NH4Cl particles:
        WHERE (CHCls>0D0) term1=kHCl*S_KelvinHCl*CHCls*dt
        ! Liquid NH4Cl particles:
        WHERE (CHCls<=0D0) term1=cHClpold*(1D0-EXP(-dt*S_KelvinHCl*kHCl/Kprim_HCl))

        ! Solid salt NH4Cl particles:
        WHERE (CHCls>0D0) term2=kHCl*dt
        ! Liquid NH4Cl particles:
        WHERE (CHCls<=0D0) term2=Kprim_HCl/S_KelvinHCl*(1D0-EXP(-dt*S_KelvinHCl*kHCl/Kprim_HCl))

 CHClm=(CHClold+SUM(term1))/(1D0+SUM(term2))
     
 ! UPDATE Cl(p):
 term1=0D0
 WHERE (CHCls>0D0) term1=cHClpold
 CHCltot=CHClold+SUM(term1) ! total vapor + particle mole concentration of solid NH4Cl particles (mol m^-3 air) 

 WHERE (CHCls>0D0) cHClp=cHClpold+dt*kHCl*(MIN(CHClm,CHCltot)-S_KelvinHCl*CHCls)
 WHERE (CHCls<=0D0) cHClp=Kprim_HCl*CHClm/S_KelvinHCl+(cHClpold-Kprim_HCl*CHClm/ &
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
!        cHClp=cHClpold

 ! Calculate the new gas phase NH3 concentration after equilibration with
 ! the particle phase according to the Newton-Raphson iteration in Jacobson, 2005 (eq. 17.126)
 cCO3p=mCO3*W*N_bins ! mol/m^3 CO32-
 cHCO3p=mHCO3*W*N_bins ! mol/m^3 HCO3-
 cOHp=mOH*W*N_bins ! mol/m^3 OH-

 charge=-fHSO4*cH2SO4p-2D0*fSO4*cH2SO4p-fCl*cHClp-fNO3*cHNO3p+cNap ! charge imbalance after condensation of HNO3, HCl and H2SO4 

 WHERE (charge>0D0) charge=0D0 ! correct for positive charge imbalance which is only possible if the acids evaporates

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

 c_p(Nlayers,2,:)=cHNO3p*1D-6*Na
 c_p(Nlayers,3,:)=cHClp*1D-6*Na
 c_p(Nlayers,4,:)=cNH3p*1D-6*Na
      END IF

      c_p(Nlayers,1,:)=cH2SO4p*1D-6*Na
      c_p(Nlayers,9:NSPEC_P-1,:)=corgp(1:NCOND-1,:)*1D-6*Na


      ! Output particle properties:
      c_p_tot(1:NSPEC_P,1:nr_bins)=sum(c_p,DIM=1)
      DO j = 1,nr_bins
        vp_dry(j)=SUM(c_p_tot(index_dry,j)/Na*MX(index_dry)/qX(index_dry)/(N_bins(j)*1D-6)) ! m^3     
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
              c_p_fixed(:,:,a-1) = c_p_fixed(:,:,a-1)+vp_fixed(a-1)/vp_dry(j)*r1*c_p(:,:,j)
            ELSE
              N_bins_fixed(a-1) = N_bins_fixed(a-1)+r1*N_bins(j)
              N_bins_fixed(a) = N_bins_fixed(a)+r2*N_bins(j)
              c_p_fixed(:,:,a-1) = c_p_fixed(:,:,a-1)+vp_fixed(a-1)/vp_dry(j)*r1*c_p(:,:,j)
              c_p_fixed(:,:,a) = c_p_fixed(:,:,a)+vp_fixed(a)/vp_dry(j)*r2*c_p(:,:,j)
            END IF
          ELSE ! evaporation
            IF (a == 0) THEN           
            ELSEIF (a == 1) THEN
              r1 = (vp_dry(j))/(vp_fixed(1)) ! Fraction of particles in size bin a
              N_bins_fixed(1) = N_bins_fixed(1)+r1*N_bins(j)
              c_p_fixed(:,:,1) = c_p_fixed(:,:,1)+vp_fixed(1)/vp_dry(j)*r1*c_p(:,:,j)
            ELSE
              r1 = (vp_dry(j)-vp_fixed(a-1))/(vp_fixed(a)-vp_fixed(a-1)) ! Fraction of particles in size bin i
              r2 = 1-r1 ! Fraction of particles in previous size bin (i-1)
              N_bins_fixed(a) = N_bins_fixed(a)+r1*N_bins(j)
              N_bins_fixed(a-1) = N_bins_fixed(a-1)+r2*N_bins(j)
              c_p_fixed(:,:,a) = c_p_fixed(:,:,a)+vp_fixed(a)/vp_dry(j)*r1*c_p(:,:,j)
              c_p_fixed(:,:,a-1) = c_p_fixed(:,:,a-1)+vp_fixed(a-1)/vp_dry(j)*r2*c_p(:,:,j)
            END IF
          END IF
        END DO
        c_p=c_p_fixed
        N_bins=N_bins_fixed
!        DO j = 1,nr_bins
!          IF (N_bins(j) <= 0.) WRITE(*,'(A14,I5)'), 'N_bins <= 0 at', tr
!          IF (N_bins(j) <= 0.) WRITE(*,*) N_bins 
!        END DO
      END IF


      c_p_tot(1:NSPEC_P,1:nr_bins)=sum(c_p,DIM=1)

      DO j = 1,nr_bins

        IF (N_bins(j)<=0.1D0) THEN ! This can occur with the full-stationary method if some particles grow to more than the next particle size bin within one condensation time step
          N_bins(j)=1D0
!          c_p(:,:,j) = c_p_old(:,:,j)/N_bins_old(j)
!          c_p_tot(:,j)=SUM(c_p(:,:,j),DIM=1);

          c_p(:,:,j)=c_p_backg(:,:,j);
          c_p_tot(:,j)=SUM(c_p(:,:,j),DIM=1);
          
!          c_p(:,:,j)=c_p_old(:,:,j)*N_bins(j)/N_bins_old(j)
!          c_p_tot(:,j)=sum(c_p(:,:,j),DIM=1)
!          c_p(:,:,j)=c_p_old(:,:,j)*N_bins_old(j)/SUM(N_bins_old)
!          c_p_tot(:,j)=sum(c_p(:,:,j),DIM=1)
        END IF
       
        vp_dry(j)=SUM(c_p_tot(index_dry,j)/Na*MX(index_dry)/qX(index_dry)/(N_bins(j)*1D-6)) ! m^3     
        vp_wet(j)=SUM(c_p_tot(:,j)/Na*MX/qX/(N_bins(j)*1D-6)) ! m^3
        V_bins(j)=N_bins(j)*vp_wet(j)
        dens(j)=SUM(c_p_tot(:,j)*MX*1D6)/Na/V_bins(j) ! Total particle density
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
      ! Variables used to calculate friction velocity
      !REAL(dp) :: flux_horiz,u
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
      rac(4,:) = (/ 0.,0.,0.,0.,0. /) ! barren
      rac(5,:) = (/ 200.,140.,120.,50.,120. /) ! shrubs
      rac(6,:) = (/ 0.,0.,0.,0.,0. /) ! water
      rac(7,:) = (/ 100.,100.,100.,100.,100. /) ! urban
      rac(8,:) = (/ 0.,0.,0.,0.,0. /) ! snow and ice
      ! resistance for uptake by soil, leaf litter, and so on at the ground SO2
      rgsS(1,:) = (/ 500.,500.,500.,100.,500. /) ! evergreen
      rgsS(2,:) = (/ 500.,500.,500.,100.,500. /) ! deciduous
      rgsS(3,:) = (/ 0.,0.,0.,100.,0. /) ! wetland
      rgsS(4,:) = (/ 1000.,1000.,1000.,1000.,1000. /) ! barren
      rgsS(5,:) = (/ 400.,400.,400.,50.,400. /) ! shrubs
      rgsS(6,:) = (/ 0.,0.,0.,0.,0. /) ! water
      rgsS(7,:) = (/ 400.,400.,400.,100.,500. /) ! urban
      rgsS(8,:) = (/ 0.,0.,0.,0.,0. /) ! ! snow and ice
      ! restistance for uptake by soil, leaf litter, and so on at the ground, O3
      rgsO(1,:) = (/ 200.,200.,200.,3500.,200. /) ! evergreen
      rgsO(2,:) = (/ 200.,200.,200.,3500.,200. /) ! deciduous
      rgsO(3,:) = (/ 1000.,800.,1000.,3500.,1000. /) ! wetland
      rgsO(4,:) = (/ 400.,400.,400.,400.,400. /) ! barren
      rgsO(5,:) = (/ 200.,200.,200.,3500.,200. /) ! shrubs
      rgsO(6,:) = (/ 2000.,2000.,2000.,2000.,2000. /) ! water
      rgsO(7,:) = (/ 300.,300.,300.,600.,300. /) ! urban
      rgsO(8,:) = (/ 2000.,2000.,2000.,2000.,2000. /) ! snow and ice
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
      !flux_horiz = dens_air*((-UMOF/dens_air)**2.+(-VMOF/dens_air)**2.)**0.5 ! Vertical turbulent flux of horizontal momentum
      !u = (flux_horiz/dens_air)**0.5;            ! Friction velocity (m/s)
        
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

!     IF (landuse_index == 6) THEN ! surface resistance for water (19.5.1 Seinfeld and Pandis)
!       kg = 0.0013*windspeed ! Hicks and Liss 1976, gas phase mass transfer coeff
!       H_T0 = (/ 1.23D0,1.1D-2,1D-2,1.9D-3,2.1D5,1D5,2.5D0,310D0,49D0,2.1D5,1.1D0,62D0 /) ! Henry's law coeff. at T=298K from table 7.2 Seinfeld and Pandis (mol L⁻¹ atm⁻¹)
!       H_exp = (/ -3020.,2300.,-2500.,-1480.,-8650.,-6800.,-6500.,-5600.,-4800.,-8650.,-2020.,-3400. /) ! Henry's law exponent used to calculate the temp dep henry's law coeff (table 19.4 Seinfeld and Pandis)
!       H = H_T0*exp(H_exp*(1./298.-1./temp)) ! temp dependent henry's law coeff, table 19.4 (mol L⁻¹ atm⁻¹)
!       H_dimless = H*Rg*0.009869*temp ! Dimensionless Henry's law coeff, see App A in Jacobson for unit conversion
!       Sc_CO2 = 2073.1-125.62*19.85+3.6276*19.85**2.-0.043219*19.85**3. ! The Schmidth nr of CO2 in seawater at 293 K (19.85 degrees Celcius)
!       Diff_water = (/ 1D-9,1D-9,1D-9,1D-9,1D-9,1D-9,1D-9,1D-9,1D-9,1D-9,1D-9,1D-9 /)! gas diffusivity in water (typical values of any gas in water, fundamentals of environmental engineering)
!       Sc_water = 1D-6/Diff_water ! Sc nr
!       Sc_ratio = Sc_CO2/Sc_water ! Sc nr ratio at the temp of interet
!       ! Estimations of liquid-phase mass transfer coeff (Liss and Merlivat, 1986)
!       IF (windspeed <= 3.6) THEN
!         kl = (0.17*windspeed*Sc_ratio**(2./3.))/3600.*1D-2 ! (m/s)
!       ELSE IF (windspeed > 3.6 .and. windspeed <= 13.) THEN
!         kl = (0.612*Sc_ratio**(2./3.)+(2.85*windspeed-10.26)*Sc_ratio**0.5)/3600.*1D-2
!       ELSE
!         kl = (0.612*Sc_ratio**(2./3.)+(5.9*windspeed-49.9)*Sc_ratio**0.5)/3600.*1D-2
!       END IF
!       rc = 1./kg+1./(kl*H_dimless) ! (s/m)

 ! Parameterization from Zhang et al 2003 where they set rc_ground(O3)=2000s/m and rc_ground(SO2)=20s/m and scale the other gases with coefficients from Zhang et al 2002
 !    alpha = 1D0
 !    bet = 1D0
 !    alpha(1:12) = (/ 1.,0.,0.,0.,10.,1.,0.8,0.8,2.,1.,10.,1./) ! table 1 in zhang et al 2002
 !    bet(1:12) = (/ 0.,1.,0.8,0.8,10.,1.,0.2,0.2,2.,1.,10.,0./) ! table 1 in Zhang et al 2002
 !    rc = 1./(alpha/20.+bet/2000.) ! s/m       
 !    ELSE
     ! HENRYWIN calculations gives H_eff>1D12 M/atm for HOM, all  H_eff>1D10  gives 0 resistance.
       H_eff(1:12) = (/ 1D5,1D-2,1D-2,2D-3,1D14,1D5,6D3,220D0,1D5,1D14,1D14,2D4/) ! effective henry's lay const, table 19.4, H2SO4 and hcl assumed same as nitric acid (both strong acids)
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
       IF (rain <= 0.1) THEN 
        rlui(1) = (1./5000.+1./(3*rlui(1)))**(-1) ! Upper canopy resistance for so2 when surface is covered by rainwater
        rlui(2) = (1./1000.+1./(3*rlui(1)))**(-1) ! Upper canopy resistance for o3 when surface is covered by rainwater
        rlui(3:10) = (1/(3*rlui(3:10))+1D-7*H_eff(3:10)+f0(3:10)/rlui(2))**(-1) ! Upper canopy resistance for rest of species when surface is covered by rainwater 
       END IF     
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
 !    END IF   

      ! dry deposition velocity                                         
      vd_gas = 1./(ra+rb+rc) ! [m/s]

      
      WHERE (vd_gas >= 0.1) vd_gas = 0.1

    END SUBROUTINE dry_dep_gases

END MODULE dynamics
