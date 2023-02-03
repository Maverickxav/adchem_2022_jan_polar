PROGRAM adchem1D_new
! Updated Photolysis rates and new corrected actinic flux input data with the sum of diffuse upward and downward radiation included (Photolysis rates evaluatd against the TUV model with good agreement)
! Updated dry deposition losses during precipitation events.
! Updated the in-cloud wet scavening based on EMEP model parameterization and scavenging of only particles activated as cloud droplets.
! Updated Kz above the PBLH to 0.1 m^2/s based on Flexpart 10.4 description
! Don't update IO3- concentrations after multiphase chemistry module (KPP_Proceed)
  
    ! Chemistry modules:
    USE second_Main
    USE second_Precision,  ONLY : dp    ! KPP Numerical type
    USE second_Parameters ! before: USE adchem_Parameters, ONLY : NSPEC ! Number of chemical species
    USE second_Monitor
    USE second_Global, ONLY : NPHOT,NKVALUES ! Number of photolysis parameters
    USE reaction_rates_20220118, ONLY : getKVALUES, getJVALUES

    USE thermodynamicsDMS 
    USE dynamicsDMS_atm ! Aerosol dynamics
    USE constants ! Set case specifics (module also includes physical constants)
    USE diffusivity ! calculate k-diffusivity coefficient
   ! USE multilayer
    USE output

    USE Megan_version_2
    USE LAI_Hyy_month

    !! For clusterin
    Use ClusterIn_plugin
    Use acdc_datatypes
    include 'cluster_chem_use.inc'

    IMPLICIT NONE

    REAL(dp) :: start_time,finish_time
    
    ! Variables for atmospheric diffusion calculations:
    REAL(dp), DIMENSION(Nz+1) :: Kz  ! Eddy diffusion coefficients in the vertical direction (0 m - 2000 m)
    REAL(dp), DIMENSION(Nz,Nz) :: c_diff,c_diff_O3 ! Representation of diffusion pattern c_diff(1,:) represent diffusion from surface layer
    REAL(dp)    ::  k_uper_BC=-1D-3 ! slope of conc. gradient as uper BC for all compounds except O3 (dc/dz=k_uper_BC*c)
    REAL(dp)    ::  k_uper_O3=1D-3 ! slope of conc. gradient as uper BC for all compounds except O3 (dc/dz=k_uper_BC*c)
    REAL(dp)    ::  k_uper_BC_canopy ! slope of conc. gradient at the canopy border (25 m) 
    REAL(dp), DIMENSION(Nz) :: conc_diff
    REAL(dp), DIMENSION(traj_len_min) :: k_uper_BC_O3
    REAL(dp), DIMENSION(Nz) :: dz
    REAL(dp), DIMENSION(Nz-1) :: dz2

    
    ! Variables along trajectory 
    REAL(dp), DIMENSION(Nz+1,traj_len_min) :: Ts, Ps, RHs, u_vel, v_vel, cloudLWC 
	REAL(dp), DIMENSION(Nz+1)              :: RH,RH_old=95D0
	REAL(dp), DIMENSION(Nz,traj_len_min) :: O2, N2, M, H2O
	REAL(dp)                             :: M_STP 
    REAL(dp), DIMENSION(traj_len_min) :: fric_veloc,PBLH,SHTFs,rain,windspeed,DSWF,T2m,P2m,RH2m,SST,latitude
	REAL(dp) :: LWC,Win,hs

    
	REAL(dp), DIMENSION(12,traj_len_min) :: EmissionCO, EmissionNO2, EmissionSO2, EmissionNH3, EmissionNMVOC
	
	REAL(dp), DIMENSION(traj_len_min) :: EmissionAPIN, EmissionBPIN, EmissionLIM, EmissionBCAR, &
    EmissionOTH_MT, EmissionOTH_SQT, Emission_all_bvoc, EmissionDMS,EmissionCH2Br2,EmissionCHBr3,&
	EmissionCH3I, EmissionISO, EmissionNH3_birds,EmissionSO2_ship,EmissionC3H7I,&
    EmissionC2H5I,EmissionCH2I2,EmissionCH2ICL,EmissionCH2IBr,EmissionDIBRET,&
    EmissionCH3BR,EmissionTCE,EmissionTRICLETH,EmissionCHCL3,EmissionCH2CL2,EmissionCH3CL
	REAL(dp) :: EmissionI2,EmissionHOI,Iion_aq,O3_ppb
	
    REAL(dp), DIMENSION(traj_len_min) :: EmissionC2H6, EmissionC3H8, EmissionNC4H10, EmissionC2H4, EmissionC3H6, EmissionC5H8, &
	EmissionOXYL, EmissionCH3OH, EmissionC2H5OH, EmissionHCHO, EmissionCH3CHO,EmissionC2H5CHO,EmissionMEK, EmissionGLYOX, &
	EmissionMGLYOX, EmissionC9_C12_alkanes,EmissionNOxsoil, EmissionNOx_road,EmissionNOx_nonroad_mach
	
	REAL(dp), DIMENSION(traj_len_min) :: EmissionC2H6_Pplants_Ind, EmissionNC4H10_Pplants_Ind, EmissionC2H4_Pplants_Ind, &
	EmissionC3H6_Pplants_Ind, EmissionC5H8_Pplants_Ind, EmissionOXYL_Pplants_Ind, EmissionCH3OH_Pplants_Ind, &
	EmissionC2H5OH_Pplants_Ind, EmissionHCHO_Pplants_Ind, EmissionCH3CHO_Pplants_Ind, EmissionMEK_Pplants_Ind, EmissionGLYOX_Pplants_Ind, &
	EmissionMGLYOX_Pplants_Ind
	
	REAL(dp), DIMENSION(traj_len_min) :: EmissionC2H6_Ship, EmissionNC4H10_Ship, EmissionC2H4_Ship, &
	EmissionC3H6_Ship, EmissionC5H8_Ship, EmissionOXYL_Ship, EmissionCH3OH_Ship, &
	EmissionC2H5OH_Ship, EmissionHCHO_Ship, EmissionCH3CHO_Ship, EmissionMEK_Ship, EmissionGLYOX_Ship, &
	EmissionMGLYOX_Ship
   
	REAL(dp), DIMENSION(14) ::	VOC_dist_A_PublicPower,VOC_dist_B_Industry,VOC_dist_C_OtherStationaryComb,&
	VOC_dist_D_Fugitiv,VOC_dist_E_Solvents,VOC_dist_F_RoadTransport,VOC_dist_G_Shipping,VOC_dist_H_Aviation,&
	VOC_dist_I_Offroad,VOC_dist_J_Waste,VOC_dist_K_AgriLivestock,VOC_dist_L_AgriOther,VOC_dist_M_Other
	
	REAL(dp), DIMENSION(18,Nz)   :: init_conc
	REAL(dp), DIMENSION(Nz)      :: PM_BC_init,PM_OC_init,PM_SO4_init,PM_NH4_init,PM_Na_init,PM_NaCl1_init,PM_NaCl2_init,PM_NaCl3_init,&
	PN_BC,PN_OC,PN_NaCl1,PN_NaCl2,PN_NaCl3,PN_SO4,PN_NH4,PN_Na
	REAL(dp)                     :: spm_BC,spm_OC,spm_NaCl1,spm_NaCl2,spm_NaCl3,spm_SO4,spm_NH4,spm_Na
	REAL(dp), DIMENSION(nr_bins) :: Nbins_BC,Vbins_BC,Nbins_OC,Vbins_OC,Nbins_SO4,Vbins_SO4,Nbins_NH4,Nbins_Na,Vbins_NH4,Vbins_Na,Nbins_NaCl1,Vbins_NaCl1,Nbins_NaCl2,Vbins_NaCl2,Nbins_NaCl3,Vbins_NaCl3
	REAL(dp), DIMENSION(Nz,nr_bins)   :: NbinsBC,NbinsOC,NbinsSO4,NbinsNaCl,NbinsNH4,NbinsNa


    REAL(dp) :: scale_bio
    REAL, DIMENSION(27,traj_len_min) :: input_MEGAN
   
    REAL(dp) :: PNship, PNroad,PNnonroad_mach,PNwood, PNbiomass
    INTEGER  :: index_st = 7*24*60 ! trajectory (tr) nr where air mass is closest to station
    INTEGER, DIMENSION(210) :: list_landcat
    INTEGER :: land_index
    REAL,     DIMENSION(22,365) :: EF_SMEARII ! BVOC standard emission potential from SMEAR II
    INTEGER, DIMENSION(traj_len_min) :: landuse_index, landuse_index_rev ! From GLC2000 (run convert1.m for wanted geographical area, then to get the land-use category along the trajectory run EmissionEMEOtraj  ! 1 = evergreen
                 ! 2 = deciduous
                 ! 3 = grass
                 ! 4 = desert
                 ! 5 = shrubs
                 ! 6 = sea
                 ! 7 = urban area

    ! Gas-phase specific variables (all gases)
    REAL(dp), DIMENSION(Nz,NSPEC) :: conc ! Vector with concentrations of all gases 

    REAL(dp), DIMENSION(6,NSPEC) :: E_gases
    REAL(dp), DIMENSION(NSPEC)   :: conc1
    REAL(dp), DIMENSION(NCOND+12):: vd_gas
    INTEGER, DIMENSION(11)       :: index_dep_gases                 
    INTEGER                      :: index_MCM(NSPEC)
    REAL(dp),DIMENSION(Nz+1)     :: NOx_initial,SO2_initial,O3_initial,CO_initial !vertical profile of some gases from ECMWF (50:100:2050 m)
	! Gas-phase specific variables (only the condensable species)
    REAL(dp), DIMENSION(52,NCOND) :: subgroupsUNIFAC
    REAL(dp), DIMENSION(60,NCOND) ::UNIFAC_Nannoolal
    REAL(dp), DIMENSION(NCOND)   :: A_Nannoolal, B_Nannoolal,Morg,NrC,NrO,NrH,NrN,NrS,VolX,Dorg
    REAL(dp), DIMENSION(NCOND) :: psat_org,y_org_water,Henrys_coeff ! Variables used to calculate the saturation vapor pressure (Nannoolal method) (pure vapor pressure coeff from http://umansysprop.seaes.manchester.ac.uk/)
    REAL(dp), DIMENSION(NCOND) :: indexROOH=1D0, indexALDEHY=1D0 
    INTEGER, DIMENSION(NCOND)  :: index_cond 
    
	! Also include inorganic species and water 
    REAL(dp), DIMENSION(NSPEC_P) :: MX, qX, VX, molec_radius, aX
    CHARACTER(LEN=20), DIMENSION(NCOND) :: mcm_name
    REAL(dp), DIMENSION(7,NCOND) :: mcm_prop
    REAL(dp), DIMENSION(Nz) :: NO_NO2_ratio
    
    REAL(dp), DIMENSION(510,traj_len_h*(Nz+1)) ::  actinic_flux_1h
    REAL(dp), DIMENSION(510,Nz+1) ::  actinic_flux_t1,actinic_flux_t2
	REAL(dp) :: x1, x2
	
    REAL(dp), DIMENSION(510) ::  actinic_flux1       
    REAL(dp) :: lambda1(451), lambda(511), lambda2(20), lambda3(29), lambda4(11), dlambda(510)

    ! Variables used in gas-phase chemistry
    REAL(dp) :: cond_sinkH2SO4,cond_sinkHNO3
     
    REAL(dp)      :: t_start, t_end, t_start_chem 
    REAL(dp)      :: JVALUES(NPHOT),KVALUES(NKVALUES) ! photolysis/reaction rate coefficients
    REAL(dp)      :: dt,dt_cond,time, time_day
    REAL(dp), DIMENSION(Nz)  :: c222Rn,cHOA,cHOM, q_222Rn, q_ion=2D0, n_ion=1D0;
    REAL(dp)      :: E_222Rn, k_222Rn, q_GCR, ipr, f_DMA_NH3=1D-2

    ! Particle properties in each size-bin 
    REAL(dp), DIMENSION(nr_bins+1) :: d,vp
    REAL(dp), DIMENSION(nr_bins) :: d_g, v_p, v_dep, vs=0D0, dlogDp, delta_L2, vp_c, vp_b, &
    V_bins1, N_bins1, d_p1, dp_dry1, pH1, Kprim_HNO31, Kprim_HCl1, Kprim_CH3SO3H1,Kprim_HIO31,Hprim_NH31, Kprim_NH31, fHSO41, &
    fSO41, fNO31, fCl1, fCH3SO31, fHIO31, mHCO31, mCO31, mOH1, mCOO1, W1, dens_p1, &
    V_bins_soot, PN_emission_ship, PV_emission_ship, N_bins_soot, v_wet1, &
    PN_emission_S7, PV_emission_S7,PN_emission_S9, PV_emission_S9, f_POA_marine, viscosity,surf_tens,S_Kelvin, &
    aw,E_sea_salt,E_snow_salt,exp_vwet,V_bins_dry,&
	PN_emission_agr, PV_emission_agr,&
    PN_emission_waste, PV_emission_waste,&
    PN_emission_traffic, PV_emission_traffic,&
    PN_emission_agr_burn, PV_emission_agr_burn,&
    PN_emission_wood, PV_emission_wood,&
    PN_emission_coal, PV_emission_coal,&
    PN_emission_other_burn, PV_emission_other_burn,&
    PN_emission_industry, PV_emission_industry,&
    PN_emission_power_plant, PV_emission_power_plant,&
    PN_emission_flaring, PV_emission_flaring,&
	dimer_C1, dimer_O1, dimer_N1, dimer_H1,PN_biomass_burn, PV_biomass_burn
	
	REAL(dp), DIMENSION(6,nr_bins) :: PN_marine,V_bins_sea_spray
	
	REAL(dp) :: PV1,PV2,PV3,PVtot,xPV1,xPV2,xPV3,xPM01,xPM1,xPM25
    REAL(dp), DIMENSION(traj_len_min) :: PM01,PM1,PM25

    REAL(dp), DIMENSION(10,traj_len_min) :: PN_emissions_GAINS
    REAL(dp), DIMENSION(15,10) :: PN_emissions_dist_param_GAINS
    REAL(dp), DIMENSION(traj_len_min) :: DM_emissions_biomass_burning=0D0
	REAL(dp), DIMENSION(nr_bins,traj_len_min) :: E_sea_salt_,E_snow_salt_
    
    REAL(dp) :: NMVOC_corr_factor_wood ! Scale the domestic wood emissions with a factor which depends on the temperature 
    REAL(dp) :: EF_T ! Temperature road traffic primary particle emission factor
    
 
    REAL(dp), DIMENSION(Nz,nr_bins) :: N_bins,V_bins,vp_dry, vp_wet, d_p, dp_dry, dens_p,&
    pH,Kprim_HNO3,Kprim_HCl,Kprim_CH3SO3H,Kprim_HIO3,Hprim_NH3,Kprim_NH3,fHSO4=1D0,fSO4=0D0,fNO3,fCl,fCH3SO3,fHIO3,mHCO3,mCO3,mOH,mCOO,&
    W, v_wet,v_wet_in_cloud=0D0, N_bins_old, loss_fraction ,V_bins_old, dimer_C, dimer_O, dimer_N, dimer_H
	
	REAL(dp), DIMENSION(Nz) :: cDMSaq=0D0,cDMSOaq=0D0,cMSIAaq=0D0,cHPMTFaq=0D0,cSO2aq=0D0,cH2O2aq=0D0,cOHaq=0D0,cO3aq=0D0	 
      
    ! Particle compound specific properties in each particle layer and size bin
    REAL(dp), DIMENSION(Nz,NSPEC_P,nr_bins) :: c_p,c_p_old  ! (molec cm^3) Matrix with the concentration of all species in each particle layer for each particle size 
    REAL(dp), DIMENSION(NSPEC_P,nr_bins) :: c_p_backg
	REAL(dp), DIMENSION(NSPEC_P)         :: c_p_nucl, c_p_nuclDMA, c_p_nuclIO
    REAL(dp), DIMENSION(NMLAYER) :: Vi, DiffX, dVidt
    REAL(dp), DIMENSION(NSPEC_P,nr_bins) :: c_p_L,c_p_L23
    REAL(dp), DIMENSION(NSPEC_P,nr_bins) :: c_p1
    REAL(dp), DIMENSION(Nz,NSPEC_P+1,nr_bins) :: y
    REAL(dp), DIMENSION(NSPEC_P+1,nr_bins) :: y1
    REAL(dp), DIMENSION(Nz) :: Jnucl_N,Jnucl_D, Jnucl_I,Jnucl_org, Jnucl_SAorg ! New particle formation rate # cm^-3 s^-1
    REAL(dp), DIMENSION(NCOND) :: corg
    REAL(dp), DIMENSION(NCOND,nr_bins) :: yorg
 
    ! kinetic 3-layer model parameters:
    REAL(dp), DIMENSION(nr_bins) :: vp_wet1
    REAL(dp) :: dp_b, dx_s, vp_b_new, vp_b1old, vp_s, vp_s_new, vf1, vf2,dp_c,dx_b
 
    ! Particle number size distribution parameters:
    INTEGER, PARAMETER  ::  modes=5
    INTEGER, PARAMETER  ::  modes3=3
    INTEGER, PARAMETER  ::  modes2=2
    !INTEGER, DIMENSION(Nlayers) :: indexx
    REAL(dp), DIMENSION(modes) :: s, N_modes, dm 
    REAL(dp), DIMENSION(modes3) :: s3, N_modes3, dm3
    REAL(dp), DIMENSION(modes2) :: s2, N_modes2, dm2
    !REAL                    :: N_tot_start
     
    INTEGER  :: i, ii, j, jj, iii, tr,tr_1h,k,index_core,index_core2,index_min,index_max,JD,index_gas,a
	REAL(dp) :: r1,r2

    ! Variables needed to calculate a-pin emission, to be removed later
    REAL(dp) :: E_pot, Bio_dens, T_s, beta, E_apinene, env_corr_tot

    ! MEGAN parameters
    ! All variables used in the emissions starting with the three characters 'EM_' to make life easier for implementation

    INTEGER                       :: Mon = 5                           ! Month
    INTEGER, PARAMETER            :: kz_can = 6                        ! Number of vertical layers inside the model (not only canopy)
    INTEGER, PARAMETER            :: Can_Lay = 5                       ! Number of vertical layers inside the canopy
    INTEGER, PARAMETER            :: EM_Year = 2017                    ! Year
    INTEGER                       :: EM_Cantype                        ! Canopy type index
    
    INTEGER                       :: EM_run                            ! Run parameter in emission module
    INTEGER                       :: EM_TIME_M2                        ! Time for new Megan in hour minutes second format
    INTEGER                       :: EM_DATE                           ! Date for new Megan in format YYYDDD scalar
    INTEGER                       :: EM_Julian                         ! Julian day of the year
    INTEGER                       :: EM_kz                             ! Amount of layers in the model run (not only canopy)
    INTEGER                       :: EM_Can_Lay                        ! Amount of layers inside the canopy

    REAL                          :: EM_DT                             ! Time step for the emission code
    REAL, PARAMETER               :: Rpi180 = 57.29578                 ! no comments

    REAL                          :: EM_TIME_M2_R                      ! Time in seconds
    REAL                          :: EM_Lat, EM_Long                   ! Latitude and longitude
    REAL, DIMENSION(traj_len_min) :: EM_PAR                            ! Photosynthetically active radiation
    REAL                          :: EM_PAR_day                        ! Daily average photosynthetically active radiation
    REAL                          :: EM_LAI                            ! Leave area index
    REAL                          :: EM_LAI_past                       ! LAI from the last month
    REAL                          :: EM_PRES                           ! Pressure in the canopy (Pa)
    REAL                          :: EM_SRAD                           ! Incoming short wave solar radiation in (W/m≤)
    REAL                          :: EM_SRAD_day                       ! Daily average short wave radiation (W/m2)
    REAL                          :: EM_TempK_day                      ! Daily average of temperature (K)
    REAL                          :: EM_SMOIST                         ! Soil moisture (%)
    REAL                          :: EM_LAI_year                       ! LAI factor for different years
    REAL(dp)                      :: hourangle
     
    REAL(dp), DIMENSION(kz_can)       :: EM_VPD_S, EM_VPD_C                ! Water vapour pressure deficit for sun and shade leaf temperatures

    REAL, DIMENSION(kz_can)           :: EM_z,EM_dz                        ! Array of height for each layer and layer width (m)
    REAL, DIMENSION(kz_can)           :: EM_LAD                            ! Array of vertical distribution for the LAI
    REAL, DIMENSION(kz_can)           :: EM_TempK                          ! Array of temperature (K)
    REAL, DIMENSION(kz_can)           :: EM_TempC                          ! Array of temperature (C)
    REAL, DIMENSION(kz_can)           :: EM_WIND                           ! Array of horizontal wind (m/s)
    REAL, DIMENSION(kz_can)           :: EM_RH                             ! Array of relative humidity (%)
    REAL, DIMENSION(kz_can)           :: EM_WVM                            ! Array of water vapour mixing ration (g/g)
    REAL, DIMENSION(kz_can)           :: EM_EW                             ! Array of water vapour pressure (Pa)
    REAL, DIMENSION(kz_can)           :: EM_ES, EM_ES_S, EM_ES_C           ! Array of saturation water vapour pressure: ambient, sun and shade leaf temperatures
    REAL, DIMENSION(12)               :: EM_LAI_month                      ! Array of data for different LAI per month
    REAL, DIMENSION(15)               :: EM_LAI_year_Hyy                   ! Array of data for different LAI per year (1996-2010)

    ! Output variables:
    REAL, DIMENSION(kz_can)           :: EM_Sunfrac                        ! Array of the fraction of sun leaves. i = 1 is the top canopy layer, 2 is the next layer, etc.
    REAL, DIMENSION(kz_can)           :: EM_Sun_Par                        ! Array of sun fraction - 1 above the canopy and decreases inside the canopy
    REAL, DIMENSION(kz_can)           :: EM_SunleafTK, EM_SunleafTC        ! Array of temparture for sun leaf in (K) and (C)
    REAL, DIMENSION(kz_can)           :: EM_ShadeleafTK, EM_ShadeleafTC    ! Array of temparture for shade leaf in (K) and (C)

    REAL, DIMENSION(kz_can,22)        :: EM_EMI                            ! Matrix of emissions for each layer and each compound
    REAL, DIMENSION(kz_can,22)        :: EMI                               ! Vector with BVOC emissions in normalised concentrations over the ADCHEM surface layer width (dz)
    REAL, DIMENSION(22)               :: EM_EFin                           ! Standard emission factor for the 22 different compounds in MEGAN
    REAL, DIMENSION(4,kz_can)         :: EM_Ea1tL_M2                       ! Array of emission activity of temperature per layer
    REAL, DIMENSION(4,kz_can)         :: EM_Ea1pL_M2                       ! Array of emission activity of light per layer
    REAL, DIMENSION(4,kz_can)         :: EM_Ea1NL_M2                       ! Array of companied emission activity
    REAL, DIMENSION(kz_can,22)        :: EM_ER_BT                          ! Array of output emission buffer
    REAL, DIMENSION(kz_can,22)        :: EM_ER_NT                          ! Array of output emission buffer
    REAL, DIMENSION(kz_can,22)        :: EM_ER_SB                          ! Array of output emission buffer
    REAL, DIMENSION(kz_can,22)        :: EM_ER_Grass                       ! Array of output emission buffer
    REAL, DIMENSION(20, kz_can)       :: EM_ER                             ! Array of output emission buffer
    REAL, DIMENSION(20,3)             :: EM_GAM_OTHER                      ! Array of other gamma factors
    REAL, DIMENSION(20,kz_can)        :: EM_GAM_TMP                        ! Array of temperature response factor
    
    ! Saturation vapour pressure of water in Pa
    REAL, PARAMETER               :: a0 = 6.107799961,               & ! Parameters to calculate the saturation vapour pressure for water
                                     a1 = 4.436518524E-1,            &
                                     a2 = 1.428945805E-2,            &
                                     a3 = 2.650648471E-4,            &
                                     a4 = 3.031240396E-6,            &
                                     a5 = 2.034080948E-8,            &
                                     a6 = 6.136820929E-11

    CHARACTER*16                  :: EM_VAR(22)                        ! Name of the 22 VOC's
    
	
	 ! ACDC parameters:
    REAL(dp)                      :: c_acid,c_base,c_org,diameter_acdc,diameter_acdcDMA,diameter_acdcHIO3,cDMA
    REAL(dp), DIMENSION(Nz)       :: CS_H2SO4=1D-3 ! s^-1
    REAL(dp), DIMENSION(Nz)       :: CS_air=1D-3 ! s^-1
    REAL(dp), DIMENSION(3)        :: Nuc_by_charge !  formation rate vector (for neu, neg and pos) 
    ! REAL(dp), DIMENSION(Nz,63)    :: c_clusters=0D0
	! REAL(dp), DIMENSION(Nz,31)    :: c_clustersDMA=0D0
    REAL(dp), allocatable     :: c_clusters_n(:,:), c_clusters_d(:,:),c_clusters_i(:,:),c_clusters1(:), c_clusters2(:),c_clusters3(:)
	! REAL(dp), DIMENSION(Nz,31)    :: c_clustersDMA=0D0
	! REAL(dp), DIMENSION(63)       :: c_clusters1
    ! REAL(dp), DIMENSION(31)       :: c_clusters2	
    LOGICAL                       :: solve_ss
 
    
    REAL(dp)                      :: correction_Primary_p = 1D0 ! Correction factor for anthropogenic primary particle emissions
    REAL(dp)                      :: correction_gas_emissions = 1D0 ! Correction factor for anthropogenic gas emissions


	REAL(dp)                      :: X_NH3,F_NH3ocean,kg_NH3,cNHx_ocean,Sc_NH3,C_D, &
	                                 H_NH3_sea,S_sea,pKa_NH3_sea,pKa_DMA_sea,pH_ocean,H_prim_NH3_sea,&
									 cMeNO3_ocean,cEtNO3_ocean,H_MeNO3,H_EtNO3,X_MeNO3,X_EtNO3,&
									 Sc_MeNO3,Sc_EtNO3,kg_MeNO3,kg_EtNO3,F_MeNO3ocean,F_EtNO3ocean,F_NOxsnow,F_NOxocean,&
									 cDMA_ocean, H_DMA_sea, X_DMA, Sc_DMA, kg_DMA, F_DMAocean, H_prim_DMA_sea

	
	REAL(dp), DIMENSION(nr_bins+1) :: vp_fixed
	REAL(dp), DIMENSION(nr_bins) :: N_bins_fixed
    REAL(dp), DIMENSION(NSPEC_P,nr_bins) :: c_p_fixed ! (molec cm^-3) 
	REAL(dp), DIMENSION(nr_bins)  :: S_c,c_H,ns,vp0,dp0,A_Kohler,B_Kohler,alpha1,alpha2,fCCN,fCS,b_ko,d_ko,lambda_ko1,lambda_ko2
    REAL(dp)                      :: S_super,CCN,dp_drops,m_drops,dp_max
    REAL(dp)            :: cw,Hion=1D-3,OHion
	REAL(dp), DIMENSION(52) :: alpha_aq,D_aq,dX_aq,VX_aq,c_speed_aq,CS_AQ,MX_aq
    REAL(dp), DIMENSION(nr_bins) :: Cunningh,Diff_p,m_p,speed_p,gasmeanfp_aq,Kn_aq,f_cor_aq,Deff_aq
    REAL(dp) :: dyn_visc,l_gas,fHCl
	REAL(dp), DIMENSION(NSPEC) :: conc_old,dconc
	
	REAL(dp) :: pH_bulk1=1D-4
	REAL(dp), DIMENSION(Nz) :: pH_bulk=1D-4
	REAL(dp), DIMENSION(8) :: yi_bulk1=1D0
	REAL(dp), DIMENSION(Nz,8) :: yi_bulk=1D0
	REAL(dp), DIMENSION(NSPEC_P) :: c_p_bulk

    !!!!!! for clusterin
    integer, parameter                                :: clustering_systems =3
    character(len=255)                                ::  filename 

    type(clustering_mod)::chem_1, chem_2, chem_3
    REAL(DP) :: c_dma , comp_evap(nz,NSPEC_P)
    integer :: n_vapor_syst
    CHARACTER(LEN=11),allocatable          :: names_vapor_syst_readin(:,:)
    integer:: kk
    real(dp) :: n_evap(nz), comp_evap_via_condensation(nz,NSPEC_P),n_evap_not_inuse(nz)
    real(dp), allocatable  :: Mx_chem1(:), Mx_chem2(:),Mx_chem3(:), qX_chem1(:), qx_chem2(:),qx_chem3(:)
    real(dp) ::  real_tmp1, real_tmp2, real_tmp3,  volume_chem1,volume_chem2,volume_chem3
    logical :: l_condensation_evap, l_coagulation_loss

    !!! for testing use constant hio2
    REAL(DP) :: cHIO2(nz)
    ! REAL(dp), ALLOCATABLE :: c_p_clust1(:,:),c_p_clust2(:,:), v_clust1(:),v_clust2(:)
    
    l_condensation_evap=.False.
    l_coagulation_loss=.True.

  
    comp_evap=0D0
    !! condensation update

    allocate(chem_1%names_vapor(n_clustering_vapors))
    allocate(chem_1%conc_vapor(n_clustering_vapors))
    allocate(chem_1%nmols_evap(n_clustering_vapors))
    
    allocate(chem_2%names_vapor(n_clustering_vapors))
    allocate(chem_2%conc_vapor(n_clustering_vapors))
    allocate(chem_2%nmols_evap(n_clustering_vapors))
    
    allocate(chem_3%names_vapor(n_clustering_vapors))
    allocate(chem_3%conc_vapor(n_clustering_vapors))
    allocate(chem_3%nmols_evap(n_clustering_vapors))
    
    allocate(chem_1%nconc_evap(nz))
    allocate(chem_2%nconc_evap(nz))
    allocate(chem_3%nconc_evap(nz))
    allocate(names_vapor_syst_readin(n_clustering_vapors,clustering_systems))
    
    allocate(Mx_chem1(n_clustering_vapors))
    allocate(Mx_chem2(n_clustering_vapors))
    allocate(Mx_chem3(n_clustering_vapors))
    
    allocate(qX_chem1(n_clustering_vapors))
    allocate(qX_chem2(n_clustering_vapors))
    allocate(qX_chem3(n_clustering_vapors))
    
    call get_system_size_1(neq_syst=chem_1%neq_syst)
    call get_system_size_2(neq_syst=chem_2%neq_syst)
    call get_system_size_3(neq_syst=chem_3%neq_syst)
    
    allocate(c_clusters_n(nz,chem_1%neq_syst))
    allocate(c_clusters_d(nz,chem_2%neq_syst))
    allocate(c_clusters_i(nz,chem_3%neq_syst))
    
    allocate(c_clusters1(chem_1%neq_syst))
    allocate(c_clusters2(chem_2%neq_syst))
    allocate(c_clusters3(chem_3%neq_syst))
    
    c_clusters1=0d0; c_clusters2=0d0; c_clusters3=0d0; c_clusters_n=0d0;c_clusters_d=0d0;c_clusters_i=0d0;
    n_evap=0D0; ;comp_evap=0D0; comp_evap_via_condensation=0D0
    
    use_clustering_plugin=.True.
    clust_evap=.True.
    clust_firstcall=.True.


    
    
!**************************************************************************************************************************************************************************
      !       Cl2, Cl,  ClO,ClO2, HCl,HOCl,ClNO,ClNO2,ClNO3,Br2,Br, BrO,HBr,HOBr,BrNO2,BrNO3,BrCl,I2,   I,   IO,  OIO, I2O2,HI,  HOI, HIO3,INO2
      MX_aq=(/71.,35.5,51.5,67.5,36.5,52.5,65.5,81.5, 97.5,160.,80.,96.,81.,97., 126.,142., 115.5,254.,127.,143.,159.,286.,128.,144.,176.,173.,&
             189.,162.5,207.,48.,17.,34. ,62.,78., 94.,  80., 108., 64.,33.,62.,100.,30., 46.,32.,80.,17.,45.,107.,63.,98.,97.,159.11/)*1D-3     
            !INO3,ICl,  IBr, O3, OH, H2O2,DMS,DMSO,DMSO2,MSIA,HPMTF,SO2,HO2,NO3,RO2, HCHO,NO2,O2,SO3,NH3,DMA,HOOCH2SCO,HNO3,H2SO4,MSA, HIO2

       ! 5.3D-4 Kreidenweis et al. [2003] O3
       D_aq=1D-5
       D_aq(30)=1.5D-5 ! O3
       D_aq(31)=2D-5 ! OH
       D_aq(32)=1.5D-5 ! H2O2
       D_aq(39)=1.5D-5 ! HO2
       D_aq(43)=1.5D-5 ! NO2
       D_aq(36)=1D-5 ! MSIA
       D_aq(37)=1D-5 ! HPMTF
       D_aq(46)=2.5D-5 ! NH3
	   D_aq(47)=1.5D-5 ! DMA
       D_aq(49)=(0.22/1.87)*1D-4 ! HNO3
       D_aq(50)=1D-5    ! H2SO4
       D_aq(51)=1D-5 ! MSA

       VX_aq=MX_aq/1500.0/Na             ! molecule volume (m^3 molec)
       dX_aq=(VX_aq*6D0/pi)**(1D0/3D0)   ! Molecule diameter (m)

       alpha_aq=1D-1 ! If not otherwise specified keep the mass accommodation coefficient at unity!
                     ! Cl2,  Cl,  ClO,ClO2,HCl, HOCl,ClNO,ClNO2,ClNO3
       alpha_aq(1:9)=(/0.08,0.05,0.064,0.05,0.1036,0.5,0.01,0.01,0.1/)
                         !Br2,Br, BrO, HBr, HOBr,BrNO2,BrNO3,BrCl
       alpha_aq(10:17)=(/0.08,0.05,0.06,0.0481,0.5,0.01,0.8,0.33/)
                          !I2,   I,   IO,  OIO, I2O2,HI,  HOI, HIO3,INO2,INO3,ICl,  IBr
       alpha_aq(18:29)=(/0.0126,0.05,0.558,1.0,0.123,0.057,0.5,0.0126,0.123,0.123,0.0126,0.0126/) ! Estimated for Iodine species (Bräuer et al., 2013)
       alpha_aq(30)=4D-4 ! O3
       alpha_aq(31)=0.05 ! OH
       alpha_aq(33)=0.001 ! DMS
       alpha_aq(50)=1D0   ! H2SO4
       alpha_aq(51)=1D0   ! MSA
       alpha_aq(52)=1D0   ! MSA
	   
	   
dz=z(2:Nz+1)-z(1:Nz)
dz2=(dz(2:Nz)+dz(1:Nz-1))/2D0


    CALL cpu_time(start_time)

    t_start = 0. ! start time (s)
    t_end = t_start + (days*24+hours)*3600 ! amount of s
    t_start_chem = (0D0*3600.) ! starting time chemistry (s)
!    t_start_chem = t_end ! starting time chemistry (s)
!    t_end = 60.
    dt = 3D1 ! model main time step (s) 
    EM_DT = dt ! Time step for the MEGAN model
    tr = 1       
	tr_1h = 1

    !--------------------------------------------------------------!
    ! Open output files                                            !
    !--------------------------------------------------------------!
    OPEN(300,FILE='adchem1d_namelist')
    READ(300,nml=adchem1d_simpar)
    ! these are saved every hour 
      OPEN(230,FILE='output/CS_H2SO4_st_'//sim//date,STATUS='REPLACE',ACTION='WRITE')
    OPEN(231,FILE='output/Jnucl_N_st_'//sim//date,STATUS='REPLACE',ACTION='WRITE')
    OPEN(300,FILE='output/Jnucl_D_st_'//sim//date,STATUS='REPLACE',ACTION='WRITE')
	OPEN(204,FILE='output/c_p_st_'//sim//date,STATUS='REPLACE',ACTION='WRITE')
	OPEN(212,FILE='output/conc_st_'//sim//date,STATUS='REPLACE',ACTION='WRITE')
    OPEN(213,FILE='output/N_bins_st_'//sim//date,STATUS='REPLACE',ACTION='WRITE')
    OPEN(215,FILE='output/c_p_tot_st_'//sim//date,STATUS='REPLACE',ACTION='WRITE')
    OPEN(216,FILE='output/dens_p_st_'//sim//date,STATUS='REPLACE',ACTION='WRITE')
    OPEN(217,FILE='output/psat_st_'//sim//date,STATUS='REPLACE',ACTION='WRITE')
    OPEN(222,FILE='output/d_p_dry_st_'//sim//date,STATUS='REPLACE',ACTION='WRITE')
    
	! These are saved only once at the beginning
    OPEN(218,FILE='output/index_cond_'//sim//date,STATUS='REPLACE',ACTION='WRITE')
    OPEN(220,FILE='output/molar_mass_'//sim//date,STATUS='REPLACE',ACTION='WRITE')
    OPEN(223,FILE='output/species_name_'//sim//date,STATUS='REPLACE',ACTION='WRITE')
    OPEN(232,FILE='output/RH_st_'//sim//date,STATUS='REPLACE',ACTION='WRITE')
    OPEN(233,FILE='output/TEMP_st_'//sim//date,STATUS='REPLACE',ACTION='WRITE')
    OPEN(234,FILE='output/pH_st_'//sim//date,STATUS='REPLACE',ACTION='WRITE')
    

 	
	OPEN(100,file='input/meteorology/cloudLWC_'//date)
    READ(100,*) cloudLWC ! g/m^3
	OPEN(101,file='input/meteorology/fricition_velocity_'//date)
    READ(101,*) fric_veloc ! m/s
	OPEN(105,file='input/meteorology/PBLH_'//date)
    READ(105,*) PBLH ! Planetary boundary layer height (m)
    OPEN(106,file='input/meteorology/surface_SHTF_'//date)
    READ(106,*) SHTFs ! Sensible heat flux at surface (W/m^2)
    OPEN(107,file='input/meteorology/precipitation_'//date)
    READ(107,*) rain ! rain (mm/h)
	!WHERE (rain>10.0) rain=10.0 ! Don't allow precipitation greater than 20 mm/h
	
	OPEN(108,file='input/meteorology/lat_'//date)
    READ(108,*) latitude ! temp in all layers (K)
    OPEN(109,file='input/meteorology/temp_'//date)
    READ(109,*) Ts ! temp in all layers (K)
    OPEN(110,file='input/meteorology/u_vel_'//date)
    READ(110,*) u_vel ! u velocity in all layers (m/s)
    OPEN(126,file='input/meteorology/v_vel_'//date)
    READ(126,*) v_vel ! v velocity in all layers (m/s)
    OPEN(127,file='input/meteorology/press_'//date)
    READ(127,*) Ps ! v velocity in all layers (m/s)
    OPEN(128,file='input/meteorology/RH_'//date)
    READ(128,*) RHs ! v velocity in all layers (m/s)
    OPEN(134,file='input/meteorology/surface_DSWF_v2_'//date)
    READ(134,*) DSWF ! downward short wave rad fux (W/m²)
	OPEN(137,file='input/meteorology/windspeed10m_'//date)
    READ(137,*) windspeed ! Wind speed at ~10 m
	OPEN(138,file='input/meteorology/SST_'//date)
    READ(138,*) SST ! Sea surface temp (K)
	
	
    CLOSE(100)
    CLOSE(101)
    CLOSE(102)
    CLOSE(103)
    CLOSE(104)
    CLOSE(105)
    CLOSE(106)
    CLOSE(107)
    CLOSE(108)
    CLOSE(109)
    CLOSE(110)
    CLOSE(126)
    CLOSE(127)
    CLOSE(128)
    CLOSE(129)
    CLOSE(134)
	CLOSE(137)
	CLOSE(138)

    RHs=RHs*1D2 ! % 

    DO i = 1,traj_len_min
      WHERE (RHs(:,i) > 99.9D0) RHs(:,i) = 99.9D0
	  WHERE (RHs(:,i) < 60.0) RHs(:,i) = 60.0
	  WHERE (cloudLWC(:,i)>1D-2) RHs(:,i) = 99.9D0 ! grid cell with clouds
	  WHERE (cloudLWC(:,i)>0.5) cloudLWC(:,i)=0.5 ! Limit the maximum cloud liquid water content
	END DO
    
	!DO i = 3,traj_len_min-2
	!RHs(:,i)=(RHs(:,i-2)+RHs(:,i-1)+RHs(:,i)+RHs(:,i+1)+RHs(:,i+2))/5D0 ! Moving average of RH to prevent rappid shift
	!END DO
	
	
	T2m=Ts(1,:)
	P2m=Ps(1,:)
	RH2m=RHs(1,:)
	!windspeed=sqrt(u_vel(2,:)**2D0+v_vel(2,:)**2D0) ! Wind speed at ~10 m



    ! load emission files and land-use category 
    OPEN(112,file='input/emissions/COemission_'//date)
    READ(112,*) EmissionCO ! emission of CO in (kg m^-2 s^-1)
	OPEN(113,file='input/emissions/NOxemission_'//date)
	READ(113,*) EmissionNO2 ! emission of NO2 in (kg m^-2 s^-1)
	OPEN(114,file='input/emissions/SO2emission_'//date)
    READ(114,*) EmissionSO2 ! emission of SO2 in (kg m^-2 s^-1) 
	OPEN(115,file='input/meteorology/land_cat_'//date)
    READ(115,*) landuse_index ! Used for dry deposition calculations
    OPEN(116,file='input/emissions/NH3emission_'//date)
    READ(116,*) EmissionNH3 ! emission of NH3 in (kg m^-2 s^-1)
	OPEN(117,file='input/emissions/NMVOCemission_'//date)
    READ(117,*) EmissionNMVOC ! emission of NMVOC in (kg m^-2 s^-1)
	OPEN(118,file='input/emissions/DMSemission_'//date)
    READ(118,*) EmissionDMS ! emission of DMS in (kg m^-2 s^-1)
	OPEN(129,file='input/emissions/NOxsoilemission_'//date)
    READ(129,*) EmissionNOxsoil ! emission of NO2 from soils in (kg m^-2 s^-1)
	OPEN(130,file='input/emissions/CH2Br2emission_'//date)
    READ(130,*) EmissionCH2Br2 ! (kg m^-2 s^-1)
	OPEN(131,file='input/emissions/CHBr3emission_'//date)
    READ(131,*) EmissionCHBr3 ! (kg m^-2 s^-1)
	OPEN(132,file='input/emissions/CH3Iemission_'//date)
    READ(132,*) EmissionCH3I ! (kg m^-2 s^-1)
	OPEN(133,file='input/emissions/NH3birdemission_'//date)
    READ(133,*) EmissionNH3_birds ! (kg m^-2 s^-1)
	OPEN(134,file='input/emissions/init_conc_v2_'//date)
    READ(134,*) init_conc ! Initial concentrations from CAMS
	IF (sea_spray_index == 2) THEN
	OPEN(135,file='input/emissions/Seaspray_Salter_'//date)
	ELSEIF (sea_spray_index == 3) THEN
	OPEN(135,file='input/emissions/Seaspray_Sofiev_'//date)
	ELSEIF (sea_spray_index == 4) THEN
	OPEN(135,file='input/emissions/Seaspray_Salter2_'//date)
	END IF
    READ(135,*) E_sea_salt_ 
	OPEN(139,file='input/emissions/Snow_salt_'//date)
	READ(139,*) E_snow_salt_ 
	

    !!! read ofr clusterin
    do kk=1, clustering_systems
   
        WRITE (filename, fmt='(a,I1a)') './ACDC_versions/Clusterin_multiple_chem_including_HIO3/cluster_chem_spec_',kk,'.txt'
        open(400,FILE=filename,action='read')
        read (400,*) n_vapor_syst
             i = 0
             do while (i .lt. n_vapor_syst)
                 i = i+1
                 read (400,*) names_vapor_syst_readin(i,kk)
              end do
              close(400)
              write(*,*) 'l567',' ', names_vapor_syst_readin(:,kk)
              
     end do
     
    chem_1%names_vapor=names_vapor_syst_readin(:,1)
    chem_2%names_vapor=names_vapor_syst_readin(:,2)
    chem_3%names_vapor=names_vapor_syst_readin(:,3)
    
    
    
    !  write(*,*) 'L572', ' ', chem_1%names_vapor
    !  write(*,*) 'L573', ' ', chem_2%names_vapor
		
	! Test to increase ship emissions x10
	EmissionSO2(7,:)=EmissionSO2(7,:)*1D1
	EmissionNO2(7,:)=EmissionNO2(7,:)*1D1
	EmissionNH3(7,:)=EmissionNH3(7,:)*1D1
	EmissionNMVOC(7,:)=EmissionNMVOC(7,:)*1D1
	EmissionNH3(7,:)=EmissionNH3(7,:)*1D1
	
	
	where (landuse_index==7) landuse_index=1 ! don't consider urban areas in dry dep. scheme!
    where (landuse_index==6 .AND. EmissionDMS<1D-13) landuse_index=8 ! Correct land use index over ocean to sea ice if the DMS emissions are very low
	
	
	
	!OPEN(139,file='input/emissions/Forest_fire_emission_'//date)
    !READ(139,*)  DM_emissions_biomass_burning     ! Dry biomass burned (g m^-2 h^-1)
    	
	CLOSE(112); CLOSE(113); CLOSE(114); CLOSE(115); CLOSE(116);CLOSE(117);CLOSE(118); CLOSE(129); CLOSE(130); 
	CLOSE(131); CLOSE(132); CLOSE(133); CLOSE(134); CLOSE(135)
	! Convert emissions from kg m^-2 s^-1 to molecules/cm^-2 h^-1:
	
	EmissionNOx_road=EmissionNO2(6,:)*1D3 ! g m^-2 s^-1
	EmissionNOx_nonroad_mach=EmissionNO2(9,:)*1D3 ! g m^-2 s^-1
	
	EmissionCO=EmissionCO*1E-4*1E3/28D0*Na*3600D0
	EmissionNH3=EmissionNH3*1E-4*1E3/17D0*Na*3600D0
	EmissionNO2=EmissionNO2*1E-4*1E3/46D0*Na*3600D0
	EmissionSO2=EmissionSO2*1E-4*1E3/64D0*Na*3600D0 
	EmissionNOxsoil=EmissionNOxsoil*1E-4*1E3/46D0*Na*3600D0

	EmissionSO2_ship=EmissionSO2(7,:) ! G=shipping
	!EmissionNO2(8,:)=EmissionNO2(8,:)*0.0 ! Don't consider emissions from Aviation
	!EmissionSO2(8,:)=EmissionSO2(8,:)*0.0 ! Don't consider emissions from Aviation
	EmissionDMS=EmissionDMS*1E-4*1E3/62D0*Na*3600D0
    EmissionCH2Br2=EmissionCH2Br2*1E-4*1E3/173.83D0*Na*3600D0
    EmissionCHBr3=EmissionCHBr3*1E-4*1E3/252.731D0*Na*3600D0 ! Corrected Pontus 20210105
    EmissionCH3I=EmissionCH3I*1E-4*1E3/141.94D0*Na*3600D0 ! Corrected Pontus 20210105
    
    ! Added Pontus 20210105:
    EmissionTCE=EmissionCHBr3*3.2D6/1.34D7 ! C2Cl4 cm^-2 s^-1, 3.2D6 Bräuer et al., 2013
    EmissionTRICLETH=EmissionCHBr3*4D6/1.34D7 ! C2HCL3, cm^-2 s^-1, 4D6 Bräuer et al., 2013
    EmissionCHCL3=EmissionCHBr3*6.4D7/1.34D7 ! cm^-2 s^-1, 6.4D7 Bräuer et al., 2013
    EmissionCH2CL2=EmissionCHBr3*3.2D7/1.34D7 ! cm^-2 s^-1, 3.2D7 Bräuer et al., 2013
    EmissionCH3CL=EmissionCHBr3*9.1D7/1.34D7  ! cm^-2 s^-1, 9.1D7 Bräuer et al., 2013

    EmissionC3H7I=EmissionCH3I*4.92D5/1.48D7 ! cm^-2 s^-1, Scaling based on Hoffmann, et al PNAS 2016 supplement
    EmissionC2H5I=EmissionCH3I*1.66D6/1.48D7 ! cm^-2 s^-1, Scaling based on Hoffmann, et al PNAS 2016 supplement
    EmissionCH2I2=EmissionCH3I*1.13D7/1.48D7 ! cm^-2 s^-1, Scaling based on Hoffmann, et al PNAS 2016 supplement
    EmissionCH2ICL=EmissionCH3I*9.18D6/1.48D7 ! cm^-2 s^-1, Scaling based on Hoffmann, et al PNAS 2016 supplement
    EmissionCH2IBr=EmissionCH3I*5.27D6/1.48D7 ! cm^-2 s^-1, Scaling based on Hoffmann, et al PNAS 2016 supplement
    EmissionDIBRET=EmissionCHBr3*4.23D6/1.34D7  ! cm^-2 s^-1, Scaling based on Hoffmann, et al PNAS 2016 supplement
    EmissionCH3BR=EmissionCHBr3*9.7D6/1.34D7  ! cm^-2 s^-1, Scaling based on Hoffmann, et al PNAS 2016 supplement
    
    EmissionC2H6=EmissionDMS*1D7/6.18D9 ! cm^-2 s^-1,  Scaling based on Bräuer et al., 2013
    EmissionC3H8=EmissionDMS*1.6D7/6.18D9 ! cm^-2 s^-1, !   Scaling based on Hoffmann et al. 2016
    EmissionCH3CHO=EmissionDMS*3.6D9/6.18D9 ! Acetaldehyde,  Scaling based on Bräuer et al., 2013
    EmissionC2H5CHO=EmissionDMS*5.47D9/6.18D9 ! Propionaldehyde,  Scaling based on Bräuer et al., 2013
    EmissionISO=EmissionDMS*1.41D8/6.18D9 ! Isoprene  Scaling based on Hoffmann et al. 2016
    
	
	
                                !   C2H6   NC4H10 C2H4  C3H6   C5H8  OXYL   CH3OH C2H5OH HCHO   CH3CHO MEK  GLYOX MGLYOX UNREAC
   VOC_dist_A_PublicPower=       (/12.589, 39.790,8.174,10.767,0.000,18.632,0.000,3.912, 5.586, 0.207,0.089,0.000,0.000,0.255/)/1D2 ! A_PublicPower = S2
   VOC_dist_B_Industry=           (/4.996, 35.610,9.044,2.089, 0.000,18.323,0.561,3.034, 24.134,0.059,1.347,0.000,0.000,0.805/)/1D2 ! B_Industry=S3
   VOC_dist_C_OtherStationaryComb=(/12.589,39.790,8.174,10.767,0.000,18.632,0.000,3.912, 5.586, 0.207,0.089,0.000,0.000,0.255/)/1D2 !C_OtherStationaryComb = S2
   VOC_dist_D_Fugitiv=            (/0.444, 44.052,0.244,0.678, 0.008,17.904,6.101,16.416,0.011, 0.000,9.965,0.000,0.000,4.176/)/1D2 !D_Fugitive=S6
   VOC_dist_E_Solvents=           (/0.444, 44.052,0.244,0.678, 0.008,17.904,6.101,16.416,0.011, 0.000,9.965,0.000,0.000,4.176/)/1D2 !E_Solvents = S6
   VOC_dist_F_RoadTransport=      (/4.832, 36.698,6.796,10.896,0.000,35.051,0.000,0.000, 2.700, 2.606,0.421,0.000,0.000,0.000/)/1D2 !F_RoadTransport = S7
   VOC_dist_G_Shipping =          (/3.775, 47.416,6.636,10.608,0.000,24.676,0.000,0.000, 3.115, 3.261,0.235,0.146,0.117,0.014/)/1D2 ! G_Shipping 
   VOC_dist_I_Offroad=            (/ 3.775,47.42,6.636,10.61,  0.000,24.68,0.000,0.000,  3.115, 3.261,0.235,0.146,0.117,0.014/)/1D2 !I_Offroad =S8
   VOC_dist_J_Waste=              (/ 25.72,36.78,5.237,1.830,  1.153,7.881,0.427,2.439,  16.060,0.000,0.093,0.000,0.000,2.383/)/1D2 !J_Waste = S9
   VOC_dist_K_AgriLivestock=      (/ 0.000,0.000,0.000,0.000,  0.000,0.000,0.000,0.000,  0.000, 0.000,0.000,0.000,0.000,100.000/)/1D2 !K_AgriLivestock=S10
   VOC_dist_L_AgriOther=          (/ 0.000,0.000,0.000,0.000,  0.000,0.000,0.000,0.000,  0.000, 0.000,0.000,0.000,0.000,100.000/)/1D2 !L_AgriOther=S10
   
   VOC_dist_M_Other =             VOC_dist_C_OtherStationaryComb ! M_Other=S11
   VOC_dist_H_Aviation=0D0 ! H_Aviation=Aviation
   
   EmissionC2H6=EmissionNMVOC(4,:)*VOC_dist_D_Fugitiv(1)+&
				EmissionNMVOC(5,:)*VOC_dist_E_Solvents(1)+EmissionNMVOC(6,:)*VOC_dist_F_RoadTransport(1)+&
				EmissionNMVOC(8,:)*VOC_dist_H_Aviation(1)+&
				EmissionNMVOC(9,:)*VOC_dist_I_Offroad(1)+EmissionNMVOC(10,:)*VOC_dist_J_Waste(1)+&
				EmissionNMVOC(11,:)*VOC_dist_K_AgriLivestock(1)+EmissionNMVOC(12,:)*VOC_dist_L_AgriOther(1)
   
   EmissionC2H6_Pplants_Ind=EmissionNMVOC(1,:)*VOC_dist_A_PublicPower(1)+EmissionNMVOC(2,:)*VOC_dist_B_Industry(1)+&
                EmissionNMVOC(3,:)*VOC_dist_C_OtherStationaryComb(1)
	
   EmissionC2H6_Ship=EmissionNMVOC(7,:)*VOC_dist_G_Shipping(1)
	
   EmissionNC4H10=EmissionNMVOC(4,:)*VOC_dist_D_Fugitiv(2)+&
				EmissionNMVOC(5,:)*VOC_dist_E_Solvents(2)+EmissionNMVOC(6,:)*VOC_dist_F_RoadTransport(2)+&
				EmissionNMVOC(8,:)*VOC_dist_H_Aviation(2)+&
				EmissionNMVOC(9,:)*VOC_dist_I_Offroad(2)+EmissionNMVOC(10,:)*VOC_dist_J_Waste(2)+&
				EmissionNMVOC(11,:)*VOC_dist_K_AgriLivestock(2)+EmissionNMVOC(12,:)*VOC_dist_L_AgriOther(2)
				
   EmissionNC4H10_Pplants_Ind=EmissionNMVOC(1,:)*VOC_dist_A_PublicPower(2)+EmissionNMVOC(2,:)*VOC_dist_B_Industry(2)+&
                EmissionNMVOC(3,:)*VOC_dist_C_OtherStationaryComb(2)
   
   EmissionNC4H10_Ship=EmissionNMVOC(7,:)*VOC_dist_G_Shipping(2)  
   
   EmissionC2H4=EmissionNMVOC(4,:)*VOC_dist_D_Fugitiv(3)+&
				EmissionNMVOC(5,:)*VOC_dist_E_Solvents(3)+EmissionNMVOC(6,:)*VOC_dist_F_RoadTransport(3)+&
				EmissionNMVOC(8,:)*VOC_dist_H_Aviation(3)+&
				EmissionNMVOC(9,:)*VOC_dist_I_Offroad(3)+EmissionNMVOC(10,:)*VOC_dist_J_Waste(3)+&
				EmissionNMVOC(11,:)*VOC_dist_K_AgriLivestock(3)+EmissionNMVOC(12,:)*VOC_dist_L_AgriOther(3)
   EmissionC2H4_Pplants_Ind=EmissionNMVOC(1,:)*VOC_dist_A_PublicPower(3)+EmissionNMVOC(2,:)*VOC_dist_B_Industry(3)+&
                EmissionNMVOC(3,:)*VOC_dist_C_OtherStationaryComb(3)
   
   EmissionC2H4_Ship=EmissionNMVOC(7,:)*VOC_dist_G_Shipping(3)
   
   EmissionC3H6=EmissionNMVOC(4,:)*VOC_dist_D_Fugitiv(4)+&
				EmissionNMVOC(5,:)*VOC_dist_E_Solvents(4)+EmissionNMVOC(6,:)*VOC_dist_F_RoadTransport(4)+&
				EmissionNMVOC(8,:)*VOC_dist_H_Aviation(4)+&
				EmissionNMVOC(9,:)*VOC_dist_I_Offroad(4)+EmissionNMVOC(10,:)*VOC_dist_J_Waste(4)+&
				EmissionNMVOC(11,:)*VOC_dist_K_AgriLivestock(4)+EmissionNMVOC(12,:)*VOC_dist_L_AgriOther(4)
   EmissionC3H6_Pplants_Ind=EmissionNMVOC(1,:)*VOC_dist_A_PublicPower(4)+EmissionNMVOC(2,:)*VOC_dist_B_Industry(4)+&
                EmissionNMVOC(3,:)*VOC_dist_C_OtherStationaryComb(4)
   
   EmissionC3H6_Ship=EmissionNMVOC(7,:)*VOC_dist_G_Shipping(4)
   
   EmissionC5H8=EmissionNMVOC(4,:)*VOC_dist_D_Fugitiv(5)+&
				EmissionNMVOC(5,:)*VOC_dist_E_Solvents(5)+EmissionNMVOC(6,:)*VOC_dist_F_RoadTransport(5)+&
				EmissionNMVOC(8,:)*VOC_dist_H_Aviation(5)+&
				EmissionNMVOC(9,:)*VOC_dist_I_Offroad(5)+EmissionNMVOC(10,:)*VOC_dist_J_Waste(5)+&
				EmissionNMVOC(11,:)*VOC_dist_K_AgriLivestock(5)+EmissionNMVOC(12,:)*VOC_dist_L_AgriOther(5)
   EmissionC5H8_Pplants_Ind=EmissionNMVOC(1,:)*VOC_dist_A_PublicPower(5)+EmissionNMVOC(2,:)*VOC_dist_B_Industry(5)+&
                EmissionNMVOC(3,:)*VOC_dist_C_OtherStationaryComb(5)				

   EmissionC5H8_Ship=EmissionNMVOC(7,:)*VOC_dist_G_Shipping(5)
   
   EmissionOXYL=EmissionNMVOC(4,:)*VOC_dist_D_Fugitiv(6)+&
				EmissionNMVOC(5,:)*VOC_dist_E_Solvents(6)+EmissionNMVOC(6,:)*VOC_dist_F_RoadTransport(6)+&
				EmissionNMVOC(8,:)*VOC_dist_H_Aviation(6)+&
				EmissionNMVOC(9,:)*VOC_dist_I_Offroad(6)+EmissionNMVOC(10,:)*VOC_dist_J_Waste(6)+&
				EmissionNMVOC(11,:)*VOC_dist_K_AgriLivestock(6)+EmissionNMVOC(12,:)*VOC_dist_L_AgriOther(6)
				
   EmissionOXYL_Pplants_Ind=EmissionNMVOC(1,:)*VOC_dist_A_PublicPower(6)+EmissionNMVOC(2,:)*VOC_dist_B_Industry(6)+&
                EmissionNMVOC(3,:)*VOC_dist_C_OtherStationaryComb(6)				

   EmissionOXYL_Ship=EmissionNMVOC(7,:)*VOC_dist_G_Shipping(6)
   
   EmissionCH3OH=EmissionNMVOC(4,:)*VOC_dist_D_Fugitiv(7)+&
				EmissionNMVOC(5,:)*VOC_dist_E_Solvents(7)+EmissionNMVOC(6,:)*VOC_dist_F_RoadTransport(7)+&
				EmissionNMVOC(8,:)*VOC_dist_H_Aviation(7)+&
				EmissionNMVOC(9,:)*VOC_dist_I_Offroad(7)+EmissionNMVOC(10,:)*VOC_dist_J_Waste(7)+&
				EmissionNMVOC(11,:)*VOC_dist_K_AgriLivestock(7)+EmissionNMVOC(12,:)*VOC_dist_L_AgriOther(7)
   
   EmissionCH3OH_Pplants_Ind=EmissionNMVOC(1,:)*VOC_dist_A_PublicPower(7)+EmissionNMVOC(2,:)*VOC_dist_B_Industry(7)+&
                EmissionNMVOC(3,:)*VOC_dist_C_OtherStationaryComb(7)				

   EmissionCH3OH_Ship=EmissionNMVOC(7,:)*VOC_dist_G_Shipping(7)
   
   EmissionC2H5OH=EmissionNMVOC(4,:)*VOC_dist_D_Fugitiv(8)+&
				EmissionNMVOC(5,:)*VOC_dist_E_Solvents(8)+EmissionNMVOC(6,:)*VOC_dist_F_RoadTransport(8)+&
				EmissionNMVOC(8,:)*VOC_dist_H_Aviation(8)+&
				EmissionNMVOC(9,:)*VOC_dist_I_Offroad(8)+EmissionNMVOC(10,:)*VOC_dist_J_Waste(8)+&
				EmissionNMVOC(11,:)*VOC_dist_K_AgriLivestock(8)+EmissionNMVOC(12,:)*VOC_dist_L_AgriOther(8)
   
   EmissionC2H5OH_Pplants_Ind=EmissionNMVOC(1,:)*VOC_dist_A_PublicPower(8)+EmissionNMVOC(2,:)*VOC_dist_B_Industry(8)+&
                EmissionNMVOC(3,:)*VOC_dist_C_OtherStationaryComb(8)				
	
   EmissionC2H5OH_Ship=EmissionNMVOC(7,:)*VOC_dist_G_Shipping(8)
   	
   EmissionHCHO=EmissionNMVOC(4,:)*VOC_dist_D_Fugitiv(9)+&
				EmissionNMVOC(5,:)*VOC_dist_E_Solvents(9)+EmissionNMVOC(6,:)*VOC_dist_F_RoadTransport(9)+&
				EmissionNMVOC(8,:)*VOC_dist_H_Aviation(9)+&
				EmissionNMVOC(9,:)*VOC_dist_I_Offroad(9)+EmissionNMVOC(10,:)*VOC_dist_J_Waste(9)+&
				EmissionNMVOC(11,:)*VOC_dist_K_AgriLivestock(9)+EmissionNMVOC(12,:)*VOC_dist_L_AgriOther(9)
   
   EmissionHCHO_Pplants_Ind=EmissionNMVOC(1,:)*VOC_dist_A_PublicPower(9)+EmissionNMVOC(2,:)*VOC_dist_B_Industry(9)+&
                EmissionNMVOC(3,:)*VOC_dist_C_OtherStationaryComb(9)
	
   EmissionHCHO_Ship=EmissionNMVOC(7,:)*VOC_dist_G_Shipping(9)
   	
   EmissionCH3CHO=EmissionNMVOC(4,:)*VOC_dist_D_Fugitiv(10)+&
				EmissionNMVOC(5,:)*VOC_dist_E_Solvents(10)+EmissionNMVOC(6,:)*VOC_dist_F_RoadTransport(10)+&
				EmissionNMVOC(8,:)*VOC_dist_H_Aviation(10)+&
				EmissionNMVOC(9,:)*VOC_dist_I_Offroad(10)+EmissionNMVOC(10,:)*VOC_dist_J_Waste(10)+&
				EmissionNMVOC(11,:)*VOC_dist_K_AgriLivestock(10)+EmissionNMVOC(12,:)*VOC_dist_L_AgriOther(10)
   
   
   EmissionCH3CHO_Pplants_Ind=EmissionNMVOC(1,:)*VOC_dist_A_PublicPower(10)+EmissionNMVOC(2,:)*VOC_dist_B_Industry(10)+&
                EmissionNMVOC(3,:)*VOC_dist_C_OtherStationaryComb(10)			
   
   EmissionCH3CHO_Ship=EmissionNMVOC(7,:)*VOC_dist_G_Shipping(10)
   			
   EmissionMEK=EmissionNMVOC(4,:)*VOC_dist_D_Fugitiv(11)+&
				EmissionNMVOC(5,:)*VOC_dist_E_Solvents(11)+EmissionNMVOC(6,:)*VOC_dist_F_RoadTransport(11)+&
				EmissionNMVOC(8,:)*VOC_dist_H_Aviation(11)+&
				EmissionNMVOC(9,:)*VOC_dist_I_Offroad(11)+EmissionNMVOC(10,:)*VOC_dist_J_Waste(11)+&
				EmissionNMVOC(11,:)*VOC_dist_K_AgriLivestock(11)+EmissionNMVOC(12,:)*VOC_dist_L_AgriOther(11)
   
   EmissionMEK_Pplants_Ind=EmissionNMVOC(1,:)*VOC_dist_A_PublicPower(11)+EmissionNMVOC(2,:)*VOC_dist_B_Industry(11)+&
                EmissionNMVOC(3,:)*VOC_dist_C_OtherStationaryComb(11)
   
   EmissionMEK_Ship=EmissionNMVOC(7,:)*VOC_dist_G_Shipping(11)
   		
   EmissionGLYOX=EmissionNMVOC(4,:)*VOC_dist_D_Fugitiv(12)+&
				EmissionNMVOC(5,:)*VOC_dist_E_Solvents(12)+EmissionNMVOC(6,:)*VOC_dist_F_RoadTransport(12)+&
				EmissionNMVOC(8,:)*VOC_dist_H_Aviation(12)+&
				EmissionNMVOC(9,:)*VOC_dist_I_Offroad(12)+EmissionNMVOC(10,:)*VOC_dist_J_Waste(12)+&
				EmissionNMVOC(11,:)*VOC_dist_K_AgriLivestock(12)+EmissionNMVOC(12,:)*VOC_dist_L_AgriOther(12)
   
   EmissionGLYOX_Pplants_Ind=EmissionNMVOC(1,:)*VOC_dist_A_PublicPower(12)+EmissionNMVOC(2,:)*VOC_dist_B_Industry(12)+&
                EmissionNMVOC(3,:)*VOC_dist_C_OtherStationaryComb(12)			
	
   EmissionGLYOX_Ship=EmissionNMVOC(7,:)*VOC_dist_G_Shipping(12)
   	
   EmissionMGLYOX=EmissionNMVOC(4,:)*VOC_dist_D_Fugitiv(13)+&
				EmissionNMVOC(5,:)*VOC_dist_E_Solvents(13)+EmissionNMVOC(6,:)*VOC_dist_F_RoadTransport(13)+&
				EmissionNMVOC(8,:)*VOC_dist_H_Aviation(13)+&
				EmissionNMVOC(9,:)*VOC_dist_I_Offroad(13)+EmissionNMVOC(10,:)*VOC_dist_J_Waste(13)+&
				EmissionNMVOC(11,:)*VOC_dist_K_AgriLivestock(13)+EmissionNMVOC(12,:)*VOC_dist_L_AgriOther(13)

   EmissionMGLYOX_Pplants_Ind=EmissionNMVOC(1,:)*VOC_dist_A_PublicPower(13)+EmissionNMVOC(2,:)*VOC_dist_B_Industry(13)+&
                EmissionNMVOC(3,:)*VOC_dist_C_OtherStationaryComb(13)
   
   EmissionMGLYOX_Ship=EmissionNMVOC(7,:)*VOC_dist_G_Shipping(13)
   				
   ! Convert emissions from kg m^-2 s^-1 to molecules cm^-2 h^-1:
   EmissionC2H6=EmissionC2H6*1E-4*1E3/30D0*3600D0*Na
   EmissionNC4H10=EmissionNC4H10*1E-4*1E3/58D0*3600D0*Na  
   EmissionC2H4=EmissionC2H4*1E-4*1E3/26D0*3600D0*Na
   EmissionC3H6=EmissionC3H6*1E-4*1E3/42D0*3600D0*Na
   EmissionC5H8=EmissionC5H8*1E-4*1E3/68D0*3600D0*Na
   EmissionOXYL=EmissionOXYL*1E-4*1E3/106D0*3600D0*Na
   EmissionCH3OH=EmissionCH3OH*1E-4*1E3/32D0*3600D0*Na
   EmissionC2H5OH=EmissionC2H5OH*1E-4*1E3/46D0*3600D0*Na
   EmissionHCHO=EmissionHCHO*1E-4*1E3/30D0*3600D0*Na
   EmissionCH3CHO=EmissionCH3CHO*1E-4*1E3/44D0*3600D0*Na
   EmissionMEK=EmissionMEK*1E-4*1E3/72D0*3600D0*Na
   EmissionGLYOX=EmissionGLYOX*1E-4*1E3/58D0*3600D0*Na
   EmissionMGLYOX=EmissionMGLYOX*1E-4*1E3/72D0*3600D0*Na
   EmissionC9_C12_alkanes=0.2*EmissionNC4H10
   
   EmissionC2H6_Pplants_Ind=EmissionC2H6_Pplants_Ind*1E-4*1E3/30D0*3600D0*Na
   EmissionNC4H10_Pplants_Ind=EmissionNC4H10_Pplants_Ind*1E-4*1E3/58D0*3600D0*Na  
   EmissionC2H4_Pplants_Ind=EmissionC2H4_Pplants_Ind*1E-4*1E3/26D0*3600D0*Na
   EmissionC3H6_Pplants_Ind=EmissionC3H6_Pplants_Ind*1E-4*1E3/42D0*3600D0*Na
   EmissionC5H8_Pplants_Ind=EmissionC5H8_Pplants_Ind*1E-4*1E3/68D0*3600D0*Na
   EmissionOXYL_Pplants_Ind=EmissionOXYL_Pplants_Ind*1E-4*1E3/106D0*3600D0*Na
   EmissionCH3OH_Pplants_Ind=EmissionCH3OH_Pplants_Ind*1E-4*1E3/32D0*3600D0*Na
   EmissionC2H5OH_Pplants_Ind=EmissionC2H5OH_Pplants_Ind*1E-4*1E3/46D0*3600D0*Na
   EmissionHCHO_Pplants_Ind=EmissionHCHO_Pplants_Ind*1E-4*1E3/30D0*3600D0*Na
   EmissionCH3CHO_Pplants_Ind=EmissionCH3CHO_Pplants_Ind*1E-4*1E3/44D0*3600D0*Na
   EmissionMEK_Pplants_Ind=EmissionMEK_Pplants_Ind*1E-4*1E3/72D0*3600D0*Na
   EmissionGLYOX_Pplants_Ind=EmissionGLYOX_Pplants_Ind*1E-4*1E3/58D0*3600D0*Na
   EmissionMGLYOX_Pplants_Ind=EmissionMGLYOX_Pplants_Ind*1E-4*1E3/72D0*3600D0*Na
   
!	OPEN(137,file='input/emissions/Primary_PN_emissions_'//date)
!   READ(137,*) PN_emissions_GAINS ! PN emissions in # m^-2 h^-1
!	PN_emissions_GAINS(5,:)=PN_emissions_GAINS(5,:)*0.1 ! Correct down power plant primary particle emissions
!	PN_emissions_GAINS=PN_emissions_GAINS*1D0/3600D0/dz(1) ! PN emissions in # m^-3 s^-1
	
!    OPEN(138,file='input/emissions/Primary_PN_emissions_dist_'//date)
!    READ(138,*)  PN_emissions_dist_param_GAINS     		
!	CLOSE(137); CLOSE(138)

IF (index_MEGAN==1) THEN
!    OPEN(122,file='input/biogenic_emissions/input_MEGAN_'//date)
!    READ(122,*) input_MEGAN 
!    CLOSE(122)
    
ELSE

    ! load emission files from MatLab program BVOC_emission.m (biogenic emission from LPJ-GUESS)
 !   OPEN(119,file='input/biogenic_emissions/apin_'//date)
 !   READ(119,*) EmissionAPIN ! emission of alpha-pinene in (molecules/cm^2h)
 !   OPEN(120,file='input/biogenic_emissions/bpin_'//date)
 !   READ(120,*) EmissionBPIN ! emission of beta-pinene in (molecules/cm^2h)
 !   OPEN(121,file='input/biogenic_emissions/lim_'//date)
 !   READ(121,*) EmissionLIM ! emission of limonene in (molecules/cm^2h)
 !   OPEN(133,file='input/biogenic_emissions/oth_mon_'//date)
 !   READ(133,*) EmissionOTH_MT ! emission of other monoterpenes (molecules/cm^2h)
 !   OPEN(124,file='input/biogenic_emissions/isoprene_'//date)
 !   READ(124,*) EmissionISO ! emission of other sesqiterpenes (molecules/cm^2h)
 !   OPEN(125,file='input/biogenic_emissions/all_bvoc_'//date)
 !   READ(125,*) Emission_all_bvoc ! emission of other sesqiterpenes (molecules/cm^2h)
 !   CLOSE(119)
 !   CLOSE(120)
 !   CLOSE(121)
 !   CLOSE(133)
 !   CLOSE(124)
 !   CLOSE(125) 
END IF

	JD=INT(input_MEGAN(4,1)) ! Julian day at the start of the simulation: 
	

OPEN(111,file = 'input/compounds/MCM33_allvoc_ELVOC_names_20180912.dat') ! mcm-name of each condensable compound
READ(111,'(A20)') mcm_name
OPEN(123,file = 'input/compounds/MCM33_allvoc_ELVOC_compound_prop_20180912.dat') ! column 1: molar mass (g/mol), column 2 and 3: Nannoolal coeff. (a and b resp)
READ(123,*) mcm_prop
OPEN(136,file = 'input/compounds/Activity_coeff_org_water_20180913.dat')
READ(136,*) y_org_water
CLOSE(111)
CLOSE(123)
CLOSE(136)

y_org_water=1D0

! Match the selected species from the MCM subset we want to have in the condensation module with all species in the MCM (from second_Monitor) to find what index the wanted species have in the gas-phase vector
    DO i = 1,NSPEC
       WRITE(223,*) SPC_NAMES(i)
       DO j = 1,NCOND
        IF (mcm_name(j) .eq. SPC_NAMES(i)) THEN ! SPC_NAMES from second-Monitor (same order as second_Parameters)
            index_cond(j) = i
        END IF
       END DO
    END DO      
	
    A_Nannoolal = mcm_prop(2,:) 
    B_Nannoolal = mcm_prop(3,:)
    NrC=mcm_prop(4,:)
    NrO=mcm_prop(7,:)
    NrH=mcm_prop(5,:)
    NrN=mcm_prop(6,:)
	NrS=NrN*0D0
	NrS(1)=1D0 ! MSA
    VolX=NrC*15.9D0+NrO*6.11D0+NrH*2.31D0+NrN*4.54D0+NrS*22.9D0 ! dimensionless diffusion volumes ! dimensionless diffusion volumes
    
    ! Initial dimer C, O, N, H atomic composition:
    dimer_C=0D0
    dimer_O=0D0
    dimer_N=0D0
    dimer_H=0D0    
    Morg = mcm_prop(1,:)*1D-3  ! Molar mass (kg/mol) of all condensable organic compounds

    MX = (/ MSO4,MNO3,MCl,MNH4,MNa,MEC,MH2O,MH2O,MMSA,MHIO3,MDMA,MHIO2,Morg(1:NCOND)/) ! (kg/mol) Vector with mole masses of all species in the particle phase
    qX(1:8)=(/ qSO4,qNO3,qCl,qNH4,qNa,qEC,qH2O,qH2O /) ! (kg/m^3) Vector with species specific densities in the particle phase
    qX(9)=qMSA ! MSA
	qX(10)=qHIO3 ! HIO3
	qX(11)=qDMA ! DMA
	qX(12)=qHIO2 ! HIO3
	qX(13:NSPEC_P)=qHC
    VX=MX/qX/Na*1D6 ! (cm^3 molec) Molecule volumes 
    molec_radius=1./2.*(VX*6./pi)**(1./3.)*1D-2;    ! (m) Estimated radius of molecules
    viscosity = 1D5 ! (Pa s) Initial particle viscosity (add particle phase water content dependence (RH dependence))  
    
    aX=1.                       ! Mass accommodation coefficients of condensable compounds
!	aX(2)=1d-1 ! HNO3 
!	aX(3)=1d-1 ! HCl Schweitzer et al. (2000)
!	aX(10)=0.0126D0 ! HIO3 Pechtl et al. (2005) 
	
    y=1D0
    pH=5.6

    !----------------------------------------------------------!
    ! Initialization (to be moved to own module later?)        !
    !----------------------------------------------------------!

  
    ! Particle initialization (1.07 nm growth ACDC)
    !d(1)=MAXVAL((/ 1.022D-9,delta_surf*2+1D-11 /)) ! Minimum size set by kinetic multi-layer model  
    d(1)=1.00001D-9 !1.022D-9!MAXVAL((/ 1.022D-9,delta_surf*2+1D-11 /)) ! Minimum size set by kinetic multi-layer model  
	vp(1)=(pi*d(1)**3)/6
    DO i=2,nr_bins+1
       vp(i)=vp(i-1)*1.32!1.265!1.31922 1.36!
    END DO
    
    ! 3 nm growth in ACDC:
    !d(1)=MAXVAL((/ 2.90D-9,delta_surf*2+1D-11 /)) ! Minimum size set by kinetic multi-layer model  
    !vp(1)=(pi*d(1)**3)/6
    !DO i=2,nr_bins+1
    !   vp(i)=vp(i-1)*1.226
    !END DO
    d = (vp*6./pi)**(1D0/3D0)
    dp_dry(1,:) = (d(1:nr_bins)+d(2:nr_bins+1))/2D0  ! Arithmetic mean diameter in each
    d_g = (d(1:nr_bins)*d(2:nr_bins+1))**0.5;  ! Geometric mean diameter in each size bin
	dp_dry1=dp_dry(1,:)
	dlogDp=log10(d(2:Nr_bins+1))-log10(d(1:nr_bins))	
    v_p = dp_dry(1,:)**3*pi/6D0    ! Single particle volume

	PM_BC_init=init_conc(13,:); ! TEST to decrease BC
	PM_OC_init=init_conc(14,:); ! TEST to decrease OA
	PM_SO4_init=init_conc(15,:);
	PM_NaCl1_init=init_conc(16,:); ! <0.5 um in diameter 
	PM_NaCl2_init=init_conc(17,:); ! >0.5 <5 um in diameter
	PM_NaCl3_init=init_conc(18,:); ! > 5 um in diameter
	
	! Limit the initial PM10: 
	where (PM_BC_init>1D-10) PM_BC_init=1D-10
	where (PM_OC_init>1D-9) PM_OC_init=1D-9
	where (PM_SO4_init>1D-8) PM_SO4_init=1D-8
	where (PM_NaCl1_init>1D-9) PM_NaCl1_init=1D-9
	where (PM_NaCl2_init>1D-8) PM_NaCl2_init=1D-8	
    where (PM_NaCl3_init>1D-8) PM_NaCl3_init=1D-8
		
	PM_NH4_init=PM_SO4_init*(5D-1*MX(4))/MX(1);
	PM_Na_init=PM_SO4_init*(5D-1*MX(5))/MX(1);
	
	! Specify the initial BC size distribution 
	N_modes3 = (/ 2D1, 1D1, 0D0 /)
    dm3 = (/ 6D-8, 1.2D-7, 1D-6 /)
    s3 = (/ 1.6, 1.6, 1.6 /)
	CALL dNdlogDp(d_g,dp_dry1,dlogDp,modes3,s3,N_modes3,dm3,Nbins_BC,Vbins_BC) 
    Nbins_BC=Nbins_BC/SUM(Nbins_BC); ! Normalize
	spm_BC=sum(Nbins_BC*v_p)*qEC! Mean single particle mass (kg) 
	PN_BC=PM_BC_init/spm_BC! m^-3 s^-1
	
	! Specify the initial OC size distribution 
	N_modes3 = (/ 2D1, 1D1, 1D-4 /)
    dm3 = (/ 6D-8, 1.2D-7, 7D-6 /)
    s3 = (/ 1.6, 1.6, 1.6 /)
	CALL dNdlogDp(d_g,dp_dry1,dlogDp,modes3,s3,N_modes3,dm3,Nbins_OC,Vbins_OC) 
    Nbins_OC=Nbins_OC/SUM(Nbins_OC); ! Normalize
	spm_OC=sum(Nbins_OC*v_p)*qHC! Mean single particle mass (kg) 
	PN_OC=PM_OC_init/spm_OC! m^-3 
	
	! Specify the initial SO4 size distribution 
	N_modes3 = (/ 2D1, 1D1, 1D-4 /)
    dm3 = (/ 6D-8, 1.2D-7, 7D-6 /)
    s3 = (/ 1.6, 1.6, 1.6 /)
	CALL dNdlogDp(d_g,dp_dry1,dlogDp,modes3,s3,N_modes3,dm3,Nbins_SO4,Vbins_SO4) 
    Nbins_SO4=Nbins_SO4/SUM(Nbins_SO4); ! Normalize
	spm_SO4=sum(Nbins_SO4*v_p)*qSO4 ! Mean single particle mass (kg) 
	PN_SO4=PM_SO4_init/spm_SO4! m^-3 
    
	! Specify the initial NH4 size distribution 
	N_modes3 = (/ 2D1, 1D1, 1D-4 /)
    dm3 = (/ 6D-8, 1.2D-7, 7D-6 /)
    s3 = (/ 1.6, 1.6, 1.6 /)
	CALL dNdlogDp(d_g,dp_dry1,dlogDp,modes3,s3,N_modes3,dm3,Nbins_NH4,Vbins_NH4) 
    Nbins_NH4=Nbins_NH4/SUM(Nbins_NH4); ! Normalize
	spm_NH4=sum(Nbins_NH4*v_p)*qNH4 ! Mean single particle mass (kg) 
	PN_NH4=PM_NH4_init/spm_NH4! m^-3 
	
	! Specify the initial Na+ size distribution connected to HSO4
	N_modes3 = (/ 2D1, 1D1, 1D-4 /)
    dm3 = (/ 6D-8, 1.2D-7, 7D-6 /)
    s3 = (/ 1.6, 1.6, 1.6 /)
	CALL dNdlogDp(d_g,dp_dry1,dlogDp,modes3,s3,N_modes3,dm3,Nbins_Na,Vbins_Na) 
    Nbins_Na=Nbins_Na/SUM(Nbins_Na); ! Normalize
	spm_Na=sum(Nbins_Na*v_p)*qNa ! Mean single particle mass (kg) 
	PN_Na=PM_Na_init/spm_Na! m^-3 
	
	
	! Specify the initial NaCl1 size distribution 
	N_modes3 = (/ 1D1, 0D0, 0D0 /)
    dm3 = (/ 1D-7, 1.5D-6, 7D-6 /)
    s3 = (/ 1.6, 1.6, 1.6 /)
	
	CALL dNdlogDp(d_g,dp_dry1,dlogDp,modes3,s3,N_modes3,dm3,Nbins_NaCl1,Vbins_NaCl1) 
    Nbins_NaCl1=Nbins_NaCl1/SUM(Nbins_NaCl1); ! Normalize
	spm_NaCl1=sum(Nbins_NaCl1*v_p)*qNa ! Mean single particle mass (kg) 
	PN_NaCl1=PM_NaCl1_init/spm_NaCl1! m^-3 
	
	
	! Specify the initial NaCl2 size distribution 
	N_modes3 = (/ 0D0, 1D1, 0D0 /)
    dm3 = (/ 1.5D-7, 1.5D-6, 7D-6 /)
    s3 = (/ 1.6, 1.6, 1.6 /)
	
	CALL dNdlogDp(d_g,dp_dry1,dlogDp,modes3,s3,N_modes3,dm3,Nbins_NaCl2,Vbins_NaCl2) 
    Nbins_NaCl2=Nbins_NaCl2/SUM(Nbins_NaCl2); ! Normalize
	spm_NaCl2=sum(Nbins_NaCl2*v_p)*qNa ! Mean single particle mass (kg) 
	PN_NaCl2=PM_NaCl2_init/spm_NaCl2! m^-3 

	! Specify the initial NaCl3 size distribution 
	N_modes3 = (/ 0D0, 0D0, 1D1 /)
    dm3 = (/ 1.5D-7, 1.5D-6, 7D-6 /)
    s3 = (/ 1.6, 1.6, 1.6 /)
	
	CALL dNdlogDp(d_g,dp_dry1,dlogDp,modes3,s3,N_modes3,dm3,Nbins_NaCl3,Vbins_NaCl3) 
    Nbins_NaCl3=Nbins_NaCl3/SUM(Nbins_NaCl3); ! Normalize
	spm_NaCl3=sum(Nbins_NaCl3*v_p)*qNa ! Mean single particle mass (kg) 
	PN_NaCl3=PM_NaCl3_init/spm_NaCl3! m^-3 

	DO j=1,Nz
	NbinsBC(j,:)=PN_BC(j)*Nbins_BC
	NbinsOC(j,:)=PN_OC(j)*Nbins_OC
	NbinsSO4(j,:)=PN_SO4(j)*Nbins_SO4
	NbinsNH4(j,:)=PN_NH4(j)*Nbins_NH4
	NbinsNa(j,:)=PN_Na(j)*Nbins_Na
	NbinsNaCl(j,:)=PN_NaCl1(j)*Nbins_NaCl1+PN_NaCl2(j)*Nbins_NaCl2+PN_NaCl3(j)*Nbins_NaCl3
	N_bins(j,:)=NbinsBC(j,:)+NbinsOC(j,:)+NbinsSO4(j,:)+NbinsNaCl(j,:)+NbinsNH4(j,:)+NbinsNa(j,:)
    V_bins(j,:)=N_bins(j,:)*dp_dry1**3*pi/6D0
	END DO

    c_p = 1D-100
    DO j=1,Nz
    c_p(j,1,:)=NbinsSO4(j,:)*(dp_dry1**3*pi/6D0)*qX(1)*1D-6/MX(1)*Na
	c_p(j,4,:)=NbinsNH4(j,:)*(dp_dry1**3*pi/6D0)*qX(4)*1D-6/MX(4)*Na
	c_p(j,5,:)=MX(5)/(MX(5)+MX(3))*NbinsNaCl(j,:)*(dp_dry1**3*pi/6D0)*qX(5)*1D-6/MX(5)*Na+&
	NbinsNa(j,:)*(dp_dry1**3*pi/6D0)*qX(5)*1D-6/MX(5)*Na
	c_p(j,3,:)=MX(3)/(MX(5)+MX(3))*NbinsNaCl(j,:)*(dp_dry1**3*pi/6D0)*qX(3)*1D-6/MX(3)*Na
	c_p(j,6,:)=NbinsBC(j,:)*(dp_dry1**3*pi/6D0)*qX(6)*1D-6/MX(6)*Na
	c_p(j,NSPEC_P,:)=NbinsOC(j,:)*(dp_dry1**3*pi/6D0)*qX(NSPEC_P)*1D-6/MX(NSPEC_P)*Na
	END DO
    c_p(:,7,:)=0.1*SUM(c_p,DIM=2)

! Start with clean air, see fig 8.15 (modeled dist) in Seinfeld and Pandis
!    dlogDp=log10(d(2:Nr_bins+1))-log10(d(1:nr_bins))
!    dm2 = (/ 1D-8, 1D-7 /) ! Mode diameter size distribution (m)
!    s2 = (/ 2D0, 2D0/) ! Standard deviation for each mode
!    N_modes2 = (/ 1D7, 1D8/) ! Number concentration modes (#/m^3)
!
!    dp_dry1=dp_dry(1,:)
!    CALL dNdlogDp(d_g,dp_dry1,dlogDp,modes2,s2,N_modes2,dm2,N_bins1,V_bins1) 
!    DO j=1,Nz
!      N_bins(j,:)=N_bins1
!      V_bins(j,:)=V_bins1
!    END DO

    ! Specify the characteristic ship primary particle emissions (Jonsson et al. 2011), table 1 (MM CA -cargo) is the largest mode. Assume that 0.99*0.9 of total particle nr is in that mode. 0.99*0.1 of the particle are in mode around 10 nm (see fig 1) and 0.01 of the particle are in mode around 100 nm.
    N_modes3 = (/ 3.6*0.99*0.1, 3.6*0.99*0.9, 3.6*0.01 /)
    dm3 = (/ 1D-8, 3.9D-8, 1D-7 /)
    s3 = (/ 1.59, 1.59, 1.5 /)
    CALL dNdlogDp(d_g, dp_dry1,dlogDp,modes3,s3,N_modes3,dm3,PN_emission_ship,PV_emission_ship) 
    PN_emission_ship=PN_emission_ship/SUM(PN_emission_ship); ! Normalize

    ! Specify the characteristic road traffic primary particle emissions (Kristensson et al.,2004) 
	! + sub 7 nm particle mode based on Hietikko et al., 2018  
    N_modes3 = (/ 35., 1700., 150. /)
    dm3 = (/ 3D-9, 2D-8, 77D-9 /)
    s3 = (/ 1.25, 1.9, 1.75 /)
 
    CALL dNdlogDp(d_g,dp_dry1,dlogDp,modes3,s3,N_modes3,dm3,PN_emission_S7,PV_emission_S7) 
    PN_emission_S7=PN_emission_S7/SUM(PN_emission_S7); ! Normalize
	
	! Specify the characteristic non-road machinery primary particle emissions  
    N_modes3 = (/ 1., 200., 50. /)
    dm3 = (/ 3D-9, 2D-8, 77D-9 /)
    s3 = (/ 1.25, 1.9, 1.75 /)
 
	CALL dNdlogDp(d_g,dp_dry1,dlogDp,modes3,s3,N_modes3,dm3,PN_emission_S9,PV_emission_S9) 
    PN_emission_S9=PN_emission_S9/SUM(PN_emission_S9); ! Normalize
	
  ! Specify the characteristic primary particle emissions from biomass burning   
    N_modes3 = (/ 0., 0., 1000./)
    dm3 = (/ 1D-8, 5D-8, 80D-9 /)
    s3 = (/ 1.6, 1.6, 1.6 /)
    CALL dNdlogDp(d_g, dp_dry1,dlogDp,modes3,s3,N_modes3,dm3,PN_biomass_burn,PV_biomass_burn) 
    PN_biomass_burn=PN_biomass_burn/SUM(PN_biomass_burn)


!! Primary particle number emission size distributions from GAINS model:
!dm=PN_emissions_dist_param_GAINS(1:5,1)
!where (dm>3D-7) dm=3D-7 ! Correct down the diameter of the accumulation model primary particle emissions 
!where (dm>1.5D-7) dm=dm/1.5 ! Correct down the diameter of the accumulation model primary particle emissions 
!s=PN_emissions_dist_param_GAINS(6:10,1)
!N_modes=PN_emissions_dist_param_GAINS(11:15,1)+1D-100
!where (dm<5D-9) N_modes=1D-100 ! Don't allow primary particle modes below 5 nm in diameter
!CALL dNdlogDp(d_g,dp_dry1,dlogDp,modes,s,N_modes,dm,PN_emission_agr,PV_emission_agr)
!PN_emission_agr=PN_emission_agr/SUM(PN_emission_agr) ! Normalize

!dm=PN_emissions_dist_param_GAINS(1:5,2)
!where (dm>3D-7) dm=3D-7 ! Correct down the diameter of the accumulation model primary particle emissions 
!where (dm>1.5D-7) dm=dm/1.5 ! Correct down the diameter of the accumulation model primary particle emissions 
!s=PN_emissions_dist_param_GAINS(6:10,2)
!N_modes=PN_emissions_dist_param_GAINS(11:15,2)+1D-100
!where (dm<5D-9) N_modes=1D-100 ! Don't allow primary particle modes below 5 nm in diameter
!CALL dNdlogDp(d_g,dp_dry1,dlogDp,modes,s,N_modes,dm,PN_emission_agr_burn,PV_emission_agr_burn)
!PN_emission_agr_burn=PN_emission_agr_burn/SUM(PN_emission_agr_burn) ! Normalize

!dm=PN_emissions_dist_param_GAINS(1:5,3)
!where (dm>3D-7) dm=3D-7 ! Correct down the diameter of the accumulation model primary particle emissions 
!where (dm>1.5D-7) dm=dm/1.5 ! Correct down the diameter of the accumulation model primary particle emissions 
!s=PN_emissions_dist_param_GAINS(6:10,3)
!N_modes=PN_emissions_dist_param_GAINS(11:15,3)+1D-100
!where (dm<5D-9) N_modes=1D-100 ! Don't allow primary particle modes below 5 nm in diameter
!CALL dNdlogDp(d_g,dp_dry1,dlogDp,modes,s,N_modes,dm,PN_emission_waste,PV_emission_waste)
!PN_emission_waste=PN_emission_waste/SUM(PN_emission_waste) ! Normalize

!dm=PN_emissions_dist_param_GAINS(1:5,4)
!where (dm>3D-7) dm=3D-7 ! Correct down the diameter of the accumulation model primary particle emissions 
!where (dm>1.5D-7) dm=dm/1.5 ! Correct down the diameter of the accumulation model primary particle emissions 
!s=PN_emissions_dist_param_GAINS(6:10,4)
!N_modes=PN_emissions_dist_param_GAINS(11:15,4)+1D-100
!where (dm<5D-9) N_modes=1D-100 ! Don't allow primary particle modes below 5 nm in diameter
!CALL dNdlogDp(d_g,dp_dry1,dlogDp,modes,s,N_modes,dm,PN_emission_traffic,PV_emission_traffic)
!PN_emission_traffic=PN_emission_traffic/SUM(PN_emission_traffic) ! Normalize

!dm=PN_emissions_dist_param_GAINS(1:5,5)
!where (dm>3D-7) dm=3D-7 ! Correct down the diameter of the accumulation model primary particle emissions 
!where (dm>1.5D-7) dm=dm/1.5 ! Correct down the diameter of the accumulation model primary particle emissions 
!s=PN_emissions_dist_param_GAINS(6:10,5)
!N_modes=PN_emissions_dist_param_GAINS(11:15,5)+1D-100
!where (dm<5D-9) N_modes=1D-100 ! Don't allow primary particle modes below 5 nm in diameter
!CALL dNdlogDp(d_g,dp_dry1,dlogDp,modes,s,N_modes,dm,PN_emission_power_plant,PV_emission_power_plant)
!PN_emission_power_plant=PN_emission_power_plant/SUM(PN_emission_power_plant) ! Normalize

!dm=PN_emissions_dist_param_GAINS(1:5,6)
!where (dm>3D-7) dm=3D-7 ! Correct down the diameter of the accumulation model primary particle emissions 
!where (dm>1.5D-7) dm=dm/1.5 ! Correct down the diameter of the accumulation model primary particle emissions 
!s=PN_emissions_dist_param_GAINS(6:10,6)
!N_modes=PN_emissions_dist_param_GAINS(11:15,6)+1D-100
!where (dm<5D-9) N_modes=1D-100 ! Don't allow primary particle modes below 5 nm in diameter
!CALL dNdlogDp(d_g,dp_dry1,dlogDp,modes,s,N_modes,dm,PN_emission_industry,PV_emission_industry)
!PN_emission_industry=PN_emission_industry/SUM(PN_emission_industry) ! Normalize

!dm=PN_emissions_dist_param_GAINS(1:5,7)
!where (dm>3D-7) dm=3D-7 ! Correct down the diameter of the accumulation model primary particle emissions 
!where (dm>1.5D-7) dm=dm/1.5 ! Correct down the diameter of the accumulation model primary particle emissions 
!s=PN_emissions_dist_param_GAINS(6:10,7)
!N_modes=PN_emissions_dist_param_GAINS(11:15,7)+1D-100
!where (dm<5D-9) N_modes=1D-100 ! Don't allow primary particle modes below 5 nm in diameter
!CALL dNdlogDp(d_g,dp_dry1,dlogDp,modes,s,N_modes,dm,PN_emission_flaring,PV_emission_flaring)
!PN_emission_flaring=PN_emission_flaring/SUM(PN_emission_flaring) ! Normalize

!dm=PN_emissions_dist_param_GAINS(1:5,8)
!where (dm>3D-7) dm=3D-7 ! Correct down the diameter of the accumulation model primary particle emissions 
!where (dm>1.5D-7) dm=dm/1.5 ! Correct down the diameter of the accumulation model primary particle emissions 
!s=PN_emissions_dist_param_GAINS(6:10,8)
!N_modes=PN_emissions_dist_param_GAINS(11:15,8)+1D-100
!where (dm<5D-9) N_modes=1D-100 ! Don't allow primary particle modes below 5 nm in diameter
!CALL dNdlogDp(d_g,dp_dry1,dlogDp,modes,s,N_modes,dm,PN_emission_wood,PV_emission_wood)
!PN_emission_wood=PN_emission_wood/SUM(PN_emission_wood) ! Normalize

!dm=PN_emissions_dist_param_GAINS(1:5,9)
!where (dm>3D-7) dm=3D-7 ! Correct down the diameter of the accumulation model primary particle emissions 
!where (dm>1.5D-7) dm=dm/1.5 ! Correct down the diameter of the accumulation model primary particle emissions 
!s=PN_emissions_dist_param_GAINS(6:10,9)
!N_modes=PN_emissions_dist_param_GAINS(11:15,9)+1D-100
!where (dm<5D-9) N_modes=1D-100 ! Don't allow primary particle modes below 5 nm in diameter
!CALL dNdlogDp(d_g,dp_dry1,dlogDp,modes,s,N_modes,dm,PN_emission_other_burn,PV_emission_other_burn)
!PN_emission_other_burn=PN_emission_other_burn/SUM(PN_emission_other_burn) ! Normalize

!dm=PN_emissions_dist_param_GAINS(1:5,10)
!where (dm>3D-7) dm=3D-7 ! Correct down the diameter of the accumulation model primary particle emissions 
!where (dm>1.5D-7) dm=dm/1.5 ! Correct down the diameter of the accumulation model primary particle emissions 
!s=PN_emissions_dist_param_GAINS(6:10,10)
!N_modes=PN_emissions_dist_param_GAINS(11:15,10)+1D-100
!where (dm<5D-9) N_modes=1D-100 ! Don't allow primary particle modes below 5 nm in diameter
!CALL dNdlogDp(d_g,dp_dry1,dlogDp,modes,s,N_modes,dm,PN_emission_coal,PV_emission_coal)
!PN_emission_coal=PN_emission_coal/SUM(PN_emission_coal) ! Normalize
	
 ! Actinic flux at the surface, calculated with the radiative transfer model
 OPEN(135, FILE='input/photolysis/Actinic_flux_1h_v2_'//date)
 READ(135,*) actinic_flux_1h
 CLOSE(135)
  
 

 
! Construct delta_lambda array (Seinfeld and Pandis table 4.2)
! Discrete wave lengths of solar irradiance
        lambda1 = (/ (i+0.5, i = 249,699) /)
        lambda2 = (/701, 711, 721, 731, 741, 751, 761, 771, 781, 791, 801, 821, 841, 861, 881, 901, 921, 941, 961, 981/)
        lambda3 = (/1002.5, 1052.5, 1102.5, 1152.5, 1202.5, 1252.5, 1302.5, 1352.5, 1402.5, 1452.5, 1502.5, &
            1552.5, 1602.5, 1652.5, 1702.5, 1752.5, 1802.5, 1852.5, 1902.5, 1952.5, 2002.5, 2107.5, 2212.5, &
            2302.5, 2402.5, 2517.5, 2617.5, 2702.5, 2832.5/)
        lambda4 = (/ 3025, 3235, 3425, 3665, 3855, 4085, 4575, 5085, 5925, 7785, 10075 /)
        lambda = (/ lambda1, lambda2, lambda3, lambda4 /)
        dlambda = lambda(2:size(lambda))-lambda(1:size(lambda)-1)

    ! KPP-chemistry setup 
    CALL KPP_SetUp ! This only called once, in the beginning

    M = 1.d-6*Ps(1:Nz,:)*Na/(Rg*Ts(1:Nz,:)) ! conc inert air molec [molec/cm^3]
	M_STP=1.d-6*1D5*Na/(Rg*273.15)
    O2 = 0.2095*M  ! in mole/cm3
    N2 = 0.7809*M
    H2O = 6.1078*exp(-(597.3-0.57*(Ts(1:Nz,:)-273.16))*18/1.986*(1/Ts(1:Nz,:)-1/273.16)) &
        *10/(1.38d-16*Ts(1:Nz,:))*RHs(1:Nz,:)
		

    ! Initial gas-phase species concentrations [molec/cm^3]
    conc = 1D-6
    conc(:,ind_CH4) = 1900d-9*M(:,1) ! CH4
    conc(:,ind_H2O2) = 0.01d-9*M(:,1) ! H2O2
    conc(:,ind_H2) = 0.5d-6*M(:,1) ! H2
 !   conc(:,ind_SO2) = 0.001d-9*M(:,1) ! SO2
    conc(:,ind_DMS) = 1d-15*M(:,1) ! SO2
    conc(:,ind_HNO3) = 1d-15*M(:,1) ! HNO3
	if (latitude(1)>0D0) then ! NH Arctic sim with high [O3] on Svalbard
    conc(:,ind_O3) = init_conc(1,:)*2D0; ! O3
	else
	conc(:,ind_O3) = init_conc(1,:); ! O3
	end if
	conc(:,ind_CO) = init_conc(2,:); ! CO 
	conc(:,ind_NO) = init_conc(3,:); ! NO 
    conc(:,ind_NO2) = init_conc(4,:); ! NO2
    conc(:,ind_SO2) = init_conc(5,:); ! SO2
	conc(:,ind_PAN) = init_conc(6,:); ! PAN
	conc(:,ind_HNO3) = init_conc(7,:); ! PAN
	conc(:,ind_H2O2) = init_conc(8,:); ! H2O2
	conc(:,ind_HCHO) = init_conc(9,:); ! HCHO
	conc(:,ind_C2H6) = init_conc(10,:); ! C2H6
	conc(:,ind_C3H8) = init_conc(11,:); ! C3H8
	conc(:,ind_C5H8) = init_conc(12,:); ! C5H8
	

	!conc(:,ind_CO) = 100d-9*M(:,1) ! CO
    conc(:,ind_HCl) = 1D0 ! HCl
    conc(:,ind_MSA) = 1D0 ! MSA
	conc(:,ind_HIO3) = 1D0 ! HIO3
	conc(:,ind_HIO2) = 1D0 ! HIO2
	!cHCl=1D0   ! Initial HCl concentration
    conc(:,ind_NH3) = 1D3   ! Initial NH3 concentration
    conc(:,ind_DMA) = 1D1   ! Initial DMA concentration
    
	!conc(:,ind_HOA) = 1D3        ! Initial HOA (POA) concentration in the gas-phase
	
    

    CALL fraction_POA_marine(dp_dry1,f_POA_marine)

    ! Concentration of each particle compounds in each layer and size bin (molecules/cm³)
 
    ! c_p = 1D-100
    ! !DO j=1,Nz
    ! !c_p(j,1,:)=MX(1)/(MX(1)+MX(4)+MX(3)+MX(5)+0.1*MX(NSPEC_P))*V_bins(j,:)*qX(1)*1D-6/MX(1)*Na 
	! !c_p(j,3,:)=MX(3)/(MX(1)+MX(4)+MX(3)+MX(5)+0.1*MX(NSPEC_P))*V_bins(j,:)*qX(3)*1D-6/MX(3)*Na 
	! !c_p(j,4,:)=MX(4)/(MX(1)+MX(4)+MX(3)+MX(5)+0.1*MX(NSPEC_P))*V_bins(j,:)*qX(4)*1D-6/MX(4)*Na
	! !c_p(j,5,:)=MX(5)/(MX(1)+MX(4)+MX(3)+MX(5)+0.1*MX(NSPEC_P))*V_bins(j,:)*qX(5)*1D-6/MX(5)*Na
	! c_p(j,1,:)=MX(1)/(MX(1)+MX(3)+MX(5)+0.1*MX(NSPEC_P))*V_bins(j,:)*qX(1)*1D-6/MX(1)*Na 
	! c_p(j,3,:)=MX(3)/(MX(1)+MX(3)+MX(5)+0.1*MX(NSPEC_P))*V_bins(j,:)*qX(3)*1D-6/MX(3)*Na 
	! !c_p(j,4,:)=MX(4)/(MX(1)+MX(4)+MX(3)+MX(5)+0.1*MX(NSPEC_P))*V_bins(j,:)*qX(4)*1D-6/MX(4)*Na
	! c_p(j,5,:)=MX(5)/(MX(1)+MX(3)+MX(5)+0.1*MX(NSPEC_P))*V_bins(j,:)*qX(5)*1D-6/MX(5)*Na
	
	! !c_p(j,7,:)=MX(7)/(MX(1)+MX(4)+MX(3)+MX(5)+MX(7)+0.1*MX(NSPEC_P))*V_bins(j,:)*qX(7)*1D-6/MX(7)*Na
    ! !c_p(j,NSPEC_P,:)=0.1*MX(NSPEC_P)/(MX(1)+MX(4)+MX(3)+MX(5)+0.1*MX(NSPEC_P))*V_bins(j,:)*qX(NSPEC_P)*1D-6/MX(NSPEC_P)*Na
    ! c_p(j,NSPEC_P,:)=0.1*MX(NSPEC_P)/(MX(1)+MX(4)+MX(3)+MX(5)+0.1*MX(NSPEC_P))*V_bins(j,:)*qX(NSPEC_P)*1D-6/MX(NSPEC_P)*Na
    
	! END DO
    ! c_p(:,7,:)=0.1*SUM(c_p,DIM=2)
	

    c_p_nucl=0.
	! c_p_nuclDMA=0.
    !! IF NH3 + H2SO4 clusters are formed
    c_p_nucl(1)=MX(1)/(MX(1)+MX(4))*V_bins(1,1)/N_bins(1,1)*qX(1)/MX(1)*Na
    c_p_nucl(4)=MX(4)/(MX(1)+MX(4))*V_bins(1,1)/N_bins(1,1)*qX(4)/MX(4)*Na
	!c_p_nucl(7)=MX(7)/(MX(1)+MX(4)+MX(7))*V_bins(1,1)/N_bins(1,1)*qX(7)/MX(7)*Na ! Estimate initial particle water content
    c_p_nucl(7)=0.1*SUM(c_p_nucl)
	
	c_p_nuclDMA=0.
	! IF DMA + H2SO4 clusters are formed
    c_p_nuclDMA(1)=MX(1)/(MX(1)+MX(11))*V_bins(1,2)/N_bins(1,2)*qX(1)/MX(1)*Na
    c_p_nuclDMA(11)=MX(11)/(MX(1)+MX(11))*V_bins(1,2)/N_bins(1,2)*qX(11)/MX(11)*Na
	c_p_nuclDMA(7)=0.1*SUM(c_p_nuclDMA)
	
    ! IF the particle number concentration in any size bin go below a minimum value the concentration
	! c_p_nuclIO=0.
	! ! IF HIO3 + HIO2 clusters are formed
    ! c_p_nuclIO(10)=MX(10)/(MX(10)+MX(12))*V_bins(1,2)/N_bins(1,2)*qX(1)/MX(1)*Na
    ! c_p_nuclIO(12)=MX(11)/(MX(10)+MX(12))*V_bins(1,2)/N_bins(1,2)*qX(11)/MX(11)*Na
	! c_p_nucIO(7)=0.1*SUM(c_p_nuclIO)
	
    ! IF the particle number concentration in any size bin go below a minimum value the concentration
    ! is corrected up to the particle composition of c_p_backg.
    c_p_backg=0D0
	c_p_backg(1,:)=MX(1)/(MX(1)+MX(4)+MX(3)+MX(5)+0.1*MX(NSPEC_P))*V_bins(1,:)/N_bins(1,:)*qX(1)*1D-6/MX(1)*Na 
	c_p_backg(3,:)=MX(3)/(MX(1)+MX(4)+MX(3)+MX(5)+0.1*MX(NSPEC_P))*V_bins(1,:)/N_bins(1,:)*qX(3)*1D-6/MX(3)*Na 
	c_p_backg(4,:)=MX(4)/(MX(1)+MX(4)+MX(3)+MX(5)+0.1*MX(NSPEC_P))*V_bins(1,:)/N_bins(1,:)*qX(4)*1D-6/MX(4)*Na
	c_p_backg(5,:)=MX(5)/(MX(1)+MX(4)+MX(3)+MX(5)+0.1*MX(NSPEC_P))*V_bins(1,:)/N_bins(1,:)*qX(5)*1D-6/MX(5)*Na 
	!c_p_backg(7,:)=MX(7)/(MX(1)+MX(4)+MX(3)+MX(5)+MX(7)+0.1*MX(NSPEC_P))*V_bins(1,:)/N_bins(1,:)*qX(7)*1D-6/MX(7)*Na 
	c_p_backg(NSPEC_P,:)=0.1*MX(NSPEC_P)/(MX(1)+0.1*MX(4)+MX(3)+MX(5)+0.1*MX(NSPEC_P))*V_bins(1,:)/N_bins(1,:)*qX(NSPEC_P)*1D-6/MX(NSPEC_P)*Na
    c_p_backg(7,:)=0.1*SUM(c_p_backg,DIM=1)

    DO i=1,Nz
    DO j = 1,nr_bins
        vp_dry(i,j) = SUM(c_p(i,index_dry,j)/Na*MX(index_dry)/qX(index_dry)/(N_bins(i,j)*1D-6)) ! m^3
        vp_wet(i,j) = SUM(c_p(i,:,j)/Na*MX/qX/(N_bins(i,j)*1D-6)) ! m^3
    END DO
    END DO
    dp_dry = (vp_dry*6D0/pi)**(1D0/3D0) ! Dry particle diameters
    d_p = (vp_wet*6D0/pi)**(1D0/3D0) ! Wet particle diameters

    ! Initialise the aerosol particle water content:
    DO ii=1,Nz
      N_bins1=N_bins(ii,:);pH1=pH(ii,:);Kprim_HNO31=Kprim_HNO3(ii,:)
      Kprim_HCl1=Kprim_HCl(ii,:);Hprim_NH31=Hprim_NH3(ii,:);Kprim_NH31=Kprim_NH3(ii,:)
      fHSO41=fHSO4(ii,:);fSO41=fSO4(ii,:);fNO31=fNO3(ii,:);fCl1=fCl(ii,:);mHCO31=mHCO3(ii,:)
      mCO31=mCO3(ii,:);mOH1=mOH(ii,:);mCOO1=mCOO(ii,:)
      W1=W(ii,:);y1=y(ii,:,:)
	  Kprim_CH3SO3H1=Kprim_CH3SO3H(ii,:)
	  Kprim_HIO31=Kprim_HIO3(ii,:)
	  fCH3SO31=fCH3SO3(ii,:)
	  fHIO31=fHIO3(ii,:)
	  c_p1=c_p(ii,:,:)
      DO i = 1,2
        surf_tens = (76.1-0.155*(Ts(ii,tr)-273.15))*1D-3 ! Surface tension of pure water dyn m^-1
        S_Kelvin = exp(2D0*surf_tens*MX(7)/(d_p(ii,:)/2D0*Rg*Ts(ii,tr)*qX(7))) 
        aw = RHs(ii,tr)/1D2/S_Kelvin  ! Water activity
       
	    CALL thermodyn_AIOMFAC_inorg(Ts(ii,tr),c_p1,N_bins1,y1,&
        conc(ii,ind_NH3),conc(ii,ind_HNO3),conc(ii,ind_HCl),aw,pCO2,pH1,Kprim_HNO31,&
		Kprim_HCl1,Kprim_CH3SO3H1,Kprim_HIO31,Hprim_NH31,Kprim_NH31,&
        fHSO41,fSO41,fNO31,fCl1,fCH3SO31,fHIO31,mHCO31,mCO31,mOH1,W1)
		    
        DO j = 1,nr_bins
          vp_dry(ii,j) = SUM(c_p(ii,index_dry,j)/Na*MX(index_dry)/qX(index_dry)/(N_bins(ii,j)*1D-6)) ! m^3
          vp_wet(ii,j) = SUM(c_p(ii,:,j)/Na*MX/qX/(N_bins(ii,j)*1D-6)) ! m^3
          V_bins(ii,j) = N_bins(ii,j)*vp_wet(ii,j)
          dens_p(ii,j) = SUM(c_p(ii,:,j)*MX*1D6)/Na/V_bins(ii,j) ! Total particle density
        END DO
        dp_dry(ii,:) = (vp_dry(ii,:)*6D0/pi)**(1D0/3D0) ! Dry particle diameters
        d_p(ii,:) = (vp_wet(ii,:)*6D0/pi)**(1D0/3D0) ! Dry particle diameters
      END DO
      N_bins(ii,:)=N_bins1; c_p(ii,:,:)=c_p1; pH(ii,:)=pH1; Kprim_HNO3(ii,:)=Kprim_HNO31
      Kprim_HCl(ii,:)=Kprim_HCl1; Hprim_NH3(ii,:)=Hprim_NH31; Kprim_NH3(ii,:)=Kprim_NH31
      fHSO4(ii,:)=fHSO41; fSO4(ii,:)=fSO41; fNO3(ii,:)=fNO31; fCl(ii,:)=fCl1; mHCO3(ii,:)=mHCO31
      mCO3(ii,:)=mCO31; mOH(ii,:)=mOH1; mCOO(ii,:)=mCOO1
      W(ii,:)=W1;y(ii,:,:)=y1
	  Kprim_CH3SO3H(ii,:)=Kprim_CH3SO3H1
	  Kprim_HIO3(ii,:)=Kprim_HIO31
	  fCH3SO3(ii,:)=fCH3SO31
	  fHIO3(ii,:)=fHIO31
	
     cw=SUM(c_p(ii,7,:))/Na*MH2O ! Inorganic bulk liquid water content	kg/cm^3 = L/cm^3	
    END DO


 DO j=1,Nz
 ! Initial particle aqueous phase chemical composition:
   conc(j,ind_Clion) = SUM(c_p(j,3,:)) ! Total Cl- particle conc
   conc(j,ind_Brion) = conc(j,ind_Clion)*1D-4 !0.0015 ! Typical sea water Br/Cl ratio von Glasow and Crutsen, 2004 and Millero et al., 2008 0.0015
   conc(j,ind_NO3ion) = SUM(c_p(j,2,:)) ! Total NO3- particle conc
   conc(j,ind_NH4ion) = SUM(c_p(j,4,:)) ! Total NH4+ particle conc
   conc(j,ind_SO42ion) = SUM(fSO4(j,:)*c_p(j,1,:)) ! Total SO42- particle conc
   conc(j,ind_HSO4ion) = SUM(fHSO4(j,:)*c_p(j,1,:)) ! Total HSO4- particle conc
   conc(j,ind_CH3SO3ion) = SUM(c_p(j,9,:)) ! Total MSA particle conc
   conc(j,ind_CH3SO3Haq) = 0D0
 END DO
 

    ! Parameters needed to estimate a-pinene emission as done in Tunved study, to be replaced by LPJ-GUESS 
    E_pot=1.5  ! Monterpene emission potential in ug/(g(dry leaf)*h) for coniferous trees
    Bio_dens=950   ! average foliar biomass density (in g/m^2) of pine and spruce in northern parts of the boreal zone
    T_s=303.15    ! Leaf temperature at standard conditions (K)
    beta=0.09    ! Empirical coefficient (1/C) 

    !r_elvoc_index=(/ind_R_ELVOC_O5,ind_R_ELVOC_O6,ind_R_ELVOC_O7,ind_R_ELVOC_O8,ind_R_ELVOC_O9,ind_R_ELVOC_O10,&
    !ind_R_ELVOC_O11,ind_R_ELVOC_O12,ind_R_ELVOC_O13,ind_R_HOM_O5OH,,ind_R_HOM_O6OH,ind_R_HOM_O7OH,ind_R_HOM_O8OH,&
    !ind_R_HOM_O9OH,ind_R_HOM_O10OH,ind_R_HOM_O11OH,ind_R_HOM_O12OH,ind_R_HOM_O13OH/)

    !---------------------------------------------------------------------------!
    !---------------------------------------------------------------------------!  
    ! Start model                                                                !
    !---------------------------------------------------------------------------!
    !---------------------------------------------------------------------------!
    
    time = t_start

      
    loop: DO WHILE (time <= t_end)
    ! write(*,*) 'Start here time', Time 
	conc(:,ind_DUMMY)=0D0
	
	x1=REAL(tr,8)-time/60D0
	IF (tr<traj_len_min) THEN
	RH=RHs(:,tr)*x1+RHs(:,tr+1)*(1D0-x1)
	ELSE
	RH=RHs(:,tr)
	END IF
	
	WHERE ((RH_old-RH)>2D0 .AND. RH_old<90D0) RH=RH_old-2D0 ! Limit the RH change per time step 
	WHERE ((RH_old-RH)<-2D0  .AND. RH_old<90D0) RH=RH_old+2D0 ! Limit the RH change per time step
	WHERE ((RH_old-RH)>1D0 .AND. RH_old<98D0) RH=RH_old-1D0 ! Limit the RH change per time step 
	WHERE ((RH_old-RH)<-1D0  .AND. RH_old<98D0) RH=RH_old+1D0 ! Limit the RH change per time step
	WHERE ((RH_old-RH)>0.3D0 .AND. RH_old>=98D0) RH=RH_old-0.3D0 ! Limit the RH change per time step 
	WHERE ((RH_old-RH)<-0.3D0  .AND. RH_old>=98D0) RH=RH_old+0.3D0 ! Limit the RH change per time step
	WHERE ((RH_old-RH)>0.1D0 .AND. RH_old>=99.8D0) RH=RH_old-0.1D0 ! Limit the RH change per time step 
	WHERE ((RH_old-RH)<-0.1D0  .AND. RH_old>=99.8D0) RH=RH_old+0.1D0 ! Limit the RH change per time step
	WHERE (RH>99.9D0) RH=99.9D0
	
  ! Don't allow the particle number concentration in any size bin drop below 0.1 m^-3:

   DO j = 1,nz
   DO i=1,nr_bins
           IF (N_bins(j,i)<1D-1) THEN
           N_bins(j,i)=1D0
           c_p(j,:,i)=c_p_backg(:,i)*1D0;
           vp_dry(j,i) = SUM(c_p(j,index_dry,i)/Na*MX(index_dry)/qX(index_dry)/(N_bins(j,i)*1D-6)) ! m^3
           vp_wet(j,i) = SUM(c_p(j,:,i)/Na*MX/qX/(N_bins(j,i)*1D-6)) ! m^3
           V_bins(j,i) = N_bins(j,i)*vp_wet(j,i)
           dens_p(j,i) = SUM(c_p(j,:,i)*MX*1D6)/Na/V_bins(j,i) ! Total particle density
           dp_dry(j,i) = (vp_dry(j,i)*6D0/pi)**(1D0/3D0) ! Dry particle diameter
           d_p(j,i) = (vp_wet(j,i)*6D0/pi)**(1D0/3D0) ! Wet particle diameter
           END IF
       END DO
    END DO


    IF (PBLH(tr) <= 100.) THEN
     PBLH(tr) = 100.
    END IF
    IF (PBLH(tr) >= 2000.) THEN 
     PBLH(tr) = 2000.
    END IF

    !---------------------------------------------------------!
    ! Mixing (turbulent diffusion) in the vertical direction  !
    !---------------------------------------------------------!
    IF (diff_coeff_index == 1) THEN ! Calculate diffusivity coefficient
       CALL diffusivity_coeff(fric_veloc(tr),T2m(tr),RH2m(tr),P2m(tr),SHTFs(tr),Ps(:,tr),Ts(:,tr),cloudLWC(:,tr),PBLH(tr),RH,u_vel(:,tr), v_vel(:,tr),Kz,tr)     
   !    CALL diffusivity_coeff(u_vel(:,tr),v_vel(:,tr),Ts(:,tr),Ps(:,tr),RHs(:,tr),RH2m(tr),T2m(tr),Kz)
    ELSE ! Assume fix diffusivity coefficient 
       Kz=5D0 ! m^2/s (assumed to be fixed to start with)
    END IF   

    DO jj=1,NSPEC
        conc_diff=conc(:,jj)    
    IF (jj==ind_O3) THEN
        CALL diffusion1D_variable_z(Kz,conc_diff,k_uper_O3,dt)
    ELSE
        CALL diffusion1D_variable_z(Kz,conc_diff,k_uper_BC,dt)
    END IF
        conc(:,jj)=conc_diff
    END DO
        
    !CALL diffusion1D_variable_z(Kz,cNH3,k_uper_BC,dt)
	CALL diffusion1D_variable_z(Kz,c222Rn,k_uper_BC,dt)
    !CALL diffusion1D_variable_z(Kz,cHCl,k_uper_BC)

! Problem if k_uper_BC < 0
k_uper_BC=0D0 
   ! Update all particle concentrations:
      DO i=1,nr_bins
         conc_diff=N_bins(:,i)    
         CALL diffusion1D_variable_z(Kz,conc_diff,k_uper_BC,dt)
         N_bins(:,i)=conc_diff
         DO jj=1,NSPEC_P
            conc_diff=c_p(:,jj,i)
            CALL diffusion1D_variable_z(Kz,conc_diff,k_uper_BC,dt)
            c_p(:,jj,i)=conc_diff
         END DO 
      END DO
	  
! Mixing of ACDC clusters, added by Pontus 2021-05-26

  DO i=1,chem_1%neq_syst !63
         conc_diff=c_clusters_n(:,i)    
         CALL diffusion1D_variable_z(Kz,conc_diff,k_uper_BC,dt)
         c_clusters_n(:,i)=conc_diff
  END DO
  DO i=1,chem_2%neq_syst !31
         conc_diff=c_clusters_d(:,i)    
         CALL diffusion1D_variable_z(Kz,conc_diff,k_uper_BC,dt)
         c_clusters_d(:,i)=conc_diff
  END DO
  DO i=1,chem_3%neq_syst !31
         conc_diff=c_clusters_i(:,i)    
         CALL diffusion1D_variable_z(Kz,conc_diff,k_uper_BC,dt)
         c_clusters_i(:,i)=conc_diff
  END DO
	
    DO ii=1,Nz
      DO j = 1,nr_bins
        vp_dry(ii,j) = SUM(c_p(ii,index_dry,j)/Na*MX(index_dry)/qX(index_dry)/(N_bins(ii,j)*1D-6)) ! m^3
        vp_wet(ii,j) = SUM(c_p(ii,:,j)/Na*MX/qX/(N_bins(ii,j)*1D-6)) ! m^3
        V_bins(ii,j) = N_bins(ii,j)*vp_wet(ii,j)
        dens_p(ii,j) = SUM(c_p(ii,:,j)*MX*1D6)/Na/V_bins(ii,j) ! Total particle density
      END DO
        dp_dry(ii,:) = (vp_dry(ii,:)*6D0/pi)**(1D0/3D0) ! Dry particle diameters
        d_p(ii,:) = (vp_wet(ii,:)*6D0/pi)**(1D0/3D0) ! Dry particle diameters
    END DO  

    !--------------------------------------------------------!
    ! Gas-phase emission                                     !
    !--------------------------------------------------------!
    E_gases = 0.

 !************** Call MEGAN model and calculate BVOC emissions or use LPJ-Guess emissions: *******************************************************************************

IF (index_MEGAN==1) THEN
EM_PAR(tr)=DSWF(tr) * 4.766 * 0.6! Estimated incoming photosynthetic active radiation [umol/m2/s] for MEGAN

EM_Lat=input_MEGAN(1,tr) 
EM_Long=input_MEGAN(2,tr)
time_day = input_MEGAN(3,tr) ! Local time of the day in hours (0-24)
JD=INT(input_MEGAN(4,tr));    ! Julian Day

! For SMEAR II MT standard emissions is 984 (ug/m^2/h) according to MEGAN
! But according to Tarvainen et al (2005) the early summer 20th May to 11th June 
! emissions at SMEAR II is 5184 ng g−1 h−1 which gives 2789 (ug/m^2/h) 2.8343 
! times larger than MEGAN.
! Before sim nr >24 3.0*input_MEGAN(6:27,tr)

EM_EFin =6.0*input_MEGAN(6:27,tr) ! Standard emission factors from http://lar.wsu.edu/megan/guides.html (ug/m^2/h)

EM_LAI_year = input_MEGAN(5,tr) ! (Currently use LAI values from LPJ-GUESS)

EM_SRAD = DSWF(tr)                ! Incoming short wave solar radiation in (W/m^2)
IF (tr>24*60) THEN
    EM_TempK_day    = SUM(Ts(1,tr-60*24:tr))/(60.*24.+1.)             ! Daily average temperature Kelvin
    EM_SRAD_day     = SUM(DSWF(tr-60*24:tr))/(60.*24.+1.)             ! Daily average short wave radiation (W/m2)
    EM_PAR_day      = SUM(EM_PAR(tr-60*24:tr))/(60.*24.+1.)            ! Incoming PAR daily average          
  ELSE
    EM_TempK_day    = SUM(Ts(1,1:tr))/(REAL(tr))                      ! Daily average temperature Kelvin
    EM_SRAD_day     = SUM(DSWF(1:tr))/(REAL(tr))                      ! Daily average short wave radiation (W/m2)
    EM_PAR_day      = SUM(EM_PAR(1:tr))/(REAL(tr))                    ! Incoming PAR daily average          
END IF


IF (landuse_index(tr)==1) THEN
EM_Cantype = 1 ! Needle leaf trees
IF (EM_EFin(1)>100.) THEN
EM_EFin(1) = 100. ! For boreal forest of Europe the emissions of Isoprene is much smaller than in US
END IF
ELSEIF (landuse_index(tr)==2) THEN
EM_Cantype = 2 ! Broad leaf trees 
ELSEIF (landuse_index(tr)==3) THEN  ! grass
EM_Cantype = 3 ! grass
ELSEIF (landuse_index(tr)==4) THEN ! desert (LAI = 0)
EM_Cantype = 4 ! shrubs
ELSEIF (landuse_index(tr)==5) THEN ! shrubs
EM_Cantype = 4  ! shrubs
ELSEIF (landuse_index(tr)==6) THEN ! sea (LAI = 0)
EM_Cantype = 4  ! shrubs
ELSEIF (landuse_index(tr)==7) THEN ! Urban areas
EM_Cantype = 4  ! shrubs
ELSEIF (landuse_index(tr)==8) THEN ! Ice and snow
EM_Cantype = 4  ! shrubs
END IF  

          ! First put all emissions to zero at each time step
          EM_EMI = 0D0
         
          ! Input
          EM_Julian       = JD                             ! Julian day
          EM_DATE         = ((1000*EM_Year)+EM_Julian)     ! in format YYYDDD scalar
          EM_Time_M2_R    = time_day*3600.                 ! Time (s)
          EM_Time_M2      = Timeformatter(EM_Time_M2_R)    ! Function in megan module: integer output - real input
          
          EM_kz           = kz_can                         ! Number of vertical layers (not only canopy)
          EM_Can_Lay      = Can_Lay                        ! Number of vertical layers inside the canopy
          
          EM_z            = z(2:EM_kz+1)                  ! Array of height of the layers (m)            
          EM_TempK        = Ts(1:EM_kz,tr)                ! Temperature in Kelvin
          EM_TempC        = Ts(1:EM_kz,tr) - 273.15       ! Temperature in Celsius
          EM_RH           = RH(1:EM_kz)               ! Relative humidity in %
          EM_dz(1)=EM_z(1)
          EM_dz(2:EM_kz)= EM_z(2:EM_kz)-EM_z(1:EM_kz-1)
         
          ! Calculate in-canopy wind speed profile:
          
          CALL in_canopy_wind_profile(Ts(1,tr),Ps(1,tr),fric_veloc(tr),landuse_index(tr),EM_Wind,EM_z)

          EM_LAD=0D0
          EM_LAD(1)       = 0.15!0.01                     ! Leaf area density as [0,1] in the canopy layer 1: 0-3 m
          EM_LAD(2)       = 0.45!0.13                     ! Leaf area density as [0,1] in the canopy layer 2: 3-9 m 
          EM_LAD(3)       = 0.40                          ! Leaf area density as [0,1] in the canopy layer 3: 9-18 m
          EM_LAD(4)       = 0.0                           ! Leaf area density as [0,1] in the canopy layer 4: 18-30 m
          EM_LAD(5)       = 0.0                           ! Leaf area density as [0,1] in the canopy layer 5: 30-45 m

          EM_Pres         = Ps(1,tr)                       ! Pressure in the canopy (Pa)
          EM_SMOIST       = 0.3                            ! (IF > 0.2 no effect on emissions) Soil moisture (%)
  
          ! LAI_Megan is for Hyytiälä monthly distribution
          EM_LAI = LAI_Megan(Mon) !EM_LAI_year * LAI_Megan(Mon)/6.5728     ! One sided LAI (for needle leaf trees it is the projected area of the 3D needle)
          
          IF (Mon > 1) THEN
             EM_LAI_past     = LAI_Megan(Mon-1) !EM_LAI_year * LAI_Megan(Mon-1)/6.5728
          ELSE
             EM_LAI_past     = LAI_Megan(12) !EM_LAI_year * LAI_Megan(12)/6.5728
          END IF

          ! Output values set to zero for each run
          EM_SUN_PAR      = 0.                             ! Array of sun fraction - 1 above the canopy and decreases inside the canopy
          EM_SunleafTK    = 0.                             ! Array of temparture for sun leaf
          EM_ShadeleafTK  = 0.                             ! Array of temparture for shade leaf
          EM_Sunfrac      = 0.                             ! Array of the fraction of sun leaves. i = 1 is the top canopy layer, 2 is the next layer, etc.
          EM_EMI          = 0.                             ! Emissions in molecules per cm3 per second
          EM_Ea1pL_M2     = 0.                             ! Emission factor light
          EM_Ea1tL_M2     = 0.                             ! Emission factor temperature
          EM_Ea1NL_M2     = 0.                             ! Emission factor compined
          EM_GAM_TMP      = 0.                             ! Gamma factor for temperature, Non isoprene
          EM_GAM_OTHER    = 0.                             ! Gamma factors for other parameters

          !  Layer list NOTE megan used the canopy for added accuracy so "cover" appears in every part of the canopy, cause it just changes a few variables
          ! in the calculations in comparison to the other plant functional types, (g/s)
          EM_ER_Grass     = 0.0                            ! emission rate from grass
          EM_ER_SB        = 0.0                            ! emission rate from Shrubland(s) , see note above
          EM_ER_NT        = 0.0                            ! emission rate from Needle Trees, - || -
          EM_ER_BT        = 0.0                            ! emission rate from Broadleaf Trees, - || -


          DO jj = 1,EM_kz
             EM_ES(jj)    = (a0 + a1 * EM_TempC(jj)**1 + a2 * EM_TempC(jj)**2 + a3 * EM_TempC(jj)**3                 &
                                + a4 * EM_TempC(jj)**4 + a5 * EM_TempC(jj)**5 + a6 * EM_TempC(jj)**6)* 100

             ! Water vapour pressure in Pa
             EM_EW(JJ)  = EM_RH(jj) * EM_ES(jj) / 100

             ! Water vapour mixing ratio (g/kg) (based on Megan conversion function)
             EM_WVM(jj) = 18.016 / 28.97 * EM_EW(jj) / (EM_Pres - EM_EW(jj))*1000

             ! Limit EM_EW to maximum value of ES
             IF (EM_EW(jj) .GT. EM_ES(jj)) THEN
                EM_EW(jj) = EM_ES(jj)
             ENDIF
          ENDDO

    Call EMISSION_M2(EM_Julian, EM_Can_Lay, EM_Lat, EM_Long, EM_Date, EM_Time_M2, EM_PAR(tr), EM_PAR_day,                                       &
                              EM_TempK, EM_TempK_day, EM_Pres, EM_WVM, EM_WIND, EM_SMOIST, EM_LAD, EM_LAI, EM_LAI_past, EM_Cantype,                  &

                              EM_ER, EM_VAR, EM_ER_Grass, EM_ER_SB, EM_ER_NT, EM_ER_BT, EM_Ea1pL_M2, EM_Ea1NL_M2, EM_Ea1tL_M2,                                   &
                              EM_GAM_TMP, EM_GAM_OTHER, EM_z, EM_kz, EM_EMI, EM_sun_par, EM_SunleafTK, EM_ShadeleafTK, EM_Sunfrac, EM_EFin)

! Sum the in-canopy emissions from MEGAN and distribute them over the entire ADCHEM surface layer (surface grid cell) width (dz). 
EMI=EM_EMI*3.6E3 ! molec cm^-3 h^-1


!EMI=0.
!DO j=1,EM_kz
!    EMI=EMI+EM_dz(j)*EM_EMI(j,:)/dz*3600. ! molec cm^-3 h^-1
!END DO 

! 1 =  Isoprene, 2 =  MBO (2methyl-3buten-2ol),  3 =  MYRC, 4 =  Sabinene
! 5 =  Limonen, 6 =  3-Carene, 7 =  Ocimene, 8 =  Beta-pinene, 9 =  Alpha-pinene
! 10 = FARN, 11 = Betacarophylene, 12 = Methanol, 13 = Aceton, 14 = Acetaldehyde
! 15 = Formaldehyde, 16 = Methan, 17 = NO, 18 = Other monoterpene, 19 = Other sesquiterpenes
! 20 = CO, 21 = Cineole, 22 = Linalool

! According to measurements in Hyytiälä the fraction of a-pinene and carene is dominating the emissions. 
! Thus, we will assume that all extra BVOCs from MEGAN not included in MCM will go to a-pinenen
! that represent carnene and a-pinene species. The limonene fraction at Hyytiälä is very low ~2%
! beta-pinene is about 9 % and more than 80 % of the MT are a-pinene and carene.

 E_gases(:,ind_APINENE) = EMI(:,9)  ! Emission of apinene from MEGAN (molecules cm⁻³ h⁻¹)
 E_gases(:,ind_CARENE)=EMI(:,3)+EMI(:,4)+EMI(:,6)+EMI(:,7)+EMI(:,18) ! All other MT emissions from MEGAN assumed to be Carene
 E_gases(:,ind_BPINENE) = EMI(:,8) ! Emission of bpinene from MEGAN (molecules cm⁻³ h⁻¹)
 E_gases(:,ind_LIMONENE) = EMI(:,5) ! Emission of limonene from MEGAN (molecules cm⁻³ h⁻¹)
 E_gases(:,ind_C5H8) = EMI(:,1) ! Emission of isoprene from MEGAN (molec/cm³/h)
 E_gases(:,ind_BCARY) = (SUM(EMI(:,2:9),DIM=2)+EMI(:,18))*0.05D0;
 
! ! Estimated emissions of small BVOCs based on CAMS average global emissions scaled with the a-pinene emissions from MEGAN: 
E_gases(:,ind_CH3OH) = EMI(:,9)*16.4874
E_gases(:,ind_CH3COCH3) = EMI(:,9)*2.9652
E_gases(:,ind_CH3CHO) = EMI(:,9)*1.6357
E_gases(:,ind_HCHO) = EMI(:,9)*0.5997
E_gases(:,ind_C3H6) = EMI(:,9)*1.7198
E_gases(:,ind_C2H4) = EMI(:,9)*4.1579
E_gases(:,ind_C2H5OH) = EMI(:,9)*1.5645
E_gases(:,ind_HCOOH) = EMI(:,9)*0.2876
E_gases(:,ind_CH3CO2H) = EMI(:,9)*0.2205

       
    ELSE

!      E_gases(1,ind_APINENE) = EmissionAPIN(tr)/(dz(1)*1D2) ! Emission of apinene from LPJ_GUESS (molecules cm⁻³ h⁻¹)
!      E_gases(1,ind_BPINENE) = EmissionBPIN(tr)/(dz(1)*1D2) ! Emission of bpinene from LPJ_GUESS (molecules cm⁻³ h⁻¹)
!      E_gases(1,ind_LIMONENE) = EmissionLIM(tr)/(dz(1)*1D2) ! Emission of apinene from LPJ_GUESS (molecules cm⁻³ h⁻¹)
!      E_gases(1,ind_C5H8) = EmissionISO(tr)/(dz(1)*1D2) ! Emission of isoprene from LPJ-GUESS (molec/cm³/h)
!      E_gases(1,ind_CARENE) = EmissionOTH_MT(tr)/(dz(1)*1D2) ! Emission of other mon assumed to be carene
      
 END IF

! *******************************************************************************************************************************
    ! EMEP emissions (molecules cm^-3 h^-1)
	E_gases(1,ind_CO) = (sum(EmissionCO(4:6,tr))+sum(EmissionCO(8:12,tr)))/(dz(1)*100.)
	E_gases(4,ind_CO) = EmissionCO(7,tr)/(dz(4)*100.) ! Ships
	E_gases(6,ind_CO) = sum(EmissionCO(1:3,tr))/(dz(6)*100.) ! Power Plants and Industry
	
    E_gases(1,ind_NO2) = (0.6*sum(EmissionNO2(4:6,tr))+0.6*sum(EmissionNO2(8:12,tr))+EmissionNOxsoil(tr))/(dz(1)*100.)
	E_gases(4,ind_NO2) = 0.6*EmissionNO2(7,tr)/(dz(4)*100.) ! Ships
	E_gases(6,ind_NO2) = 0.6*sum(EmissionNO2(1:3,tr))/(dz(6)*100.) ! Power Plants and industries emit  at higher altitudes
	
	E_gases(1,ind_NH3) = (sum(EmissionNH3(4:6,tr))+sum(EmissionNH3(8:12,tr))+EmissionNH3_birds(tr))/(dz(1)*100.)
	E_gases(4,ind_NH3) = EmissionNH3(7,tr)/(dz(4)*100.) ! Ships
	E_gases(6,ind_NH3) = sum(EmissionNH3(1:3,tr))/(dz(6)*100.) ! Power Plants and Industry
	
	E_gases(1,ind_DMA) = f_DMA_NH3*(sum(EmissionNH3(4:6,tr))+sum(EmissionNH3(8:12,tr))+EmissionNH3_birds(tr))/(dz(1)*100.)
	E_gases(4,ind_DMA) = f_DMA_NH3*EmissionNH3(7,tr)/(dz(4)*100.) ! Ships
	E_gases(6,ind_DMA) = f_DMA_NH3*sum(EmissionNH3(1:3,tr))/(dz(6)*100.) ! Power Plants and Industry

    IF (index_MEGAN==1) THEN
    E_gases(:,ind_NO2) =  E_gases(:,ind_NO2)+EMI(:,17) ! Anthropogenic + biogenic NO2 emissions
    END IF
    E_gases(1,ind_NO) = (0.4*sum(EmissionNO2(4:6,tr))+0.4*sum(EmissionNO2(8:12,tr)))/(dz(1)*100.)
	E_gases(4,ind_NO) = 0.4*EmissionNO2(7,tr)/(dz(4)*100.) ! Ships
	E_gases(6,ind_NO) = 0.4*sum(EmissionNO2(1:3,tr))/(dz(6)*100.) ! Power Plants and industries emit at higher altitudes
	
    E_gases(1,ind_SO2) = (sum(EmissionSO2(4:6,tr))+sum(EmissionSO2(8:12,tr)))/(dz(1)*100.)
	E_gases(4,ind_SO2) = EmissionSO2(7,tr)/(dz(4)*100.) ! Ships
	E_gases(6,ind_SO2) = sum(EmissionSO2(1:3,tr))/(dz(6)*100.) ! Power Plants and industries emit  at higher altitudes
	
    E_gases(1,ind_DMS) = EmissionDMS(tr)/(dz(1)*1D2) ! Emission of DMS
    
	! Added by Pontus 2021-11-26
	! Carpenter et al. (Nature Geoscience, 2013):
	IF (landuse_index(tr) == 6) THEN
	!IF (SST(tr)>273.15) THEN
	!Iion_aq=(0.225D0*(SST(tr)-273.15)**2+19D0)*1D-9 ! mol/dm^3 Estimated from Chance et al., Environ. Sci.-Proc. Imp., 16, 1841–1859 2014
	!ELSE
	!Iion_aq=19D-9
	!END IF
	!ELSE
	!Iion_aq=0D0	
	Iion_aq=1.46D6*EXP(-9134D0/SST(tr)) ! mol/dm^3 Estimated from MacDonald et al. Atmos. Chem. Phys., 14, 5841–5852, https://doi.org/10.5194/acp-14-5841-2014, 2014.
	END IF
	
	O3_ppb=(conc(1,ind_O3)/M(1,tr))*1D9
	EmissionI2=O3_ppb*(Iion_aq**1.3)*(1.74D9-6.54D8*LOG(windspeed(tr)))*1D-4*1D-9*Na/(24D0) ! molec cm^-2 h^-1
	EmissionHOI=O3_ppb*SQRT(Iion_aq)*(3.56D5/windspeed(tr)-2.16D4)*1D-4*1D-9*Na/(24D0) ! molec cm^-2 h^-1
	IF (EmissionI2<0D0) THEN
	EmissionI2=0D0
	END IF
	IF (EmissionHOI<0D0) THEN
	EmissionHOI=0D0
	END IF
	
	! Test to set ionorganic iodine emissions to zero!
	!EmissionHOI=0D0
	!EmissionI2=0D0
	
    ! Added by Pontus 20210105:
    E_gases(1,ind_CH3BR) = EmissionCH3BR(tr)/(dz(1)*1D2) 
    E_gases(1,ind_DIBRET) = EmissionDIBRET(tr)/(dz(1)*1D2) 
    E_gases(1,ind_CHBr3) = EmissionCHBr3(tr)/(dz(1)*1D2) 
    E_gases(1,ind_CH3I) = EmissionCH3I(tr)/(dz(1)*1D2) 
    E_gases(1,ind_C3H7I) = EmissionC3H7I(tr)/(dz(1)*1D2) 
    E_gases(1,ind_C2H5I) = EmissionC2H5I(tr)/(dz(1)*1D2) 
    E_gases(1,ind_CH2ICL) = EmissionCH2ICL(tr)/(dz(1)*1D2) 
    E_gases(1,ind_CH2IBr) = EmissionCH2IBr(tr)/(dz(1)*1D2) 
    E_gases(1,ind_CH2I2) = EmissionCH2I2(tr)/(dz(1)*1D2) 
    E_gases(1,ind_IODINE2) = EmissionI2/(dz(1)*1D2)
	E_gases(1,ind_HOI) = EmissionHOI/(dz(1)*1D2)
	
    E_gases(1,ind_TCE) = EmissionTCE(tr)/(dz(1)*1D2) 
    E_gases(1,ind_TRICLETH) = EmissionTRICLETH(tr)/(dz(1)*1D2) 
    E_gases(1,ind_CHCL3) = EmissionCHCL3(tr)/(dz(1)*1D2) 
    E_gases(1,ind_CH2CL2) = EmissionCH2CL2(tr)/(dz(1)*1D2) 
    E_gases(1,ind_CH3CL) = EmissionCH3CL(tr)/(dz(1)*1D2) 
    !E_gases(1,ind_C2H6) = EmissionC2H6(tr)/(dz(1)*1D2) 
    !E_gases(1,ind_C3H8) = EmissionC3H8(tr)/(dz(1)*1D2) 
    E_gases(1,ind_CH3CHO) = EmissionCH3CHO(tr)/(dz(1)*1D2) 
    E_gases(1,ind_C2H5CHO) = EmissionC2H5CHO(tr)/(dz(1)*1D2) 
    E_gases(1,ind_C5H8) = E_gases(1,ind_C5H8)+EmissionISO(tr)/(dz(1)*1D2) 
    
    ! NMVOC EMEP emissions (to be used with chemistry including anthropogenic emissions) molec cm⁻3 h⁻1 
	!IF (index_organic_chem==1) THEN
	E_gases(1,ind_C2H6) = EmissionC2H6(tr)/(dz(1)*1D2) ! Ethane   
    E_gases(1,ind_NC4H10) = EmissionNC4H10(tr)/(dz(1)*1D2) ! Butane	
    E_gases(1,ind_C2H4) = E_gases(1,ind_C2H4)+EmissionC2H4(tr)/(dz(1)*1D2) ! Etene/ethylene
    E_gases(1,ind_C3H6) = E_gases(1,ind_C3H6)+EmissionC3H6(tr)/(dz(1)*1D2) ! ! Propene/popylene
	E_gases(1,ind_C5H8) = E_gases(1,ind_C5H8)+EmissionC5H8(tr)/(dz(1)*1D2) 
	
    E_gases(1,ind_OXYL) =   0.1*EmissionOXYL(tr)/(dz(1)*1D2)
    E_gases(1,ind_BENZENE)= 0.2*EmissionOXYL(tr)/(dz(1)*1D2) ! benzene
    E_gases(1,ind_PXYL) =   0.1*EmissionOXYL(tr)/(dz(1)*1D2) ! p-xylene
    E_gases(1,ind_TOLUENE)= 0.2*EmissionOXYL(tr)/(dz(1)*1D2) ! toluene
	E_gases(1,ind_MXYL) =   0.1*EmissionOXYL(tr)/(dz(1)*1D2) ! m-xylene
    E_gases(1,ind_TM124B) = 0.1*EmissionOXYL(tr)/(dz(1)*1D2) ! 1,2,4 trimethylbenzene
    E_gases(1,ind_TM135B) = 0.1*EmissionOXYL(tr)/(dz(1)*1D2) ! 1,3,5 trimethylbenzene
    E_gases(1,ind_TM123B) = 0.1*EmissionOXYL(tr)/(dz(1)*1D2) ! 1,2,3 trimethylbenzene
   
	E_gases(1,ind_CH3OH) = E_gases(1,ind_CH3OH)+EmissionCH3OH(tr)/(dz(1)*1D2) ! CH3OH
	E_gases(1,ind_C2H5OH) = E_gases(1,ind_C2H5OH)+EmissionC2H5OH(tr)/(dz(1)*1D2) ! C2H5OH
	E_gases(1,ind_HCHO) = E_gases(1,ind_HCHO)+EmissionHCHO(tr)/(dz(1)*1D2) ! Formaldehyde
	E_gases(1,ind_CH3CHO) = E_gases(1,ind_CH3CHO)+EmissionCH3CHO(tr)/(dz(1)*1D2) ! acetaldehyde
    E_gases(1,ind_MEK) = EmissionMEK(tr)/(dz(1)*1D2) ! MEK (Methyl ethyl ketone)
    E_gases(1,ind_GLYOX) = EmissionGLYOX(tr)/(dz(1)*1D2) ! Glyoxal
    E_gases(1,ind_MGLYOX) = EmissionMGLYOX(tr)/(dz(1)*1D2) ! Methylglyoxal
	
	E_gases(1,ind_NC9H20) =  0.25*EmissionC9_C12_alkanes(tr)/(dz(1)*1D2) ! nonane
    E_gases(1,ind_NC10H22) = 0.25*EmissionC9_C12_alkanes(tr)/(dz(1)*1D2) ! decane
	E_gases(1,ind_NC11H24) = 0.25*EmissionC9_C12_alkanes(tr)/(dz(1)*1D2) ! undecane
	E_gases(1,ind_NC12H26) = 0.25*EmissionC9_C12_alkanes(tr)/(dz(1)*1D2) ! dodecane
    
	 ! NMVOC EMEP emissions (to be used with chemistry including anthropogenic emissions) molec cm⁻3 h⁻1 
	E_gases(6,ind_C2H6) = EmissionC2H6_Pplants_Ind(tr)/(dz(6)*1D2) ! Ethane   
    E_gases(6,ind_NC4H10) = EmissionNC4H10_Pplants_Ind(tr)/(dz(6)*1D2) ! Butane	
    E_gases(6,ind_C2H4) = EmissionC2H4_Pplants_Ind(tr)/(dz(6)*1D2) ! Etene/ethylene
    E_gases(6,ind_C3H6) = EmissionC3H6_Pplants_Ind(tr)/(dz(6)*1D2) ! ! Propene/popylene
	E_gases(6,ind_C5H8) = EmissionC5H8_Pplants_Ind(tr)/(dz(6)*1D2) ! Isoprene
	
    E_gases(6,ind_OXYL) =   0.1*EmissionOXYL_Pplants_Ind(tr)/(dz(6)*1D2)
    E_gases(6,ind_BENZENE)= 0.2*EmissionOXYL_Pplants_Ind(tr)/(dz(6)*1D2) ! benzene
    E_gases(6,ind_PXYL) =   0.1*EmissionOXYL_Pplants_Ind(tr)/(dz(6)*1D2) ! p-xylene
    E_gases(6,ind_TOLUENE)= 0.2*EmissionOXYL_Pplants_Ind(tr)/(dz(6)*1D2) ! toluene
	E_gases(6,ind_MXYL) =   0.1*EmissionOXYL_Pplants_Ind(tr)/(dz(6)*1D2) ! m-xylene
    E_gases(6,ind_TM124B) = 0.1*EmissionOXYL_Pplants_Ind(tr)/(dz(6)*1D2) ! 1,2,4 trimethylbenzene
    E_gases(6,ind_TM135B) = 0.1*EmissionOXYL_Pplants_Ind(tr)/(dz(6)*1D2) ! 1,3,5 trimethylbenzene
    E_gases(6,ind_TM123B) = 0.1*EmissionOXYL_Pplants_Ind(tr)/(dz(6)*1D2) ! 1,2,3 trimethylbenzene
   
	E_gases(6,ind_CH3OH) = EmissionCH3OH_Pplants_Ind(tr)/(dz(6)*1D2) ! CH3OH
	E_gases(6,ind_C2H5OH) = EmissionC2H5OH_Pplants_Ind(tr)/(dz(6)*1D2) ! C2H5OH
	E_gases(6,ind_HCHO) = EmissionHCHO_Pplants_Ind(tr)/(dz(6)*1D2) ! Formaldehyde
	E_gases(6,ind_CH3CHO) = EmissionCH3CHO_Pplants_Ind(tr)/(dz(6)*1D2) ! acetaldehyde
    E_gases(6,ind_MEK) = EmissionMEK_Pplants_Ind(tr)/(dz(6)*1D2) ! MEK (Methyl ethyl ketone)
    E_gases(6,ind_GLYOX) = EmissionGLYOX_Pplants_Ind(tr)/(dz(6)*1D2) ! Glyoxal
    E_gases(6,ind_MGLYOX) = EmissionMGLYOX_Pplants_Ind(tr)/(dz(6)*1D2) ! Methylglyoxal
	
	
	E_gases(4,ind_C2H6) = EmissionC2H6_Ship(tr)/(dz(4)*1D2) ! Ethane   
    E_gases(4,ind_NC4H10) = EmissionNC4H10_Ship(tr)/(dz(4)*1D2) ! Butane	
    E_gases(4,ind_C2H4) = EmissionC2H4_Ship(tr)/(dz(4)*1D2) ! Etene/ethylene
    E_gases(4,ind_C3H6) = EmissionC3H6_Ship(tr)/(dz(4)*1D2) ! ! Propene/popylene
	E_gases(4,ind_C5H8) = EmissionC5H8_Ship(tr)/(dz(4)*1D2) ! Isoprene
	
    E_gases(4,ind_OXYL) =   0.1*EmissionOXYL_Ship(tr)/(dz(4)*1D2)
    E_gases(4,ind_BENZENE)= 0.2*EmissionOXYL_Ship(tr)/(dz(4)*1D2) ! benzene
    E_gases(4,ind_PXYL) =   0.1*EmissionOXYL_Ship(tr)/(dz(4)*1D2) ! p-xylene
    E_gases(4,ind_TOLUENE)= 0.2*EmissionOXYL_Ship(tr)/(dz(4)*1D2) ! toluene
	E_gases(4,ind_MXYL) =   0.1*EmissionOXYL_Ship(tr)/(dz(4)*1D2) ! m-xylene
    E_gases(4,ind_TM124B) = 0.1*EmissionOXYL_Ship(tr)/(dz(4)*1D2) ! 1,2,4 trimethylbenzene
    E_gases(4,ind_TM135B) = 0.1*EmissionOXYL_Ship(tr)/(dz(4)*1D2) ! 1,3,5 trimethylbenzene
    E_gases(4,ind_TM123B) = 0.1*EmissionOXYL_Ship(tr)/(dz(4)*1D2) ! 1,2,3 trimethylbenzene
   
	E_gases(4,ind_CH3OH) = EmissionCH3OH_Ship(tr)/(dz(4)*1D2) ! CH3OH
	E_gases(4,ind_C2H5OH) = EmissionC2H5OH_Ship(tr)/(dz(4)*1D2) ! C2H5OH
	E_gases(4,ind_HCHO) = EmissionHCHO_Ship(tr)/(dz(4)*1D2) ! Formaldehyde
	E_gases(4,ind_CH3CHO) = EmissionCH3CHO_Ship(tr)/(dz(4)*1D2) ! acetaldehyde
    E_gases(4,ind_MEK) = EmissionMEK_Ship(tr)/(dz(4)*1D2) ! MEK (Methyl ethyl ketone)
    E_gases(4,ind_GLYOX) = EmissionGLYOX_Ship(tr)/(dz(4)*1D2) ! Glyoxal
    E_gases(4,ind_MGLYOX) = EmissionMGLYOX_Ship(tr)/(dz(4)*1D2) ! Methylglyoxal
	
    
    ! !! Biomass burning (wild fires) emissions:
    ! E_gases(1,ind_CO)=E_gases(1,ind_CO)+DM_emissions_biomass_burning(tr)*127D0/1D3/28D0/1D4*Na/(dz(1)*100.) ! molec cm^-3 h^-1 
    ! E_gases(1,ind_NO)=E_gases(1,ind_NO)+DM_emissions_biomass_burning(tr)*0.9D0/1D3/30D0/1D4*Na/(dz(1)*100.) ! molec cm^-3 h^-1
    ! E_gases(1,ind_SO2)=E_gases(1,ind_SO2)+DM_emissions_biomass_burning(tr)*1.1D0/1D3/64D0/1D4*Na/(dz(1)*100.) ! molec cm^-3 h^-1
    ! E_gases(1,ind_C2H6)=E_gases(1,ind_C2H6)+DM_emissions_biomass_burning(tr)*1.79D0/1D3/30D0/1D4*Na/(dz(1)*100.) ! molec cm^-3 h^-1
    ! E_gases(1,ind_CH3OH)=E_gases(1,ind_CH3OH)+DM_emissions_biomass_burning(tr)*2.82D0/1D3/32D0/1D4*Na/(dz(1)*100.) ! molec cm^-3 h^-1
    ! E_gases(1,ind_C2H5OH)=E_gases(1,ind_C2H5OH)+DM_emissions_biomass_burning(tr)*0.055D0/1D3/46D0/1D4*Na/(dz(1)*100.) ! molec cm^-3 h^-1
    ! E_gases(1,ind_C3H8)=E_gases(1,ind_C3H8)+DM_emissions_biomass_burning(tr)*0.44D0/1D3/44D0/1D4*Na/(dz(1)*100.) ! molec cm^-3 h^-1
    ! E_gases(1,ind_C2H4)=E_gases(1,ind_C2H4)+DM_emissions_biomass_burning(tr)*1.42D0/1D3/28D0/1D4*Na/(dz(1)*100.) ! molec cm^-3 h^-1
    ! E_gases(1,ind_C3H6)=E_gases(1,ind_C3H6)+DM_emissions_biomass_burning(tr)*1.13D0/1D3/44D0/1D4*Na/(dz(1)*100.) ! molec cm^-3 h^-1
    ! E_gases(1,ind_C5H8)=E_gases(1,ind_C5H8)+DM_emissions_biomass_burning(tr)*0.15D0/1D3/68D0/1D4*Na/(dz(1)*100.) ! molec cm^-3 h^-1
    ! E_gases(1,ind_BENZENE)=E_gases(1,ind_BENZENE)+DM_emissions_biomass_burning(tr)*1.11D0/1D3/78D0/1D4*Na/(dz(1)*100.) ! molec cm^-3 h^-1
    ! E_gases(1,ind_TOLUENE)=E_gases(1,ind_TOLUENE)+DM_emissions_biomass_burning(tr)*0.48D0/1D3/92D0/1D4*Na/(dz(1)*100.) ! molec cm^-3 h^-1
    ! E_gases(1,ind_MXYL)=E_gases(1,ind_MXYL)+DM_emissions_biomass_burning(tr)*0.18D0/1D3/106D0/1D4*Na/(dz(1)*100.) ! molec cm^-3 h^-1
    ! E_gases(1,ind_HCHO)=E_gases(1,ind_HCHO)+DM_emissions_biomass_burning(tr)*1.86D0/1D3/30D0/1D4*Na/(dz(1)*100.) ! molec cm^-3 h^-1
    ! E_gases(1,ind_CH3CHO)=E_gases(1,ind_CH3CHO)+DM_emissions_biomass_burning(tr)*0.77D0/1D3/44D0/1D4*Na/(dz(1)*100.) ! molec cm^-3 h^-1
    ! E_gases(1,ind_C2H5CHO)=E_gases(1,ind_C2H5CHO)+DM_emissions_biomass_burning(tr)*0.75D0/1D3/58D0/1D4*Na/(dz(1)*100.) ! molec cm^-3 h^-1
    ! E_gases(1,ind_HCOOH)=E_gases(1,ind_HCOOH)+DM_emissions_biomass_burning(tr)*0.57D0/1D3/46D0/1D4*Na/(dz(1)*100.) ! molec cm^-3 h^-1
    ! E_gases(1,ind_CH3CO2H)=E_gases(1,ind_CH3CO2H)+DM_emissions_biomass_burning(tr)*4.41D0/1D3/60D0/1D4*Na/(dz(1)*100.) ! molec cm^-3 h^-1
    ! E_gases(1,ind_MEK)=E_gases(1,ind_MEK)+DM_emissions_biomass_burning(tr)*0.22D0/1D3/72D0/1D4*Na/(dz(1)*100.) ! molec cm^-3 h^-1
    ! E_gases(1,ind_MGLYOX)=E_gases(1,ind_MGLYOX)+DM_emissions_biomass_burning(tr)*0.73D0/1D3/72D0/1D4*Na/(dz(1)*100.) ! molec cm^-3 h^-1
    ! E_gases(1,ind_APINENE)=E_gases(1,ind_APINENE)+DM_emissions_biomass_burning(tr)*1.2D0/1D3/136D0/1D4*Na/(dz(1)*100.) ! molec cm^-3 h^-1
    ! E_gases(1,ind_BPINENE)=E_gases(1,ind_BPINENE)+DM_emissions_biomass_burning(tr)*0.8D0/1D3/136D0/1D4*Na/(dz(1)*100.) ! molec cm^-3 h^-1
    ! E_gases(1,ind_BUT1ENE)=E_gases(1,ind_BUT1ENE)+DM_emissions_biomass_burning(tr)*0.385D0/1D3/64D0/1D4*Na/(dz(1)*100.) ! molec cm^-3 h^-1
    ! E_gases(1,ind_NC4H10)=E_gases(1,ind_NC4H10)+DM_emissions_biomass_burning(tr)*0.349D0/1D3/66D0/1D4*Na/(dz(1)*100.) ! molec cm^-3 h^-1    
     
	! END IF

    !-------------------------------------------------------------!
    ! Dry deposition gases (SO2,O3,NO2,NO,HNO3,H2O2,HCHO,         !
    ! CH3OOH,HONO,H2SO4)                                          !
    !-------------------------------------------------------------!
IF (dep_g_index == 1) THEN
psat_org=10.**(A_Nannoolal-B_Nannoolal/T2m(tr)) ! Vapor pressures condensable organic comp (atm).
psat_org(635:742)=psat_org(635:742)
psat_org(635:645)=1D-100 ! RO2
psat_org(756:761)=1D-100 ! RO2
Henrys_coeff=55.56*1/(y_org_water*psat_org) ! M/atm
Dorg=1D-7*T2m(tr)**1.75D0*SQRT(1/(Mair*1D3)+1/mcm_prop(1,:))/&
     (Ps(1,tr)/1.01325D5*(VolX**(1D0/3D0)+20.1**(1D0/3D0))**2D0) ! Fuller's method /(Tang et al., ACP 2015) ! Diffusivity organic compounds (m^2 s^-1)

      CALL dry_dep_gases(T2m(tr),P2m(tr),RH2m(tr),fric_veloc(tr),PBLH(tr),SHTFs(tr),DSWF(tr),&
        landuse_index(tr),vd_gas,rain(tr),windspeed(tr),season_index,Henrys_coeff,Dorg)
  
     index_dep_gases = (/ ind_SO2,ind_O3,ind_NO2,ind_NO,ind_HNO3,ind_H2O2,ind_HCHO,&
         ind_CH3OOH,ind_HONO,ind_H2SO4,ind_HCl /) ! index in gas-phase vector of gases in dry dep subroutine
    
	IF (index_no_O3_dep==1) THEN
	vd_gas(2)=1D-100 ! Skip depostion of O3
	END IF
	
    IF (landuse_index(tr)==6 .OR. landuse_index(tr)==8) THEN
         conc(1,index_dep_gases) = conc(1,index_dep_gases)*exp(-vd_gas(1:11)*dt/dz(1))
         !cHCl(1) = cHCl(1)*exp(-vd_gas(11)*dt/dz(1)) ! Update HCl gas conc after dry dep
         
		 !cNH3(1) = cNH3(1)*exp(-vd_gas(12)*dt/dz(1)) ! Update NH3 gas conc after dry dep
		 conc(1,ind_NH3) = conc(1,ind_NH3)*exp(-vd_gas(12)*dt/dz(1)) ! Update NH3 gas conc after dry dep
         conc(1,ind_DMA) = conc(1,ind_DMA)*exp(-vd_gas(12)*dt/dz(1)) ! Update DMA gas conc after dry dep
		 !conc(1,index_cond)=conc(1,index_cond)*exp(-vd_gas(13)*dt/dz(1))   ! Update ELVOC monomer gas conc after dry dep
         conc(1,index_cond)=conc(1,index_cond)*exp(-vd_gas(13:NCOND+12)*dt/dz(1))   ! Other SVOCs etc 
		 conc(1,ind_HIO3) = conc(1,ind_HIO3)*exp(-vd_gas(10)*dt/dz(1)) ! Assume the same dry deposition rate of HIO3 as H2SO4
		 conc(1,ind_MSA) = conc(1,ind_MSA)*exp(-vd_gas(10)*dt/dz(1)) ! Assume the same dry deposition rate of MSA as H2SO4
 
    ELSE 
         conc(1,index_dep_gases) = conc(1,index_dep_gases)*exp(-vd_gas(1:11)*dt/z(4))
         conc(2,index_dep_gases) = conc(2,index_dep_gases)*exp(-vd_gas(1:11)*dt/z(4))
         conc(3,index_dep_gases) = conc(3,index_dep_gases)*exp(-vd_gas(1:11)*dt/z(4))
         
         !cHCl(1:3) = cHCl(1:3)*exp(-vd_gas(11)*dt/z(4)) ! Update HCl gas conc after dry dep
         conc(1:3,ind_NH3) = conc(1:3,ind_NH3)*exp(-vd_gas(12)*dt/z(4)) ! Update NH3 gas conc after dry dep
		 conc(1:3,ind_DMA) = conc(1:3,ind_DMA)*exp(-vd_gas(12)*dt/z(4)) ! Update DMA gas conc after dry dep
  
         conc(1,index_cond)=conc(1,index_cond)*exp(-vd_gas(13:NCOND+12)*dt/z(4))   ! Other SVOCs etc 
         conc(2,index_cond)=conc(2,index_cond)*exp(-vd_gas(13:NCOND+12)*dt/z(4))   ! Other SVOCs etc 
         conc(3,index_cond)=conc(3,index_cond)*exp(-vd_gas(13:NCOND+12)*dt/z(4))   ! Other SVOCs etc
		 
		 conc(1,ind_HIO3) = conc(1,ind_HIO3)*exp(-vd_gas(10)*dt/z(4))  ! Assume the same dry deposition rate of HIO3 as H2SO4
		 conc(2,ind_HIO3) = conc(2,ind_HIO3)*exp(-vd_gas(10)*dt/z(4))  ! Assume the same dry deposition rate of HIO3 as H2SO4
		 conc(3,ind_HIO3) = conc(3,ind_HIO3)*exp(-vd_gas(10)*dt/z(4))  ! Assume the same dry deposition rate of HIO3 as H2SO4
		 
		 conc(1,ind_MSA) = conc(1,ind_MSA)*exp(-vd_gas(10)*dt/z(4))  ! Assume the same dry deposition rate of MSA as H2SO4
		 conc(2,ind_MSA) = conc(2,ind_MSA)*exp(-vd_gas(10)*dt/z(4))  ! Assume the same dry deposition rate of MSA as H2SO4
		 conc(3,ind_MSA) = conc(3,ind_MSA)*exp(-vd_gas(10)*dt/z(4))  ! Assume the same dry deposition rate of MSA as H2SO4

		 
    END IF
END IF    

! Consider wet scaveninging of water soluble gases:
IF (rain(tr)>0D0) THEN
CALL wet_deposition_gas(rain(tr),cloudLWC(:,tr),conc,dt)
END IF

! Ocean NH3 air-sea flux of NH3 (see Wentworth et al., Atmos. Chem. Phys., 16, 1937–1953, 2016) 
IF (landuse_index(tr)==6) THEN 
! New NH3 ocean flux estimates Robin and Pontus 2021-12-02, based on Paulot et al., 10.1002/2015GB005106: 
! New CH3NO3 and C2H5NO3 flux Fisher et al. Journal of Geophysical Research: Atmospheres, 123, 12,429–12,451. https://doi.org/10.1029/2018JD029046

cNHx_ocean=0.2D-3 ! Ocean total NH conc southern ocean mol/m^3
H_NH3_sea=((17.93*SST(tr)/273.15)*EXP(4092D0/SST(tr)-9.7D0))**(-1) ! dimensionless effective gas-over-liquid Henry’s law constant
S_sea=3.5D-2 ! Ocean salinity mass 
pKa_NH3_sea=10.04-3.16D-2*(SST(tr)-273.15)+3.1D-3*S_sea
pH_ocean=8.1
H_prim_NH3_sea=H_NH3_sea/(1D0+10**(-pH_ocean+pKa_NH3_sea))
X_NH3=H_prim_NH3_sea*cNHx_ocean*Na ! Compensation point, surface ocean equilibrium [NH3(g)] molec/m^3 (~1 nmol/m^3)
!kg_NH3=windspeed(tr)/(770D0+45D0*17.03**(1D0/3D0)) ! Air-side transfer velocoty m/s

! DMA Section : Based on Carpenter (Ocean-atmosphere trace gas exchange) saying that pKa can be estimated similar to NH3
! ---------- Test 1 : Just use same approach as for MeNO3 and EtNO3
cDMA_ocean     = 9D0 * 1D-9 * 1D3 ! nmol/L -> mol/m3  :   Carpenter (Ocean-atmosphere trace gas exchange) Massachusetts DMA ocean surface conc.
H_DMA_sea      = 1D0/(3D-1*EXP(4000D0*(1D0/SST(tr)-1D0/298.15D0))*Rg*298.15) !From H(cp) to H(cc) multiply by RT !*Rg*298.15) ! Sander et al., 2015 (from Sander paper, similar to e.g. MeNo3)
pKa_DMA_sea=11.50-3.16D-2*(SST(tr)-273.15)+3.1D-3*S_sea ! Estimated based on CRC Handbook and T dependence as NH3 
H_prim_DMA_sea = H_DMA_sea/(1D0+10**(-pH_ocean+pKa_DMA_sea)) ! Use the same pKa dependance as for ammonia (according to carpenter)
X_DMA          = H_prim_DMA_sea*cDMA_ocean*Na ! Compensation point, surface ocean equilibrium conc molec/m^3

cMeNO3_ocean=200D-9 ! Southern Ocean upper concentration estimate (mol/m^3) 
cEtNO3_ocean=cMeNO3_ocean/6D0 ! Southern Ocean  mol/m^3
! Dimensionless effective gas-over-liquid Henry’s law constants:
H_MeNO3=1D0/(2D-2*EXP(4700D0*(1D0/SST(tr)-1D0/298.15D0))*Rg*298.15) ! Sander et al., 2015
H_EtNO3=1D0/(1.6D-2*EXP(5400D0*(1D0/SST(tr)-1D0/298.15D0))*Rg*298.15) ! Sander et al., 2015

X_MeNO3=H_MeNO3*cMeNO3_ocean*Na ! Compensation point, surface ocean equilibrium conc molec/m^3 
X_EtNO3=H_EtNO3*cEtNO3_ocean*Na ! Compensation point, surface ocean equilibrium conc molec/m^3 

Sc_NH3=0.57D0
Sc_MeNO3=1D0
Sc_EtNO3=1D0
Sc_DMA   = 1D0

C_D=(fric_veloc(tr)/windspeed(tr))**2D0
kg_NH3=1D-3+fric_veloc(tr)/(13.3*Sc_NH3+sqrt(1D0/C_D)-5D0+LOG(Sc_NH3)/(2D0*0.4D0))
kg_MeNO3=1D-3+fric_veloc(tr)/(13.3*Sc_MeNO3+sqrt(1D0/C_D)-5D0+LOG(Sc_MeNO3)/(2D0*0.4D0))
kg_EtNO3=1D-3+fric_veloc(tr)/(13.3*Sc_EtNO3+sqrt(1D0/C_D)-5D0+LOG(Sc_EtNO3)/(2D0*0.4D0))
kg_DMA   = 1D-3+fric_veloc(tr)/(13.3*Sc_DMA+sqrt(1D0/C_D)-5D0+LOG(Sc_DMA)/(2D0*0.4D0))

F_NH3ocean=kg_NH3*(X_NH3-conc(1,ind_NH3)*1D6) ! NH3 flux molec m^-2 s^-1
F_MeNO3ocean=kg_MeNO3*(X_MeNO3-conc(1,ind_CH3NO3)*1D6) ! CH3NO3 flux molec m^-2 s^-1
F_EtNO3ocean=kg_EtNO3*(X_EtNO3-conc(1,ind_C2H5NO3)*1D6) ! C2H5NO3 flux molec m^-2 s^-1
F_DMAocean   = kg_DMA*(X_DMA-conc(1,ind_DMA)*1D6) ! C2H5NO3 flux molec m^-2 s^-1


conc(1,ind_NH3)=MAXVAL((/conc(1,ind_NH3)+1D-6*F_NH3ocean/dz(1)*dt, 1D3/)) 
conc(1,ind_CH3NO3)=MAXVAL((/conc(1,ind_CH3NO3)+1D-6*F_MeNO3ocean/dz(1)*dt, 1D3/)) 
conc(1,ind_C2H5NO3)=MAXVAL((/conc(1,ind_C2H5NO3)+1D-6*F_EtNO3ocean/dz(1)*dt, 1D3/))
conc(1,ind_DMA)=MAXVAL((/conc(1,ind_DMA)+1D-6*F_DMAocean/dz(1)*dt, 1D3/))

! kg_MeNO3=windspeed(tr)/(770D0+45D0*17.03**(1D0/3D0)) ! Air-side transfer velocoty m/s

! F_NH3ocean=kg_MeNO3*(X_MeNO3-conc(1,ind_CH3NO3)*1D6) ! NH3 flux molec m^-2 s^-1

END IF

! Barbero, et al (2021) Journal of Geophysical Research: Atmospheres, 126, e2021JD035062. https://doi.org/10.1029/2021JD035062, from Dome C
! Jones et al. (2001) GEOPHYSICAL RESEARCH LETTERS, VOL. 28, NO. 8, PAGES 1499-1502, APRIL 15, 2001, from Neumayer  
! NOx emissions from snow covered sourceses on Antarctica
IF (landuse_index(tr)==8) THEN
! Estimated NOx emission flux from snow covered surfaces on Antarctica:
!F_NOxsnow=(DSWF(tr)/2D2)*5D12 ! molec m^-2 s^-1 Domes C
F_NOxsnow=(DSWF(tr)/8D2)*1D12 ! molec m^-2 s^-1 Neumayer estimates
conc(1,ind_NO2)=conc(1,ind_NO2)+1D-6*F_NOxsnow/dz(1)*dt
END IF
IF (landuse_index(tr)==6) THEN
! Estimated NOx emission flux from open ocean (Bräuer et al., 2013; Thompson and Zafiriou, 1983; Tian et al Ocean Sci., 16, 135–148, 2020):
!F_NOxsnow=(DSWF(tr)/1D2)*5D12 ! molec m^-2 s^-1
!F_NOxocean=(DSWF(tr)/2D2)*2D12 ! molec m^-2 s^-1
F_NOxocean=(DSWF(tr)/8D2)*1D11 ! molec m^-2 s^-1
conc(1,ind_NO2)=conc(1,ind_NO2)+1D-6*F_NOxocean/dz(1)*dt ! 1D-6*F_NOxsnow/dz(1)*dt
END IF

!-----------------------------------------------------------------
! New particle formation
!  IF (nucleation_index==1) THEN
 
!  DO j=1,Nz
!  c_acid=conc(j,ind_H2SO4)*1D6 
!  c_base=conc(j,ind_NH3)*1D6 
 
!  ipr=q_ion(j)*1D6 ! Ion production rate (ions/m^3/s)
 
!  c_org=0D0 ! total ELVOC concentration (molec m^-3)
!  solve_ss=.false. ! .true. if the steady state assumption is used
!  !CALL get_acdc_J(c_acid,c_base,c_org,CS_H2SO4(j),Ts(j,tr),dt,solve_ss,Jnucl(j),diameter_acdc)
!  c_clusters1=c_clusters(j,:)
!  CALL  get_acdc_J(c_acid,c_base,c_org,CS_H2SO4(j),Ts(j,tr),ipr,dt,solve_ss,Jnucl_N(j),diameter_acdc,Nuc_by_charge,c_clusters1)
! !  N_bins(j,1) = N_bins(j,1)+Jnucl(j)*dt
! !  c_p(j,:,1) = c_p(j,:,1)+c_p_nucl*Jnucl(j)*1D-6*dt
! ! Update cluster and vapour concentrations!
!  c_clusters(j,:)=c_clusters1
!  conc(j,ind_H2SO4)=c_acid/1D6
!  conc(j,ind_NH3)=c_base/1D6
 
!  c_acid=conc(j,ind_H2SO4)*1D6 
!  c_base=conc(j,ind_DMA)*1D6 
!  c_clusters2=c_clustersDMA(j,:)
!  CALL  get_acdc_D(c_acid,c_base,c_org,CS_H2SO4(j),Ts(j,tr),dt,solve_ss,Jnucl_D(j),diameter_acdcDMA,c_clusters2)
!  N_bins(j,1) = N_bins(j,1)+Jnucl_N(j)*dt
!  N_bins(j,2) = N_bins(j,2)+Jnucl_D(j)*dt
!  c_p(j,:,1) = c_p(j,:,1)+c_p_nucl*Jnucl_N(j)*1D-6*dt
!  c_p(j,:,2) = c_p(j,:,2)+c_p_nuclDMA*Jnucl_D(j)*1D-6*dt
!  ! Update cluster and vapour concentrations!
!  c_clusters(j,:)=c_clusters1
!  c_clustersDMA(j,:)=c_clusters2
!  conc(j,ind_H2SO4)=c_acid/1D6
!  conc(j,ind_DMA)=c_base/1D6
!  END DO        
!  END IF
 

IF(nucleation_index==1) THEN


    
    ! write(*,*) 'l2089 before nucl call',  Jnucl_N(1)*1d-6,Jnucl_D(1)*1d-6,'D:',sum(conc(:, ind_DMA),dim=1),'A:',sum(conc(:, ind_H2SO4),dim=1),'N:', &
    ! sum(conc(:, ind_NH3),dim=1)
    write(*,*) 'l2144 nbins ', sum(N_bins(1,1:10))!, sum(chem_1%conc_coag_clust), sum(chem_2%conc_coag_clust), sum(comp_evap(:,1),dim=1), sum(comp_evap(:,4),dim=1), sum(comp_evap(:,11),dim=1)
    ! cHIO2=1D5
    ! conc(:,ind_HIO3)=1D6
    write(*,*) 'l2146 [HIO2] = ', conc(1,ind_HIO2),'',conc(1,ind_HIO3),'',conc(1,ind_H2SO4),'',conc(1,ind_DMA),'',conc(1,ind_NH3)
    
    DO j=1, nz
        ipr=q_ion(j)*1D6 ! Ion production rate (ions/m^3/s)
        m_p=(dens_p(j,:)*pi*d_p(j,:)**3D0)/6D0
        
      
        c_clusters1=c_clusters_n(j,:)
        c_clusters2=c_clusters_d(j,:)
        c_clusters3=c_clusters_i(j,:)
        N_bins1=N_bins(j,:)
        

        if (clust_firstcall) then
                call clustering_subroutine(chem_1, chem_2, chem_3, clust_firstcall, n_clustering_vapors, nr_bins, conc(j,ind_H2SO4), & 
                conc(j,ind_dma), conc(j,ind_NH3),conc(j,ind_HIO3), conc(j,ind_HIO2), Jnucl_N(j), Jnucl_D(j), Jnucl_I(j), &
                diameter_acdc,diameter_acdcDMA,diameter_acdcHIO3, CS_H2SO4(j),Ts(j,tr),Ps(j,tr),ipr,dt, d, d_p(j,:),&
                m_p, N_bins1,Mx, qX, n_evap(j), comp_evap(j,:), clustering_systems, c_clusters1, c_clusters2, c_clusters3,j,&
                l_condensation_evap,l_coagulation_loss)

                !!! check the right bins for the outgrowing clusters to be placed in
                !!!! AN system
                Mx_chem1(1) =Mx(1)
                Mx_chem1(2) =Mx(4)
                qX_chem1(1) =qX(1)
                qX_chem1(2) =qX(4)

                !!AD system
                Mx_chem2(1)=Mx(1)
                Mx_chem2(2)=Mx(11)
                qX_chem2(1) =qX(1)
                qX_chem2(2) =qX(11)
                
                !!IiIo system
                Mx_chem3(1)=Mx(10)
                Mx_chem3(2)=Mx(12)
                qX_chem3(1) =qX(10)
                qX_chem3(2) =qX(12) !! assumed same as HIO3

                do kk=1, chem_1%nclust_out
                    ! volume_chem1=real(chem_1%clust_out_molec(kk,:),KIND=dp)
                   
                    volume_chem1 = sum(real(chem_1%clust_out_molec(kk,:),KIND=dp)/Na*Mx_chem1(:)/qX_chem1(:))
                    ! write(*,*) kk, Mx_chem1, qX_chem1
                    DO ii=2,nr_bins
                        IF (volume_chem1 .LT. vp(ii)) THEN
                            chem_1%ind_out_bin(kk)=ii-1
                            EXIT
                        END IF
                    END DO
                    IF (chem_1%ind_out_bin(kk) .EQ. 0) THEN
                        WRITE(*,*) 'Could not find aerosol bin for outgrown cluster of composition',chem_1%clust_out_molec(i,:)
                        STOP
       
                    END IF
                end do

                do kk=1, chem_2%nclust_out
                    volume_chem2 = sum(real(chem_2%clust_out_molec(kk,:),KIND=dp)/Na*Mx_chem2(:)/qX_chem2(:))
                    DO ii=2,nr_bins
                        IF (volume_chem2 .LT. vp(ii)) THEN
                            chem_2%ind_out_bin(kk)=ii-1
                            EXIT
                        END IF
                    END DO
                    IF (chem_2%ind_out_bin(kk) .EQ. 0) THEN
                        WRITE(*,*) 'Could not find aerosol bin for outgrown cluster of composition',chem_2%clust_out_molec(i,:)
                        STOP
                  
                    END IF
                end do
                
                do kk=1, chem_3%nclust_out
                    volume_chem3 = sum(real(chem_3%clust_out_molec(kk,:),KIND=dp)/Na*Mx_chem3(:)/qX_chem3(:))
                    DO ii=2,nr_bins
                        IF (volume_chem3 .LT. vp(ii)) THEN
                            chem_3%ind_out_bin(kk)=ii-1
                            EXIT
                        END IF
                    END DO
                    IF (chem_3%ind_out_bin(kk) .EQ. 0) THEN
                        WRITE(*,*) 'Could not find aerosol bin for outgrown cluster of composition',chem_3%clust_out_molec(i,:)
                        STOP
                  
                    END IF
                end do
                
                ! clust_firstcall=.FALSE.             
                if (l_coagulation_loss) then
                    write(*,*) 'in here'
                    DO i=1,chem_1%nclust_syst
                        chem_1%c_p_clust(1,i) = REAL(chem_1%clust_molec(i,1),KIND=dp)
                        chem_1%c_p_clust(4,i) = REAL(chem_1%clust_molec(i,2),KIND=dp)
                        chem_1%v_clust(i)=SUM(chem_1%c_p_clust(index_dry,i)/Na*MX(index_dry)/qX(index_dry)) ! m^3
                    END DO
                    
                    
                    DO i=1,chem_2%nclust_syst
                        chem_2%c_p_clust(1,i) = REAL(chem_2%clust_molec(i,1),KIND=dp)
                        chem_2%c_p_clust(11,i) = REAL(chem_2%clust_molec(i,2),KIND=dp)
                        chem_2%v_clust(i)=SUM(chem_2%c_p_clust(index_dry,i)/Na*MX(index_dry)/qX(index_dry)) ! m^3
                    END DO
                    
                    DO i=1,chem_3%nclust_syst
                        chem_3%c_p_clust(10,i) = REAL(chem_3%clust_molec(i,1),KIND=dp)
                        chem_3%c_p_clust(12,i) = REAL(chem_3%clust_molec(i,2),KIND=dp)
                        chem_3%v_clust(i)=SUM(chem_3%c_p_clust(index_dry,i)/Na*MX(index_dry)/qX(index_dry)) ! m^3
                    END DO
                end if

                write(*,*) 'sum(chem_1%clust_molec)',sum(chem_1%v_clust),sum(chem_2%v_clust),sum(chem_3%v_clust)

        else
                ! call clustering_subroutine(chem_1, chem_2, clust_firstcall, n_clustering_vapors, nr_bins, conc(j,ind_H2SO4), conc(j,ind_dma), conc(j,ind_NH3), Jnucl_N(j), Jnucl_D(j), &
                ! diameter_acdc,diameter_acdcDMA, CS_H2SO4(j),Ts(j,tr),Ps(j,tr),ipr,dt, d, d_p(j,:), m_p, N_bins1,Mx, qX,n_evap(j), comp_evap(j,:), clustering_systems, c_clusters1,  c_clusters2,j,&
                ! l_condensation_evap,l_coagulation_loss)

                call clustering_subroutine(chem_1, chem_2, chem_3, clust_firstcall, n_clustering_vapors, nr_bins, conc(j,ind_H2SO4), & 
                conc(j,ind_dma), conc(j,ind_NH3),conc(j,ind_HIO3), conc(j,ind_HIO2), Jnucl_N(j), Jnucl_D(j), Jnucl_I(j), &
                diameter_acdc,diameter_acdcDMA,diameter_acdcHIO3, CS_H2SO4(j),Ts(j,tr),Ps(j,tr),ipr,dt, d, d_p(j,:),&
                m_p, N_bins1,Mx, qX, n_evap(j), comp_evap(j,:), clustering_systems, c_clusters1, c_clusters2, c_clusters3,j,&
                l_condensation_evap,l_coagulation_loss)

        end if 

        ! write(*,*) 'l2176 after nucl call ',  Jnucl_N(j)*1d-6,Jnucl_D(j)*1d-6,sum(chem_1%conc_out_all(1:chem_1%nclust_syst)),sum(chem_2%conc_out_all(1:chem_2%nclust_syst))
        !!! scavenging to bin i of clustering molecules
      
          
       
        if (l_condensation_evap) then
            
            c_p(j,1,:) = c_p(j,1,:)+(chem_1%conc_coag_molec(:,1)+chem_2%conc_coag_molec(:,1))*1D-6 !! H2SO4
            c_p(j,4,:) = c_p(j,4,:)+chem_1%conc_coag_molec(:,2)*1D-6 !!! NH3
            c_p(j,11,:) = c_p(j,11,:)+chem_2%conc_coag_molec(:,2)*1D-6  !!! DMA
            c_p(j,10,:) = c_p(j,10,:)+chem_3%conc_coag_molec(:,1)*1D-6  !!! HIO3
            c_p(j,12,:) = c_p(j,12,:)+chem_3%conc_coag_molec(:,2)*1D-6  !!! HIO2

            N_bins1=N_bins(j,:); V_bins1=V_bins(j,:); d_p1=d_p(j,:); dp_dry1=dp_dry(j,:)
            dens_p1=dens_p(j,:); c_p1=c_p(j,:,:)

            CALL full_stationary_rebinning(N_bins1,V_bins1,d_p1,dp_dry1,MX,qX,dens_p1, &
                c_p1,c_p_backg,dt,n_evap(j),comp_evap(j,:))
          
            N_bins(j,:)=N_bins1; V_bins(j,:)=V_bins1; d_p(j,:)=d_p1; dp_dry(j,:)=dp_dry1
            dens_p(j,:)=dens_p1; c_p(j,:,:)=c_p1

        end if

        if (l_coagulation_loss) then
           
            N_bins1=N_bins(j,:); V_bins1=V_bins(j,:); c_p1=c_p(j,:,:); d_p1=d_p(j,:)
            dp_dry1=dp_dry(j,:); dens_p1=dens_p(j,:)
    
            CALL coagulation_clusterin(N_bins1,V_bins1,c_p1,c_p_backg, d_p1,dp_dry1,dt,dens_p1,MX,qX,chem_1,chem_2,chem_3,n_clustering_vapors)
           
            N_bins(j,:)=N_bins1; V_bins(j,:)=V_bins1; c_p(j,:,:)=c_p1; d_p(j,:)=d_p1
            dp_dry(j,:)=dp_dry1; dens_p(j,:)=dens_p1   

        end if


       

        !!!! update cluster and vapor concentration here
        c_clusters_n(j,:)=c_clusters1
        c_clusters_d(j,:)=c_clusters2
        c_clusters_i(j,:)=c_clusters3
    

        DO ii=1,chem_1%nclust_out
            !!! for chemistry 1 H2SO4-NH3
            
            i = chem_1%ind_out_bin(ii)

            N_bins(j,i) = N_bins(j,i)+chem_1%conc_out_all(ii)
            
            c_p(j,1,i) = c_p(j,1,i)+real(chem_1%clust_out_molec(ii,1),KIND=dp)*chem_1%conc_out_all(ii)*1D-6
            c_p(j,4,i) = c_p(j,4,i)+real(chem_1%clust_out_molec(ii,2),KIND=dp)*chem_1%conc_out_all(ii)*1D-6
            
        end do
          
        DO ii=1,chem_2%nclust_out
            !!! for chemistry 2 H2SO4-DMA

            i = chem_2%ind_out_bin(ii)

            N_bins(j,i) = N_bins(j,i)+chem_2%conc_out_all(ii)
            ! write(*,*) 
            c_p(j,1,i) = c_p(j,1,i)+real(chem_2%clust_out_molec(ii,1),KIND=dp)*chem_2%conc_out_all(ii)*1D-6
            c_p(j,11,i) = c_p(j,11,i)+real(chem_2%clust_out_molec(ii,2),KIND=dp)*chem_2%conc_out_all(ii)*1D-6
            
            ! c_p(j,7,i) = c_p(j,7,i)+init_rel_water*SUM(comp_out_all(ii,:))*c_out_all(ii)*1D-6
        end do
        
        !!! for HIO3-HIO2 chemistry 3
        DO ii=1,chem_3%nclust_out
            !!! for chemistry 2 H2SO4-DMA

            i = chem_3%ind_out_bin(ii)

            N_bins(j,i) = N_bins(j,i)+chem_3%conc_out_all(ii)
            ! write(*,*) 
            c_p(j,10,i) = c_p(j,10,i)+real(chem_3%clust_out_molec(ii,1),KIND=dp)*chem_3%conc_out_all(ii)*1D-6
            ! c_p(j,11,i) = c_p(j,11,i)+real(chem_2%clust_out_molec(ii,2),KIND=dp)*chem_2%conc_out_all(ii)*1D-6
            
            ! c_p(j,7,i) = c_p(j,7,i)+init_rel_water*SUM(comp_out_all(ii,:))*c_out_all(ii)*1D-6
        end do

    end do

    ! write(*,*) real(chem_1%clust_out_molec,KIND=dp)
    write(*,*) 'l2359 nbins ', sum(N_bins(1,1:10))!, sum(chem_1%conc_coag_clust), sum(chem_2%conc_coag_clust) , sum(comp_evap(:,1),dim=1), sum(comp_evap(:,4),dim=1), sum(comp_evap(:,11),dim=1)
    write(*,*) 'l2347 [HIO2] = ', conc(1,ind_HIO2),'',conc(1,ind_HIO3),'',conc(1,ind_H2SO4),'',conc(1,ind_DMA),'',conc(1,ind_NH3)
end if
 
    ! Dry and wet deposition of particles                         !
    !-------------------------------------------------------------!
    v_dep = 0D0
	vs = 0D0 ! Sedimentation velocity
	
    IF (dry_dep_index == 1) THEN 
      dens_p1=dens_p(1,:); d_p1=d_p(1,:)
      CALL dry_dep_particles(landuse_index(tr),d_p1,dens_p1,T2m(tr),P2m(tr),&
     RH2m(tr),fric_veloc(tr),PBLH(tr),SHTFs(tr),windspeed(tr),tr,v_dep,vs,season_index)
    END IF  

    ! Wet deposition velocity (below cloud) in s⁻¹
    v_wet = 0D0
    
	IF (wet_dep_index == 1) THEN
		
     DO j=1,Nz
        d_p1=d_p(j,:)
		IF (SUM(cloudLWC(j:Nz,tr)*dz(j:Nz))/1D3>0.01) THEN ! IF LWC collumn > 0.01 mm
        CALL wet_deposition(d_p1,rain(tr),v_wet1)
        v_wet(j,:)=v_wet1
		ELSE 
		v_wet(j,:)=0D0
		END IF
    END DO    
	
	v_wet=v_wet+v_wet_in_cloud ! Add enhanced particle wet scavening rates for activated cloud droplets in precipitating clouds  
    END IF



    N_bins_old=N_bins
    c_p_old=c_p

    WHERE (N_bins(1,:)>1D-2) N_bins(1,:)=N_bins(1,:)*exp(-v_dep/dz(1)*dt) ! calculates how the particle number concentration changes in each size bin due to dry deposition	
	
    DO j=1,Nz
      exp_vwet = -v_wet(j,:)*dt
      WHERE (exp_vwet <= -1.) exp_vwet = -1.  
      WHERE (N_bins(j,:)>1D-2) N_bins(j,:)=N_bins(j,:)*exp(exp_vwet) ! calculates how the particle number concentration changes in each size bin due to wet deposition
    END DO

    ! Update particle species concentrations after deposition: 
    loss_fraction=N_bins/N_bins_old
    DO i=1,Nz
      DO j=1,nr_bins
      c_p(i,:,j)=c_p_old(i,:,j)*loss_fraction(i,j) ! correct for deposition losses
      END DO
    END DO
	
	 ! Consider particle sedimentation from upper layer to lower altitude layer!
	 c_p_old=c_p
	 N_bins_old=N_bins
     DO ii=1,Nz-1
     N_bins(Nz+1-ii,:)=N_bins(Nz+1-ii,:)*exp(-vs/dz(Nz+1-ii)*dt) ! Sedimentation sink
	 DO j=1,nr_bins
	 c_p(Nz+1-ii,:,j)=c_p(Nz+1-ii,:,j)*exp(-vs(j)/dz(Nz+1-ii)*dt) ! Sedimentation sink
	 END DO
	 END DO
	 DO ii=1,Nz-1
	 N_bins(ii,:)=N_bins(ii,:)+(N_bins_old(ii+1,:)-N_bins(ii+1,:))*dz(ii+1)/dz(ii) ! Sedimentation source from upper layer
	 c_p(ii,:,:)=c_p(ii,:,:)+(c_p_old(ii+1,:,:)-c_p(ii+1,:,:))*dz(ii+1)/dz(ii) ! Sedimentation source from upper layer
	 END DO

	  ! ! Don't allow the particle number concentration in any size bin drop below 0.1 m^-3:
     DO j = 1,nz
        DO i=1,nr_bins
            IF (N_bins(j,i)<1D-1) THEN
            N_bins(j,i)=1D0
            c_p(j,:,i)=c_p_backg(:,i)*1D0;
            vp_dry(j,i) = SUM(c_p(j,index_dry,i)/Na*MX(index_dry)/qX(index_dry)/(N_bins(j,i)*1D-6)) ! m^3
            vp_wet(j,i) = SUM(c_p(j,:,i)/Na*MX/qX/(N_bins(j,i)*1D-6)) ! m^3
            V_bins(j,i) = N_bins(j,i)*vp_wet(j,i)
            dens_p(j,i) = SUM(c_p(j,:,i)*MX*1D6)/Na/V_bins(j,i) ! Total particle density
            dp_dry(j,i) = (vp_dry(j,i)*6D0/pi)**(1D0/3D0) ! Dry particle diameter
            d_p(j,i) = (vp_wet(j,i)*6D0/pi)**(1D0/3D0) ! Wet particle diameter
            END IF
        END DO
     END DO

    DO i=1,Nz
      DO j = 1,nr_bins
      vp_dry(i,j)=SUM(c_p(i,index_dry,j)/Na*MX(index_dry)/qX(index_dry)/(N_bins(i,j)*1D-6)) ! m^3/particle
      vp_wet(i,j)=SUM(c_p(i,:,j)/Na*MX/qX/(N_bins(i,j)*1D-6)) ! m^3/particle
      V_bins(i,j)=N_bins(i,j)*vp_wet(i,j)
      dens_p(i,j)=SUM(c_p(i,:,j)*MX*1D6)/Na/V_bins(i,j) ! Total particle density
      END DO
    END DO
    dp_dry=(vp_dry*6D0/pi)**(1D0/3D0) ! Dry particle diameters
    d_p=(vp_wet*6D0/pi)**(1D0/3D0) ! Dry particle diameters 
	
	    !----------------------------------------------------------------!
     ! Consider primary particle emission of sea salt

    !IF (landuse_index(tr) == 6) THEN ! Trajectory over water
	
	E_sea_salt=E_sea_salt_(:,tr)
	E_snow_salt=E_snow_salt_(:,tr)
	
    !d_p1=dp_dry(1,:)
	!  CALL sea_spray(windspeed(tr),SST(tr),dlogdp,d_p1,d_g,d,E_sea_salt,sea_spray_index)
	! Emissions of NaCl aerosols from sublimating blowing snow!
	!IF (landuse_index(tr)==8) THEN
	!  CALL sea_salt_ice(windspeed(tr),T2m(tr),d_p1,fric_veloc(tr),SHTFs(tr),E_snow_salt)
	!  ELSE
	!  E_snow_salt=0D0
	!END IF  
    	
	
      PN_marine(1,:) =E_snow_salt/dz(1)*dt ! particles/(m³) in each size bin from blowing snow sea spray particles
      PN_marine(2,:) =0.25*E_sea_salt/dz(2)*dt ! particles/(m³) in each size bin from sea spray 
      PN_marine(3,:) =0.2*E_sea_salt/dz(3)*dt ! particles/(m³) in each size bin from sea spray
      PN_marine(4,:) =0.25*E_sea_salt/dz(4)*dt ! particles/(m³) in each size bin from sea spray
      PN_marine(5,:) =0.2*E_sea_salt/dz(5)*dt ! particles/(m³) in each size bin from sea spray
      PN_marine(6,:) =0.1*E_sea_salt/dz(6)*dt ! particles/(m³) in each size bin from sea spray
      
	  V_bins_sea_spray = PN_marine*(dp_dry(1:6,:)**3D0*pi)/6D0  ! particle volume in each size bin

      ! Update particle numbers and chemical composition in each size bin , in the lowest 6 model layers: 
      N_bins(1:6,:) = N_bins(1:6,:)+PN_marine     
     ! Sea spray PM composition Cl:Na =1 	  
      c_p(2:6,3,:) = c_p(2:6,3,:) + MCl/(MCl+MNa)*V_bins_sea_spray(2:6,:)*qCl/MCl*Na*1D-6 ! molec cm^-3
      c_p(2:6,5,:) = c_p(2:6,5,:) + MNa/(MCl+MNa)*V_bins_sea_spray(2:6,:)*qNa/MNa*Na*1D-6  ! molec cm^-3
     ! Blowing snow PM composition Cl:Na:SO4 = 0.9:1.0:0.05 (Some ~10 % depleating in Cl- relative to Na+ Frey et al. Atmos. Chem. Phys., 20, 2549–2578, 2020)	  
      c_p(1,1,:) = c_p(1,1,:) + 0.05*MSO4/(0.9*MCl+MNa+0.05*MSO4)*V_bins_sea_spray(1,:)*qSO4/MSO4*Na*1D-6  ! molec cm^-3
	  c_p(1,3,:) = c_p(1,3,:) + 0.9*MCl/(0.9*MCl+MNa+0.05*MSO4)*V_bins_sea_spray(1,:)*qCl/MCl*Na*1D-6 ! molec cm^-3
      c_p(1,5,:) = c_p(1,5,:) + MNa/(0.9*MCl+MNa+0.05*MSO4)*V_bins_sea_spray(1,:)*qNa/MNa*Na*1D-6  ! molec cm^-3
	  
	  conc(1:6,ind_Brion) = conc(1:6,ind_Brion)+1D-4*SUM(MCl/(MCl+MNa)*V_bins_sea_spray*qX(3)*1D-6/MX(3)*Na,DIM=2) ! Millero et al., 2008 Br:Cl = 0.0015 in reference surface ocean sea water
	  conc(1:6,ind_Clion) = SUM(c_p(1:6,3,:),DIM=2)
    !END IF

    DO i=1,Nz
      DO j = 1,nr_bins
      vp_dry(i,j)=SUM(c_p(i,index_dry,j)/Na*MX(index_dry)/qX(index_dry)/(N_bins(i,j)*1D-6)) ! m^3/particle
      vp_wet(i,j)=SUM(c_p(i,:,j)/Na*MX/qX/(N_bins(i,j)*1D-6)) ! m^3/particle
      V_bins(i,j)=N_bins(i,j)*vp_wet(i,j)
      dens_p(i,j)=SUM(c_p(i,:,j)*MX*1D6)/Na/V_bins(i,j) ! Total particle density
      END DO
    END DO
    dp_dry=(vp_dry*6D0/pi)**(1D0/3D0) ! Dry particle diameters
    d_p=(vp_wet*6D0/pi)**(1D0/3D0) ! Dry particle diameters 
	

     !----------------------------------------------------------------!     

    ! Consider primary particle emission from ship traffic and road  ! 
    ! traffic and marine aerosol                                     !
    !----------------------------------------------------------------!
    !Primary particle emission from ship traffic (scaled with NOx-emission from ship and road from EMEP)
    ! Emission factor NOx = 60g(NOx)/kg(fuel), emission factor PN = 1.5*10**16 particles/kg(fuel), from Beecken et al 2015
    PNship = 0.
	EF_T=-2D13*(Ts(1,tr)-273.15)+6.98D14 ! Ripamonti et al., Tellus B, 2013 T dependent road traffic emission factor Helsinki (PN /km) approx. PN/(g NOx)
    PNroad = EmissionNOx_road(tr)*dt*EF_T/dz(1) ! PN/m^3
	PNnonroad_mach = EmissionNOx_nonroad_mach(tr)*dt*EF_T/(100D0)/dz(1) ! PN/m^3 (100 x lower PN/NOx than road vehicles)
	
    !PNwood = 0.
    IF (landuse_index(tr) == 6) THEN ! only ship emission over water
       PNship =1.5D16/18D0*((EmissionSO2_ship(tr)/(dz(1)*1D2))/Na*1D3*dt/3600.*1D6*MSO2) ! ship (PN/m³), the factor SOx: 18g(SO2)/kg(fuel)
    END IF
    
!	PNbiomass= 5D14*DM_emissions_biomass_burning(tr)*15.3D0/1D3/3600D0/(dz(1)) ! Biomass burning (PN/m^3 s^-1) assuming 5D14 #/(ug PM2.5)  

 NMVOC_corr_factor_wood=MAX(291.1500-Ts(1,tr),1D0)/9.8D0 ! Heating degree day (Simpson et al., 2012), The factor 9.8 is valid for Lund

 N_bins_soot = PN_emission_ship*PNship + PN_emission_S7*PNroad + PN_emission_S9*PNnonroad_mach 
    
	!N_bins_soot = PN_emission_ship*PNship + PN_biomass_burn*PNbiomass*dt + &
    !   (PN_emission_agr_burn*PN_emissions_GAINS(2,tr)+ &
    !   PN_emission_waste*PN_emissions_GAINS(3,tr)+ &
    !   PN_emission_power_plant*PN_emissions_GAINS(5,tr)+ &
    !   PN_emission_industry*PN_emissions_GAINS(6,tr)+ &
    !   PN_emission_flaring*PN_emissions_GAINS(7,tr)+ &
    !   NMVOC_corr_factor_wood*PN_emission_wood*PN_emissions_GAINS(8,tr)+ &
    !   PN_emission_other_burn*PN_emissions_GAINS(9,tr)+ &
    !   PN_emission_coal*PN_emissions_GAINS(10,tr))*dt+& 
	!   PN_emission_S7*PNroad+PN_emission_S9*PNnonroad_mach  ! Traffic emissions scaled according to the NOx emissions from road traffic
	!   ! PN_emission_traffic*PN_emissions_GAINS(4,tr)+ &

    N_bins(1,:) = N_bins(1,:) + N_bins_soot

    V_bins_soot = N_bins_soot*(dp_dry(1,:)**3D0*pi/6D0)               ! particle volume in each size bin
 
    ! Update particle composition after ship emission
    index_core = MAXLOC(dp_dry(1,:) - 70D-9,1,MASK = (dp_dry(1,:) - 70D-9) <=0.) ! first index in dp_dry-vector where dp_dry is bigger than 80 nm
    !index_core2 = MAXLOC(dp_dry(1,:) - 300D-9,1,MASK = (dp_dry(1,:) - 300D-9) <=0.) ! first index in dp_dry-vector where dp_dry is bigger than 300 nm

    ! Update particle composition in core:
    c_p(1,NSPEC_P,1:index_core) = c_p(1,NSPEC_P,1:index_core) + &
          0.7*V_bins_soot(1:index_core)*qX(NSPEC_P)*1D-6/MX(NSPEC_P)*Na
    c_p(1,1,1:index_core) = c_p(1,1,1:index_core) + &
          0.3*V_bins_soot(1:index_core)*qX(1)*1D-6/MX(1)*Na
    ! BC primary aerosol particles 60-300 nm in diameter
    c_p(1,6,index_core+1:nr_bins) = c_p(1,6,index_core+1:nr_bins) + &
          V_bins_soot(index_core+1:nr_bins)*qX(6)*1D-6/MX(6)*Na
          
    DO i=1,Nz
      DO j = 1,nr_bins
      vp_dry(i,j)=SUM(c_p(i,index_dry,j)/Na*MX(index_dry)/qX(index_dry)/(N_bins(i,j)*1D-6)) ! m^3/particle
      vp_wet(i,j)=SUM(c_p(i,:,j)/Na*MX/qX/(N_bins(i,j)*1D-6)) ! m^3/particle
      V_bins(i,j)=N_bins(i,j)*vp_wet(i,j)
      dens_p(i,j)=SUM(c_p(i,:,j)*MX*1D6)/Na/V_bins(i,j) ! Total particle density
      END DO
    END DO
    dp_dry=(vp_dry*6D0/pi)**(1D0/3D0) ! Dry particle diameters
    d_p=(vp_wet*6D0/pi)**(1D0/3D0) ! Dry particle diameters 

    !----------------------------------------------------------!      
    ! Brownian coagulation                                     !
    !----------------------------------------------------------!
    IF (coag_index == 1) THEN
      DO j=1,Nz
        N_bins1=N_bins(j,:); V_bins1=V_bins(j,:); c_p1=c_p(j,:,:); d_p1=d_p(j,:)
        dp_dry1=dp_dry(j,:); dens_p1=dens_p(j,:)
        ! CALL coagulation(N_bins1,V_bins1,c_p1,d_p1,dp_dry1,dt,Ts(j,tr),Ps(j,tr),dens_p1,MX,qX,chem_1,chem_2,n_clustering_vapors,.true.)
        CALL coagulation(N_bins1,V_bins1,c_p1,d_p1,dp_dry1,dt,Ts(j,tr),Ps(j,tr),dens_p1,MX,qX)
        N_bins(j,:)=N_bins1; V_bins(j,:)=V_bins1; c_p(j,:,:)=c_p1; d_p(j,:)=d_p1
        dp_dry(j,:)=dp_dry1; dens_p(j,:)=dens_p1        
      END DO
    END IF    



    !-----------------------------------------------------------!
    ! Internal loop for gas-phase chemistry                     !
    !-----------------------------------------------------------!
    
      ! Update concentration after emissions
	  
      conc(1:kz_can,:) = conc(1:kz_can,:)+E_gases/3600.*dt
	  !conc(1,ind_NH3) = conc(1,ind_NH3)+(sum(EmissionNH3(:,tr))+EmissionNH3_birds(tr))/(dz(1)*1D2)/3600.*dt
	  
	  
      !cNH3(1) = cNH3(1)+((EmissionNH3(tr)+EmissionNH3_birds(tr)+&
      !DM_emissions_biomass_burning(tr)*2.72D0/1D3/18D0/1D4*Na)/(dz(1)*1D2))/3600.*dt    
	  
	  IF (landuse_index(tr)==6 .OR. landuse_index(tr)==8) THEN
      E_222Rn=0.0 ! atom cm^-2 s^-1 (Radon emission 0 over ocean)
      ELSE
      E_222Rn=0.2 ! atom cm^-2 s^-1 Estimated for Northern Europe and Russia (see e.g. Zhang et al., ACP 7817-7838, 2011)  
      END IF
      c222Rn(1) = c222Rn(1)+E_222Rn/(dz(1)*1D2)*dt ! Update surface layer Radon concentration
      k_222Rn=1D0/(3.8*3600D0*24D0) ! Radon decay rate (s^-1)
      c222Rn=c222Rn-k_222Rn*c222Rn*dt ! Consider Radon decay
      q_222Rn=6.8D5*k_222Rn*c222Rn ! Radon induced ionization rate (cm^-3 s^-1)
      q_GCR=1.7D0 ! Estimated Galactic cosmic ray ionization rate (Kirkby et al., 2016 Nature) (cm^-3 s^-1)
      q_ion=q_222Rn+q_GCR
 
	  

	  actinic_flux_t1=actinic_flux_1h(:,(tr_1h-1)*(Nz+1)+1:(tr_1h)*(Nz+1))
	  actinic_flux_t2=actinic_flux_1h(:,(tr_1h)*(Nz+1)+1:(tr_1h+1)*(Nz+1))
	  
	  x2=(REAL(tr)-REAL(tr_1h-1)*60D0)/60D0
	  x1=1D0-x2


    !-----------------------------------------------------------!
    ! Thermodynamics                                            !
    !-----------------------------------------------------------!
    IF (thermodyn_index==1) THEN
    DO ii=1,Nz
        N_bins1=N_bins(ii,:);pH1=pH(ii,:);Kprim_HNO31=Kprim_HNO3(ii,:)
        Kprim_HCl1=Kprim_HCl(ii,:);Hprim_NH31=Hprim_NH3(ii,:);Kprim_NH31=Kprim_NH3(ii,:)
        fHSO41=fHSO4(ii,:);fSO41=fSO4(ii,:);fNO31=fNO3(ii,:);fCl1=fCl(ii,:);mHCO31=mHCO3(ii,:)
        mCO31=mCO3(ii,:);mOH1=mOH(ii,:);mCOO1=mCOO(ii,:)
        W1=W(ii,:);y1=y(ii,:,:)
		Kprim_CH3SO3H1=Kprim_CH3SO3H(ii,:)
		Kprim_HIO31=Kprim_HIO3(ii,:)
	    fCH3SO31=fCH3SO3(ii,:)
		fHIO31=fHIO3(ii,:)
		
        surf_tens = (76.1-0.155*(Ts(ii,tr)-273.15))*1D-3 ! Surface tension of pure water dyn m^-1
        S_Kelvin = exp(2D0*surf_tens*MX(7)/(d_p(ii,:)/2D0*Rg*Ts(ii,tr)*qX(7))) 
        
        !IF (cloudLWC(ii,tr)>1D-2) RH(ii)=99.99
        aw = RH(ii)/100./S_Kelvin  ! Water activity
        c_p1 = c_p(ii,:,:)
    CALL thermodyn_AIOMFAC_inorg(Ts(ii,tr),c_p1,N_bins1,y1,&
        conc(ii,ind_NH3),conc(ii,ind_HNO3),conc(ii,ind_HCl),aw,pCO2,pH1,Kprim_HNO31,&
		Kprim_HCl1,Kprim_CH3SO3H1,Kprim_HIO31,Hprim_NH31,Kprim_NH31,&
        fHSO41,fSO41,fNO31,fCl1,fCH3SO31,fHIO31,mHCO31,mCO31,mOH1,W1)
			    
        c_p(ii,:,:) = c_p1; 
	    N_bins(ii,:)=N_bins1; pH(ii,:)=pH1; Kprim_HNO3(ii,:)=Kprim_HNO31
        Kprim_HCl(ii,:)=Kprim_HCl1; Hprim_NH3(ii,:)=Hprim_NH31; Kprim_NH3(ii,:)=Kprim_NH31
        fHSO4(ii,:)=fHSO41; fSO4(ii,:)=fSO41; fNO3(ii,:)=fNO31; fCl(ii,:)=fCl1; mHCO3(ii,:)=mHCO31
        mCO3(ii,:)=mCO31; mOH(ii,:)=mOH1; mCOO(ii,:)=mCOO1
        W(ii,:)=W1;y(ii,:,:)=y1
		Kprim_CH3SO3H(ii,:)=Kprim_CH3SO3H1
		Kprim_HIO3(ii,:)=Kprim_HIO31
	    fCH3SO3(ii,:)=fCH3SO31
		fHIO3(ii,:)=fHIO31
    
	END DO
    
	END IF

    DO i=1,Nz
      DO j = 1,nr_bins
         vp_dry(i,j)=SUM(c_p(i,index_dry,j)/Na*MX(index_dry)/qX(index_dry)/(N_bins(i,j)*1D-6)) ! m^3
         vp_wet(i,j)=SUM(c_p(i,:,j)/Na*MX/qX/(N_bins(i,j)*1D-6)) ! m^3
         V_bins(i,j)=N_bins(i,j)*vp_wet(i,j)
         dens_p(i,j)=SUM(c_p(i,:,j)*MX*1D6)/Na/V_bins(i,j) ! Total particle density
      END DO
    END DO
   
    dp_dry=(vp_dry*6D0/pi)**(1D0/3D0) ! Dry particle diameters
    d_p=(vp_wet*6D0/pi)**(1D0/3D0) ! Dry particle diameters
 
DO j=1,Nz
      conc1=conc(j,:) 
	  actinic_flux1=actinic_flux_t1(:,j)*x1+actinic_flux_t2(:,j)*x2
	  !tr_3h
	  !actinic_flux_3h
	  
      ! Chemistry
IF (time >= t_start_chem) THEN
	 
          CALL getKVALUES(KVALUES,Ts(j,tr),M(j,tr),O2(j,tr),H2O(j,tr))
          !IF (j<=6) THEN
          !actinic_flux1=actinic_flux1*EM_sun_par(j)
          !ELSE
          !actinic_flux1=actinic_flux1
          !END IF
 
          CALL getJVALUES(JVALUES,Ts(j,tr),actinic_flux1)
		  
		LWC=cloudLWC(j,tr)*1D-3+1D-100 ! kg/m^3
		
		

		
		cw=SUM(c_p(j,7,:))/Na*MH2O ! Inorganic bulk liquid water content	kg/cm^3 = L/cm^3
		IF (LWC>1D-5) THEN
        !cw=EXP(RH(j)-1D2)*LWC*1D-6 ! Inorganic bulk liquid water content	kg/cm^3 = L/cm^3
		cw=LWC*1D-6 ! Inorganic bulk liquid water content	kg/cm^3 = L/cm^3
		!LWC=cw*1D6
        END IF
		!IF (LWC<1D-4) THEN ! kg/m^3
        S_super = 5D-3 ! Assumed maximum supersaturation in the clouds (fraction)
		!ELSEIF (LWC<3D-4) THEN ! kg/m^3
		!ELSE
		!S_super = 3D-3 ! Assumed maximum supersaturation in the clouds (fraction)
		!END IF
		
        surf_tens=(76.1-0.155*(Ts(j,tr)-273.15))*1D-3 ! Surface tension of water (kg s^-2) 
	    ns=(SUM(c_p(j,1:5,:), DIM=1)+SUM(c_p(j,9:NSPEC_P,:),DIM=1))/(N_bins(j,:)*1D-6)/Na ! Moles of soluble material in each single particle inorganic ions + all organic compounds but not soot    
        vp0=(c_p(j,6,:)/Na*MEC/qEC)/(N_bins(j,:)*1D-6) ! Equivalent volume of insoluble core (considered to be soot only)
        dp0=(vp0*6D0/pi)**(1D0/3D0) ! Equivalent diameter of insoluble core (considered to be soot only)
     
        A_Kohler=4*MH2O*surf_tens/(Rg*Ts(j,tr)*qH2O)
        B_Kohler=6*ns*MH2O/(pi*qH2O)
 
        ! Calculate the critical supersaturation for droplet activation according to Kokkola et al., 2008 
        !alpha1=12*(81D0*dp0**6D0+12D0*dp0**3D0*(3D0*B_Kohler/A_Kohler)**(3D0/2D0))**(1D0/2D0)
        !alpha2=(108D0*dp0**3D0+8*(3D0*B_Kohler/A_Kohler)**(3D0/2D0)+alpha1)**(1D0/3D0)
        !S_c= EXP( A_Kohler/(alpha2/6D0 + 2D0/3D0*(3D0*B_Kohler/A_Kohler)*1D0/alpha2+&
        !1D0/3D0*(3D0*B_Kohler/A_Kohler)**(1D0/2D0))+&
        !B_Kohler/((alpha2/6D0 + 2D0/3D0*(3D0*B_Kohler/A_Kohler)*1D0/alpha2+&
        !1D0/3D0*(3D0*B_Kohler/A_Kohler)**(1D0/2D0))**3D0 - dp0**3D0))-1D0
		
		b_ko=sqrt(3D0*B_Kohler/A_Kohler)
        d_ko=dp0**3D0
        lambda_ko1=(2D0*(b_ko/3D0)**3D0+d_ko)+sqrt((4D0*(b_ko/3D0)**3D0)*d_ko+d_ko**2D0)
        lambda_ko2=(2D0*(b_ko/3D0)**3D0+d_ko)-sqrt((4D0*(b_ko/3D0)**3D0)*d_ko+d_ko**2D0)
        S_c= A_Kohler*(2D0*(b_ko/3D0)**2D0+&
        (b_ko/3D0)*((lambda_ko1/2D0)**(1D0/3D0)+(lambda_ko2/2D0)**(1D0/3D0))+&
        ((lambda_ko1/2D0)**(2D0/3D0)+(lambda_ko2/2D0)**(2D0/3D0)))/&
        (9D0*(b_ko/3D0)**3D0+&
        (6D0*(b_ko/3D0)**2D0)*((lambda_ko1/2D0)**(1D0/3D0)+(lambda_ko2/2D0)**(1D0/3D0))+&
        3D0*(b_ko/3D0)*((lambda_ko1/2D0)**(2D0/3D0)+(lambda_ko2/2D0)**(2D0/3D0)) +d_ko)
		
        CCN=SUM(PACK(N_bins(j,:), S_c < S_super)) ! Number of CCN (activated droplets at 0.1 % supersaturation)
        fCCN=N_bins(j,:)/CCN
        
        m_drops=(cw*1D6)/CCN ! Average single droplet mass/volume of cloud droplets (kg or dm^3)
        dp_drops=(m_drops/qH2O*6D0/pi)**(1D0/3D0) ! Diameter of cloud droplets
        !IF (LWC>1D-5 .AND. dp_drops<MAXVAL(d_p(j,:))) THEN
		!write(*,*) 'dp_drops<max(d_p)'
		!write(*,*) MAXVAL(d_p(j,:))
		!dp_drops=MAXVAL(d_p(j,:))*1.2
		!END IF
       dyn_visc=1.8D-5*(Ts(j,tr)/298.)**0.85                ! dynamic viscosity of air
       l_gas=2D0*dyn_visc/(Ps(j,tr)*SQRT(8D0*Mair/(pi*Rg*Ts(j,tr)))) ! Gas mean free path in air (m)
       Cunningh=1D0+(2D0*l_gas/d_p(j,:))*(1.257+0.4*EXP(-1.1/(2D0*l_gas/d_p(j,:))))  ! Cunninghams correction factor
       Diff_p=Cunningh*kb*Ts(j,tr)/(3D0*pi*dyn_visc*d_p(j,:))                          ! Diffusivitys for the different particle sizes m^2/s
       m_p=(dens_p(j,:)*pi*d_p(j,:)**3D0)/6D0  ! mass of particles 
       speed_p=SQRT(8D0*kb*Ts(j,tr)/(pi*m_p)) ! speed of particles
 
       c_speed_aq=SQRT(8D0*kb*Ts(j,tr)/(pi*(MX_aq/Na)))     ! Thermal speed of molecules
      DO i=1,51
         gasmeanfp_aq=3D0*(D_aq(i)+Diff_p)/SQRT(c_speed_aq(i)**2D0+speed_p**2D0) 
         Kn_aq=2D0*gasmeanfp_aq/(d_p(j,:)+dX_aq(i)) ! Knudsen number
         f_cor_aq=(0.75*alpha_aq(i)*(1D0+Kn_aq))/(Kn_aq**2D0+Kn_aq+0.283*Kn_aq*alpha_aq(i)+0.75*alpha_aq(i)) ! Fuchs-Sutugin correction factor for transit and kinetic regime a=accomodation coefficient
         Deff_aq=(D_aq(i)+Diff_p)*f_cor_aq                   ! m^2/s

      IF (LWC>1D-5) THEN 
	   where (S_c<S_super .AND. d_p(j,:)<dp_drops) d_p(j,:)=dp_drops
	    CS_AQ(i)=SUM(N_bins(j,:)*2D0*pi*(d_p(j,:)+dX_aq(i))*Deff_aq)  ! mass transfer coefficient s^-1
       ELSE
        CS_AQ(i)=SUM(N_bins(j,:)*2D0*pi*(d_p(j,:)+dX_aq(i))*Deff_aq)  ! mass transfer coefficient s^-1
      END IF
      END DO
      
	  
      DO i=1,nr_bins
      vp_dry(j,i) = SUM(c_p(j,index_dry,i)/Na*MX(index_dry)/qX(index_dry)/(N_bins(j,i)*1D-6)) ! m^3
      vp_wet(j,i) = SUM(c_p(j,:,i)/Na*MX/qX/(N_bins(j,i)*1D-6)) ! m^3
      V_bins(j,i) = N_bins(j,i)*(d_p(j,i)**3D0)*pi/6D0
      V_bins_dry(i)=N_bins(j,i)*vp_dry(j,i)
      END DO
	  c_p_old=c_p
      WHERE (S_c<S_super) c_p(j,7,:)=(V_bins(j,:)-V_bins_dry)*1D-6*qH2O/MH2O*Na 
      WHERE (c_p(j,7,:)>1D3*c_p_old(j,7,:)) c_p(j,7,:)=c_p_old(j,7,:)*1D3  ! Updated 2022-01-04 by Pontus
	  cw=SUM(c_p(j,7,:))/Na*MH2O ! Inorganic bulk liquid water content	kg/cm^3 = L/cm^3
	  IF (LWC>1D-5 .AND. cw<LWC*1D-6 ) THEN
	  c_p(j,7,:)=c_p(j,7,:)*(LWC*1D-6)/cw
	  cw=SUM(c_p(j,7,:))/Na*MH2O
	  END IF

!	  IF (LWC>1D-5) THEN
!	  Hion=1D-2
!	  ELSE
!	  Hion=SUM(c_p(j,7,:)/Na*MH2O*10**(-pH(j,:)))/SUM(c_p(j,7,:)/Na*MH2O)
!	  IF (Hion>0.1D0) THEN 
!	  Hion=0.1D0 ! Limit the acidity
!      END IF
!	  END IF

    c_p_bulk=SUM(c_p(j,:,:),DIM=2) ! bulk species content	molec/cm^3
	pH_bulk1=pH_bulk(j)
	yi_bulk1=yi_bulk(j,:)
	CALL thermodyn_AIOMFAC_inorg_bulk(Ts(ii,tr),c_p_bulk,pCO2,pH_bulk1,yi_bulk1)
	pH_bulk(j)=pH_bulk1
	yi_bulk(j,:)=yi_bulk1
	Hion=10**(-pH_bulk1)
	  OHion=(1D-14/Hion) ! Bulk [OH-] mol/kg	
	  CS_AQ(46:47)=1D-100 ! The condensation of NH3 is treated by the aerosol dynamics module
	  CS_AQ(49:51)=1D-100 ! The condensation of HNO3, H2SO4 and MSA is treated by the aerosol dynamics  module
      CS_AQ(5)=1D-100 ! The condensation of HCl is treated by the aerosol dynamics module
	  CS_AQ(25)=1D-100 ! The condensation of HIO3 is treated by the aerosol dynamics module

     !fHCl=conc1(ind_HClaq)/(conc1(ind_Clion)+conc1(ind_HClaq)+1D-100)
     
	  conc_old=conc1
      conc1(ind_Clion) = conc_old(ind_Clion)/(conc_old(ind_Clion)+conc_old(ind_HClaq)+1D-100)*SUM(c_p(j,3,:)) ! Total Cl- particle conc
      conc1(ind_HClaq) = conc_old(ind_HClaq)/(conc_old(ind_Clion)+conc_old(ind_HClaq)+1D-100)*SUM(c_p(j,3,:)) ! Total HCl particle conc
      conc1(ind_NO3ion) = conc_old(ind_NO3ion)/(conc_old(ind_NO3ion)+conc_old(ind_HNO3aq)+1D-100)*SUM(c_p(j,2,:)) ! Total NO3- particle conc
	  conc1(ind_HNO3aq) = conc_old(ind_HNO3aq)/(conc_old(ind_NO3ion)+conc_old(ind_HNO3aq)+1D-100)*SUM(c_p(j,2,:)) ! Total NO3- particle conc
	  conc1(ind_NH4ion) = conc_old(ind_NH4ion)/(conc_old(ind_NH4ion)+conc_old(ind_NH3aq)+1D-100)*SUM(c_p(j,4,:)) ! Total N(III) particle conc
	  conc1(ind_NH3aq) = conc_old(ind_NH3aq)/(conc_old(ind_NH4ion)+conc_old(ind_NH3aq)+1D-100)*SUM(c_p(j,4,:)) ! Total NO3- particle conc
      conc1(ind_DMAion) = conc_old(ind_DMAion)/(conc_old(ind_DMAion)+conc_old(ind_DMAaq)+1D-100)*SUM(c_p(j,11,:)) ! Total DMA+ particle conc
	  conc1(ind_DMAaq) = conc_old(ind_DMAaq)/(conc_old(ind_DMAion)+conc_old(ind_DMAaq)+1D-100)*SUM(c_p(j,11,:)) ! Total DMA particle conc
	  
	  conc1(ind_SO42ion) = conc_old(ind_SO42ion)/(conc_old(ind_HSO4ion)+conc_old(ind_SO42ion)+1D-100)*SUM(c_p(j,1,:)) ! Total SO42- particle conc
      conc1(ind_HSO4ion) = conc_old(ind_HSO4ion)/(conc_old(ind_HSO4ion)+conc_old(ind_SO42ion)+1D-100)*SUM(c_p(j,1,:)) ! Total HSO4- particle conc
	  conc1(ind_CH3SO3ion) = conc_old(ind_CH3SO3ion)/(conc_old(ind_CH3SO3ion)+conc_old(ind_CH3SO3Haq)+1D-100)*SUM(c_p(j,9,:)) ! Total MSA particle conc
      conc1(ind_CH3SO3Haq) = conc_old(ind_CH3SO3Haq)/(conc_old(ind_CH3SO3ion)+conc_old(ind_CH3SO3Haq)+1D-100)*SUM(c_p(j,9,:)) ! Total MSA particle conc
      conc1(ind_IO3ion) = conc_old(ind_IO3ion)/(conc_old(ind_IO3ion)+conc_old(ind_HIO3aq)+1D-100)*SUM(c_p(j,10,:)) ! Total IO3- particle conc
      conc1(ind_HIO3aq) = conc_old(ind_HIO3aq)/(conc_old(ind_IO3ion)+conc_old(ind_HIO3aq)+1D-100)*SUM(c_p(j,10,:)) ! Total HIO3 particle conc
     
     conc_old=conc1
     
	 cond_sinkH2SO4=0D0
	 cond_sinkHNO3=0D0
	 !where (conc1<0D0) conc1=0D0 ! For safety!

	 CALL KPP_Proceed(conc1,time,time+dt,Ts(j,tr),RH(j),Ps(j,tr),O2(j,tr),N2(j,tr),M(j,tr),H2O(j,tr),&
	 cond_sinkH2SO4,cond_sinkHNO3,CS_AQ,cw,Hion,OHion,JVALUES,KVALUES)
	 !where (conc1<0D0) conc1=0D0 ! For safety!
	 
     ! Update S(VI), MSA and HIO3 particle-phase content after aqueous phase chemistry:
     dconc=conc1-conc_old

	 IF (dconc(ind_SO42ion)+dconc(ind_HSO4ion)>0D0) THEN
	 c_p(j,1,:)=c_p(j,1,:)+(dconc(ind_SO42ion)+dconc(ind_HSO4ion))*c_p(j,7,:)/SUM(c_p(j,7,:))
	 ELSE
	 c_p(j,1,:)=c_p(j,1,:)+(dconc(ind_SO42ion)+dconc(ind_HSO4ion))*c_p(j,1,:)/SUM(c_p(j,1,:))
	 END IF
	 IF (dconc(ind_CH3SO3ion)+dconc(ind_CH3SO3Haq)>0D0) THEN
     c_p(j,9,:)=c_p(j,9,:)+(dconc(ind_CH3SO3ion)+dconc(ind_CH3SO3Haq))*c_p(j,7,:)/SUM(c_p(j,7,:))
	 ELSE 
	 c_p(j,9,:)=c_p(j,9,:)+(dconc(ind_CH3SO3ion)+dconc(ind_CH3SO3Haq))*c_p(j,9,:)/SUM(c_p(j,9,:))
	 END IF
	 !IF (dconc(ind_IO3ion)+dconc(ind_HIO3aq)>0D0) THEN
     !c_p(j,10,:)=c_p(j,10,:)+(dconc(ind_IO3ion)+dconc(ind_HIO3aq))*c_p(j,7,:)/SUM(c_p(j,7,:))
	 !ELSE 
	 !c_p(j,10,:)=c_p(j,10,:)+(dconc(ind_IO3ion)+dconc(ind_HIO3aq))*c_p(j,10,:)/SUM(c_p(j,10,:))
	 !END IF
	 
     
    DO i = 1,nr_bins
        vp_dry(j,i)=SUM(c_p(j,index_dry,i)/Na*MX(index_dry)/qX(index_dry)/(N_bins(j,i)*1D-6)) ! m^3     
    END DO
     vp_fixed(1:nr_bins)=(dp_dry(j,:)**3D0)*pi/6D0 ! dry single particle particle volume
     N_bins_fixed=0D0
     c_p_fixed=0D0
     dp_max = dp_dry(j,nr_bins)*dp_dry(j,nr_bins)/dp_dry(j,nr_bins-1)
     vp_fixed(nr_bins+1)=(pi*dp_max**3D0)/6D0

     ! ! Full-stationary size distribution: 
     DO i = 1,nr_bins
        a = MINLOC(vp_fixed-vp_dry(j,i),1,mask = (vp_fixed-vp_dry(j,i)) >= 0D0)
        IF (a > i) THEN ! growth
        r1 = (vp_fixed(a)-vp_dry(j,i))/(vp_fixed(a)-vp_fixed(a-1)) ! Fraction of particles in size bin a-1
        r2 = 1D0-r1 ! Fraction of particles in next size bin (a)
        IF (a > nr_bins+1) THEN
        ELSE IF (a == nr_bins+1) THEN
        N_bins_fixed(a-1) = N_bins_fixed(a-1)+r1*N_bins(j,i)
        c_p_fixed(:,a-1) = c_p_fixed(:,a-1)+vp_fixed(a-1)/vp_dry(j,i)*r1*c_p(j,:,i)
        ELSE
        N_bins_fixed(a-1) = N_bins_fixed(a-1)+r1*N_bins(j,i)
        N_bins_fixed(a) = N_bins_fixed(a)+r2*N_bins(j,i)
        c_p_fixed(:,a-1) = c_p_fixed(:,a-1)+vp_fixed(a-1)/vp_dry(j,i)*r1*c_p(j,:,i)
        c_p_fixed(:,a) = c_p_fixed(:,a)+vp_fixed(a)/vp_dry(j,i)*r2*c_p(j,:,i)
        END IF
        ELSE ! evaporation
        IF (a == 0) THEN           
        ELSEIF (a == 1) THEN
        r1 = (vp_dry(j,i))/(vp_fixed(1)) ! Fraction of particles in size bin a
        N_bins_fixed(1) = N_bins_fixed(1)+r1*N_bins(j,i)
        c_p_fixed(:,1) = c_p_fixed(:,1)+vp_fixed(1)/vp_dry(j,i)*r1*c_p(j,:,i)
        ELSE
        r1 = (vp_dry(j,i)-vp_fixed(a-1))/(vp_fixed(a)-vp_fixed(a-1)) ! Fraction of particles in size bin i
        r2 = 1D0-r1 ! Fraction of particles in previous size bin (i-1)
        N_bins_fixed(a) = N_bins_fixed(a)+r1*N_bins(j,i)
        N_bins_fixed(a-1) = N_bins_fixed(a-1)+r2*N_bins(j,i)
        c_p_fixed(:,a) = c_p_fixed(:,a)+vp_fixed(a)/vp_dry(j,i)*r1*c_p(j,:,i)
        c_p_fixed(:,a-1) = c_p_fixed(:,a-1)+vp_fixed(a-1)/vp_dry(j,i)*r2*c_p(j,:,i)
       END IF
       END IF
     END DO

      c_p(j,:,:)=c_p_fixed
      N_bins(j,:)=N_bins_fixed

 ! Don't allow the particle number concentration in any size bin to drop below 0.1 m^-3:
      DO i = 1,nr_bins
        IF (N_bins(j,i)<=1D-1) THEN ! This can occur with the full-stationary method if some particles grow to more than the next particle size bin within one condensation time step
          N_bins(j,i)=1D0
          c_p(j,:,i)=c_p_backg(:,i)*1D0;
        END IF       
        vp_dry(j,i)=SUM(c_p(j,index_dry,i)/Na*MX(index_dry)/qX(index_dry)/(N_bins(j,i)*1D-6)) ! m^3     
        vp_wet(j,i)=SUM(c_p(j,:,i)/Na*MX/qX/(N_bins(j,i)*1D-6)) ! m^3
        V_bins(j,i)=N_bins(j,i)*vp_wet(j,i)
        dens_p(j,i)=SUM(c_p(j,:,i)*MX*1D6)/Na/V_bins(j,i) ! Total particle density
      END DO
      dp_dry(j,:)=(vp_dry(j,:)*6D0/pi)**(1D0/3D0) ! Dry particle diameters
      d_p(j,:)=(vp_wet(j,:)*6D0/pi)**(1D0/3D0) ! Dry particle diameters 


! End section with gas and aqueous phase chemistry
	     		 
END IF
conc(j,:)=conc1

! Consider enhanced in cloud scavening of aerosol particles during precipitation events (used in next model time step):
v_wet_in_cloud(j,:)=0D0

IF (rain(tr)>1D-2 .AND. LWC>1D-5) THEN
Win=1D6 ! Table S20 Simpson et al. (2012)
hs=1D3 ! m
where (S_c<S_super) v_wet_in_cloud(j,:)=Win*(rain(tr)/3.6D3)/(hs*qH2O)
ELSE
v_wet_in_cloud(j,:)=0D0
END IF
       
END DO  
         
    conc(:,ind_CH4) = 1900d-9*M(:,1) ! CH4

    IF (thermodyn_index==1) THEN
    DO ii=1,Nz
        N_bins1=N_bins(ii,:);pH1=pH(ii,:);Kprim_HNO31=Kprim_HNO3(ii,:)
        Kprim_HCl1=Kprim_HCl(ii,:);Hprim_NH31=Hprim_NH3(ii,:);Kprim_NH31=Kprim_NH3(ii,:)
        fHSO41=fHSO4(ii,:);fSO41=fSO4(ii,:);fNO31=fNO3(ii,:);fCl1=fCl(ii,:);mHCO31=mHCO3(ii,:)
        mCO31=mCO3(ii,:);mOH1=mOH(ii,:);mCOO1=mCOO(ii,:)
        W1=W(ii,:);y1=y(ii,:,:)
		Kprim_CH3SO3H1=Kprim_CH3SO3H(ii,:)
		Kprim_HIO31=Kprim_HIO3(ii,:)
	    fCH3SO31=fCH3SO3(ii,:)
		fHIO31=fHIO3(ii,:)
        surf_tens = (76.1-0.155*(Ts(ii,tr)-273.15))*1D-3 ! Surface tension of pure water dyn m^-1
        S_Kelvin = exp(2D0*surf_tens*MX(7)/(d_p(ii,:)/2D0*Rg*Ts(ii,tr)*qX(7))) 
        
        !IF (cloudLWC(ii,tr)>1D-2) RH(ii)=99.99
        aw = RH(ii)/100./S_Kelvin  ! Water activity
        c_p1 = c_p(ii,:,:)
    CALL thermodyn_AIOMFAC_inorg(Ts(ii,tr),c_p1,N_bins1,y1,&
          conc(ii,ind_NH3),conc(ii,ind_HNO3),conc(ii,ind_HCL),aw,pCO2,pH1,Kprim_HNO31,&
		  Kprim_HCl1,Kprim_CH3SO3H1,Kprim_HIO31,Hprim_NH31,Kprim_NH31,&
          fHSO41,fSO41,fNO31,fCl1,fCH3SO31,fHIO31,mHCO31,mCO31,mOH1,W1)		    
        c_p(ii,:,:) = c_p1; 
	    N_bins(ii,:)=N_bins1; pH(ii,:)=pH1; Kprim_HNO3(ii,:)=Kprim_HNO31
        Kprim_HCl(ii,:)=Kprim_HCl1; Hprim_NH3(ii,:)=Hprim_NH31; Kprim_NH3(ii,:)=Kprim_NH31
        fHSO4(ii,:)=fHSO41; fSO4(ii,:)=fSO41; fNO3(ii,:)=fNO31; fCl(ii,:)=fCl1; mHCO3(ii,:)=mHCO31
        mCO3(ii,:)=mCO31; mOH(ii,:)=mOH1; mCOO(ii,:)=mCOO1
        W(ii,:)=W1;y(ii,:,:)=y1
		Kprim_CH3SO3H(ii,:)=Kprim_CH3SO3H1
		Kprim_HIO3(ii,:)=Kprim_HIO31
	    fCH3SO3(ii,:)=fCH3SO31
		fHIO3(ii,:)=fHIO31
    
	END DO
    
	END IF

    DO i=1,Nz
      DO j = 1,nr_bins
         vp_dry(i,j)=SUM(c_p(i,index_dry,j)/Na*MX(index_dry)/qX(index_dry)/(N_bins(i,j)*1D-6)) ! m^3
         vp_wet(i,j)=SUM(c_p(i,:,j)/Na*MX/qX/(N_bins(i,j)*1D-6)) ! m^3
         V_bins(i,j)=N_bins(i,j)*vp_wet(i,j)
         dens_p(i,j)=SUM(c_p(i,:,j)*MX*1D6)/Na/V_bins(i,j) ! Total particle density
      END DO
    END DO
   
    dp_dry=(vp_dry*6D0/pi)**(1D0/3D0) ! Dry particle diameters
    d_p=(vp_wet*6D0/pi)**(1D0/3D0) ! Dry particle diameters
 
    !---------------------------------------------------------------!
    ! Condensation                                                  !
    !---------------------------------------------------------------!
	
	
    IF (condensation_index == 1) THEN
      yorg = 1. ! set this to 1 for now
    DO j=1,Nz
	psat_org=10.**(A_Nannoolal-B_Nannoolal/T2m(tr))*101325D0 ! Vapor pressures condensable organic comp (atm).
    psat_org(635:742)=psat_org(635:742)
    psat_org(635:645)=1D-100 ! RO2
    psat_org(756:761)=1D-100 ! RO2
	
     Dorg=1D-7*Ts(j,tr)**1.75D0*SQRT(1/(Mair*1D3)+1/mcm_prop(1,:))/&
          (Ps(j,tr)/1.01325D5*(VolX**(1D0/3D0)+20.1**(1D0/3D0))**2D0) ! Fuller's method /(Tang et al., ACP 2015) ! Diffusivity organic compounds (m^2 s^-1)
                  
          corg(1:NCOND)=conc(j,index_cond)  ! Vector with concentrations of condensable organic compounds
          N_bins1=N_bins(j,:); V_bins1=V_bins(j,:); d_p1=d_p(j,:); dp_dry1=dp_dry(j,:)
          dens_p1=dens_p(j,:); c_p1=c_p(j,:,:); W1=W(j,:); Kprim_HNO31=Kprim_HNO3(j,:)
          Kprim_HCl1=Kprim_HCl(j,:); Hprim_NH31=Hprim_NH3(j,:); Kprim_NH31=Kprim_NH3(j,:)
          fHSO41=fHSO4(j,:); fSO41=fSO4(j,:); fNO31=fNO3(j,:); fCl1=fCl(j,:); mHCO31=mHCO3(j,:)
          mCO31=mCO3(j,:); mOH1=mOH(j,:)
          dimer_C1=dimer_C(j,:); dimer_O1=dimer_O(j,:);dimer_N1=dimer_N(j,:);dimer_H1=dimer_H(j,:)
		  Kprim_CH3SO3H1=Kprim_CH3SO3H(j,:)
          Kprim_HIO31=Kprim_HIO3(j,:)			  
		  fCH3SO31=fCH3SO3(j,:)
		  fHIO31=fHIO3(j,:)
		  CALL condensation(N_bins1,V_bins1,Ps(j,tr),Ts(j,tr),RH(j),&
		       d_p1,dp_dry1,aX,MX,qX,molec_radius,dens_p1,corg,conc(j,ind_H2SO4),conc(j,ind_HNO3),&
               conc(j,ind_HCl),conc(j,ind_MSA),conc(j,ind_HIO3),conc(j,ind_NH3),conc(j,ind_DMA),conc(j,ind_HIO2),&
			   c_p1,W1,Kprim_HNO31,Kprim_HCl1,Kprim_CH3SO3H1,Kprim_HIO31,Hprim_NH31,Kprim_NH31,&
			   fHSO41,fSO41,fNO31,fCl1,fCH3SO31,fHIO31,mHCO31,mCO31,mOH1,psat_org,&
			   yorg,Dorg,CS_H2SO4(j),CS_air(j),c_p_backg,dt, n_evap_not_inuse(j),comp_evap_via_condensation(j,:), &
               clust_evap,.True.)

            !    write(*,*) 'outside module nconc_evap1 and 2' , chem_1%Nconc_evap(j),chem_2%Nconc_evap(j) , n_evap1(j),n_evap2(j)  
			   
          conc(j,index_cond)=corg(1:NCOND) 
          
          N_bins(j,:)=N_bins1; V_bins(j,:)=V_bins1; d_p(j,:)=d_p1; dp_dry(j,:)=dp_dry1
          dens_p(j,:)=dens_p1; c_p(j,:,:)=c_p1; dimer_C(j,:)=dimer_C1; dimer_O(j,:)=dimer_O1;
          dimer_N(j,:)=dimer_N1; dimer_H(j,:)=dimer_H1;
        !END DO		
      END DO
    !   chem_1%Nconc_evap=n_evap1
    !   chem_2%Nconc_evap=n_evap2
    !   comp_evap=comp_evap_via_condensation


    !   write(*,*) 'l2843 Nconc_evap', sum(chem_1%Nconc_evap), sum(chem_2%Nconc_evap) , sum(comp_evap(:,1),dim=1), &
    !   sum(comp_evap(:,4),dim=1), sum(comp_evap(:,11),dim=1)
    END IF
 

    !-------------------------------------------------------!
    ! Write to output-files                      !
    !-------------------------------------------------------!
	
	RH_old=RH
	
    !CALL write_output(time,conc,index_cond,N_bins,V_bins,d_p,dp_dry,c_p,v_dep,v_wet1,dens_p,&
    !index_st,mcm_prop(1,:)*1D-3,psat_org,E_gases,Kz,vd_gas,tr,&
    !PN_marine,Jnucl,CS_H2SO4,RH)
	IF (tr==1 .OR. tr==index_st) THEN
    CALL write_output(time,conc,index_cond,N_bins,V_bins,d_p,dp_dry,c_p,v_dep,v_wet1,dens_p,&
    index_st,mcm_prop(1,:)*1D-3,psat_org,E_gases,Kz,vd_gas,tr,&
    PN_marine,Jnucl_N, Jnucl_D,CS_H2SO4,RH,Ts(:,tr),pH)
    END IF
	

    time = time + dt
!    tr = tr + 1
   tr = CEILING(time/60D0)
   tr_1h = CEILING(time/(3600D0)) 
   !write(*,*) conc(1,ind_O3), conc(1,ind_OH)
   !write(*,*) sum(N_bins(1,:))*1D-6
   !write(*,*) sum(V_bins(1,:))*1D12
   !write(*,*) time

    write(*,*) 'End loop Time: ', ' ', tr, ' ', time 
    END DO loop

    !CLOSE(200)
    !CLOSE(201)
    !CLOSE(202)
    !CLOSE(203)
    CLOSE(204)
    !CLOSE(205)
    !CLOSE(206)
    !CLOSE(207)
    !CLOSE(208)
    !CLOSE(209)
    !CLOSE(210)
    !CLOSE(211)
    CLOSE(212)
    CLOSE(213)
    !CLOSE(214)
    CLOSE(215)
    CLOSE(216)
    CLOSE(217)
    CLOSE(218)
    !CLOSE(219)
    CLOSE(220)
    !CLOSE(221)
    CLOSE(222)
    CLOSE(223)
    !CLOSE(224)
    !CLOSE(225)
    !CLOSE(226)
    !CLOSE(227)
    !CLOSE(228)
    !CLOSE(229)
    CLOSE(230)
    CLOSE(231)
    CLOSE(300)
    CLOSE(232)
    CLOSE(233)
    CLOSE(234)
    !CLOSE(235)
    
    CALL cpu_time(finish_time)
    PRINT '("Simulation time = ",f6.3," hours")', (finish_time-start_time)/3600.

CONTAINS


!------------------------------------------------------------!
! Vertical diffusion                                         !
!------------------------------------------------------------!

SUBROUTINE diffusion1D_variable_z(Kz,c,k_uper,dt)
! Calculate the atmospheric diffusion equation for any 
! gas or particle without emissions and deposition: 
! For each time the atmospheric diffusion equation is sloved, it is assumed
! that each grid box consists of one unique trase substance with
! concentration equal to 1. 
REAL(dp), DIMENSION(Nz+1), INTENT(in) :: Kz
REAL(dp), INTENT(in) :: k_uper,dt
REAL(dp), DIMENSION(Nz), INTENT(inout) :: c
REAL(dp), DIMENSION(Nz) :: c_old
REAL(dp) :: dt_diff
INTEGER :: N

! run 1D dispersion
 c_diff=0D0
!DO i=1,Nz
 !  concentration in each grid cell:
! c=0D0
! c(i)=1D0 ! include trace species in one grid box

 dt_diff=MINVAL((/ MINVAL(0.1D0*(dz**2D0)/(Kz(2:Nz+1)+1D-10)),dt /))

 N =CEILING(dt/dt_diff)
 dt_diff=dt/REAL(N)
 
 c_old=c
 DO j=1,N
    c(1)=c_old(1)+dt_diff/dz(1)*(Kz(2)/dz2(1)*(c_old(2)-c_old(1)))
    c(2:Nz-1)=c_old(2:Nz-1)+dt_diff/dz(2:Nz-1)*(Kz(3:Nz)/dz2(2:Nz-1)*(c_old(3:Nz)-c_old(2:Nz-1))-&
    Kz(2:Nz-1)/dz2(1:Nz-2)*(c_old(2:Nz-1)-c_old(1:Nz-2)))
    c(Nz)=c_old(Nz)+dt_diff/dz(Nz)*(Kz(Nz+1)*k_uper-Kz(Nz)/dz2(Nz-1)*(c_old(Nz)-c_old(Nz-1)))
    c_old=c
 END DO

where (c<0D0) c=0D0
!c_diff(i,:) = c
!END DO


END SUBROUTINE diffusion1D_variable_z

!SUBROUTINE diffusion1D(Kz,c_diff,k_uper_BC)

!! Calculate the atmospheric diffusion equation for any 
!! gas or particle without emissions and deposition: 
!! For each time the atmospheric diffusion equation is sloved, it is assumed
!! that each grid box consists of one unique trase substance with
!! concentration equal to 1. 
!REAL(dp), DIMENSION(Nz+1), INTENT(in) :: Kz
!REAL(dp), DIMENSION(Nz,Nz) :: T1,dT1
!REAL(dp), DIMENSION(Nz,Nz), INTENT(out) :: c_diff
!REAL(dp), DIMENSION(Nz) :: c, R, dcdt
!REAL(dp) :: dt_diff
!REAL(dp), INTENT(in) ::  k_uper_BC ! slope of conc. gradient as uper BC (dc/dz=k_uper_BC*c)
!INTEGER :: i,N

!! tree-diagonal toeplitz matrix for the vertical direction 
!T1 = 0D0
!DO i = 1, Nz
!   T1(i,i) =-Kz(i)-Kz(i+1)
!!   I_mat(i,i)=1D0
!   IF (i>1) THEN
!   T1(i,i-1) =Kz(i)
!   END IF
!   IF (i<Nz) THEN
!   T1(i,i+1) =Kz(i+1)
!   END IF
!END DO

!T1(Nz,Nz)=T1(Nz,Nz)+Kz(Nz+1) ! or +Kz(N) ?
!T1(1,1)=T1(1,1)+Kz(1)


!! run 1D dispersion
! c_diff=0D0
!DO i=1,Nz
! !  concentration in each grid cell:
! c=0D0
! c(i)=1D0 ! include trace species in one grid box
! dT1=(1D0/(dz**2D0))*T1 ! semi-discretization

! ! rest terms:
! R=0D0
! dt_diff=MINVAL((/ MINVAL(0.25D0*(dz**2D0)/(Kz+1D-10)),dt /))
! !dt_diff=MINVAL((/ 0.25D0*(dz**2D0)/(Kz(i)+1D-10),dt /))

! N =CEILING(dt/dt_diff)
! dt_diff=dt/REAL(N)
! DO j=1,N
!  R(Nz)=k_uper_BC*Kz(Nz+1)*c(Nz)/dz ! include uper BC (dc/dz=k)
!  dcdt=MATMUL(dT1,c)+R
!  c=c+dcdt*dt_diff !
! END DO

! ! Solves the semidiscretized atmospheric diffusion equation one time step forward in time:
! !c_diff(i,:)=c+dcdt*dt ! Forward Euler explicit method
! c_diff(i,:) = c
!END DO

!END SUBROUTINE diffusion1D

!----------------------------------------------------!
! In canopy wind profile                             !
!----------------------------------------------------!
SUBROUTINE in_canopy_wind_profile(temp,press,u,landuse_index,Wind,Zin)
    
      REAL(dp), INTENT(in) :: temp,press ! meteorological parameters:
      ! temp [K],pressure [Pa],RH,u-resp v-comp of mom flux [N/m^2],boundary layer height [m], sensible heat net flux at surface [W/m^2] (GDAS) and downward short wave rad fux (W/m²),rain(mm/h)
      ! and time along trajectory in hour and trajectory time step
      INTEGER, INTENT(in)      :: landuse_index
      REAL, DIMENSION(6), INTENT(in) :: Zin ! Altitude vector in and just above canopy 
      REAL(dp), DIMENSION(8,5) :: z0
      REAL(dp)                 :: ka = 0.4 ! Von Karman constant
      REAL(dp)                 :: dens_air
      REAL(dp), INTENT(in)     :: u
      REAL(dp)                 :: Mair = 28.92D-3
      REAL(dp)                 :: D_height ! Displacement height
      REAL, DIMENSION(6), INTENT(out) :: Wind
      
   ! Rougness length for diffreent land-use catagories and seasons, also table 19.2 [m]
      z0(1,:) = (/ 0.8,0.9,0.9,0.9,0.8 /) ! evergreen, needleleaf trees
      z0(2,:) = (/ 1.05,1.05,0.95,0.55,0.75 /) ! deciduous broadleaf trees
      z0(3,:) = (/ 0.1,0.1,0.05,0.02,0.05 /) ! grass
      z0(4,:) = (/ 0.04,0.04,0.04,0.04,0.04 /) ! desert
      z0(5,:) = (/ 0.1,0.1,0.1,0.1,0.1 /) ! shrubs and interrupted woodlands
      z0(6,:) = (/ 1D-3,1D-3,1D-3,1D-3,1D-3 /) ! sea, from table 8.1 in Jacobson
      z0(7,:) = (/ 0.26,0.26,0.26,0.26,0.26 /) ! small urban area, from table 8.1 in Jacobson
      z0(8,:) = (/ 1D-3,1D-3,1D-3,1D-3,1D-3 /) ! Snow or ice

      Wind=0.      
      D_Height=Zin(6)*(2D0/3D0)

      ! Calculation of the friction velocity using vertical turbulent flux of
      ! horizontal momentum (as computed in HYSPLIT)
      dens_air = Mair*press/(Rg*temp)  ! Air density
      Wind(5:6)=u/ka*LOG((Zin(5:6)-D_height)/z0(landuse_index,season_index))
      Wind(1:4)=0.5*Wind(5)
END SUBROUTINE in_canopy_wind_profile

END PROGRAM adchem1D_new
