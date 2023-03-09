MODULE constants

    USE second_Precision,  ONLY : dp ! KPP Numerical type
    IMPLICIT NONE
    PUBLIC

    
    !!!!!!!!!!!!!!!! Set chosen trajectory:
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CHARACTER(14)      :: date
    NAMELIST /adchem1d_simpar/date
    CHARACTER(3)       :: sim = '39_' ! nr of simulation (for output files)
    INTEGER,PARAMETER  :: nr_st = 1 ! nr of stations the trajectory passes 
     
    INTEGER,PARAMETER  :: days = 7 ! Set the amount of whole days
    INTEGER,PARAMETER  :: hours = 0 ! Set the amount of extra hours that exeeds the days above
    
    INTEGER,PARAMETER :: traj_len_h = days*24+hours+1 ! Length of chosen trajectory in hours
!    INTEGER,PARAMETER :: traj_len_1h = days*8+hours+1 ! Length of chosen trajectory in hours
    INTEGER,PARAMETER :: traj_len_min = (traj_len_h-1)*60+1 ! Length of chosen trajectory in minutes
    INTEGER, PARAMETER :: Nz=40 ! Number of vertical grid cells of the 1D-column model
    REAL(dp), DIMENSION(Nz+1), PARAMETER :: z=(/0,3,9,18,30,45,63,84,108,135,165,198,234,273,315,&
    360,408,459,513,570,630,693,759,828,900,975,1053,1134,1218,1305,1395,1488,1584,1683,&
    1785,1890,1998,2109,2223,2340,2460/)     
    
    REAL(dp), DIMENSION(Nz+2), PARAMETER :: z2=(/0,3,9,18,30,45,63,84,108,135,165,198,234,273,315,&
    360,408,459,513,570,630,693,759,828,900,975,1053,1134,1218,1305,1395,1488,1584,1683,&
    1785,1890,1998,2109,2223,2340,2460,2583/)
    
    
    REAL(dp), PARAMETER :: first_z = 1.5 ! first model height where meteorological parameters are extracted
    

    ! Set particle properties
    INTEGER, PARAMETER :: nr_bins = 100
    ! Set nr of complex reactions in MCM-chemistry (see second_Main.f90)
    !INTEGER, PARAMETER :: NKVALUES = 45
 
    ! Set nr of condensable organic compounds and total number of species in the particles
    INTEGER, PARAMETER :: NCOND = 828 ! nr of condensed species extracted from mcm + POA
    INTEGER, PARAMETER :: NSPEC_P = NCOND+12 ! organic compounds + 9 inorganic + water in inorganic and organic phase
    !INTEGER, PARAMETER :: NSPEC_AQ = 18 ! Number of species in the cloud aquesous phase chemsitry solver + HPMTF + SO2 + H2O2

 
    INTEGER, PARAMETER :: Nlayers = 3    ! Fixed number of particle layers (kinetic multi-layer model)
    INTEGER, PARAMETER :: NMLAYER = Nlayers*NSPEC_P ! Number of differential equations for kinetic multilayer model
    ! Land use catagories and season from Seinfeld and Pandis, table 19.2
    INTEGER, PARAMETER :: coag_index = 1 ! if 1 coagulation is running, otherwise not.
    INTEGER, PARAMETER :: dry_dep_index = 1 ! if 1 dry deposition is running, otherwise not.
    INTEGER, PARAMETER :: dep_g_index = 1 ! if 1, dry dep of some gases running, otherwise not.
	INTEGER, PARAMETER :: index_no_O3_dep = 0 ! if 1, no dry dep of O3
    INTEGER, PARAMETER :: wet_dep_index = 1 ! if 1 wet deposition is running, otherwise not.
    INTEGER, PARAMETER :: thermodyn_index = 1    ! if 1 run the inorganic thermodynamics code and condensation of NH3, HNO3 and HCl
    INTEGER, PARAMETER :: diff_coeff_index = 1 ! if 1, calculate diff coeff, if 0, use a constant k-diff value
    INTEGER, PARAMETER :: condensation_index = 1 ! if 1 condensation is running, otherwise not
    INTEGER, PARAMETER :: fullstat_index = 1 ! if 1 use a full-stationary grid, otherwise a full-moving (full-moving version do not work with nucleation at the moment (Pontus 2015-04-28)
    INTEGER, PARAMETER :: nucleation_index = 1 ! if 1 nucleation is running
    INTEGER, PARAMETER :: index_MEGAN = 0 ! If 1 run MEGAN, else use LPJ-GUESS
	INTEGER, PARAMETER :: index_organic_chem = 0 ! If 1 run with organic gas-phase chemsitry MCM+PRAM
	
    INTEGER, PARAMETER :: index_bin_nucl = 0 ! If 1 run binary water, sulfuric acid npf
	INTEGER, PARAMETER :: sea_spray_index = 3 ! 1:Mårtensson et al., 2003, 2: Salter et al., 2015, 3: Sofiev et al., 2011, 4: Salter distribution but Sofiev PM

    INTEGER, PARAMETER :: index_oligomerization = 0 ! IF 1 run with particle phase oligomerization
 ! season index depending on date (4th column in megan-input is julian date)
    !  season_index = 1 ! midsummer with lush veg (22 june-15 aug)
    !  season_index = 2 ! autumn, not harvested (16 aug-15 sep)
    !  season_index = 3 ! late autumn, no snow, after frost (15 sep-30 okt)
    !  season_index = 4 ! winter, snow (1 nov-31 dec)
    !  season_index = 4 ! winter, snow (1 jan-30 april)                     
    INTEGER, PARAMETER :: season_index = 5 ! transitional (1 may-21 june)
    !INTEGER, PARAMETER :: index_measured_gases = 1 ! IF 1 run with the measured O3, NOx and SO2 as input for the model               
    
    !!! netcdf output related stuff
    integer, parameter :: nr_outputs=14 
    integer, parameter :: days_upwind_save=4 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Some parameters
    REAL(dp), PARAMETER :: Na = 6.022d23   ! Avogadro's constant [molec/mol]
    REAL(dp), PARAMETER :: Rg = 8.3145    ! Universal gas constant [J/(mol K)]
    REAL(dp), PARAMETER :: pi = 2*ASIN(1._dp)
 REAL(dp), PARAMETER :: kb = 1.381D-23;   ! J/K  Boltzmanns constant
 
 ! Species specific particle densities: 
 REAL(dp), PARAMETER :: qNH4=1769.   ! density of NH4 (kg/m^3)
 REAL(dp), PARAMETER :: qDMA=1769.   ! density of NH4 (kg/m^3) 
 REAL(dp), PARAMETER :: qSO4=1769.  ! density of H2SO4 (kg/m^3) 
 REAL(dp), PARAMETER :: qNO3=1725.   ! density of HNO3 (kg/m^3) 
 REAL(dp), PARAMETER :: qNa=2160.   ! density of Na (kg/m^3) 
 REAL(dp), PARAMETER :: qCl=2160.   ! density of Cl (kg/m^3) 
 REAL(dp), PARAMETER :: qHC=1200.   ! density of HC (kg/m^3) 
 REAL(dp), PARAMETER :: qEC=1800.    ! density of EC (kg/m^3) 
 REAL(dp), PARAMETER :: qH2O=1000.   ! density of water (kg/m^3)
 REAL(dp), PARAMETER :: qMetal=2000. ! density of Metal (kg/m^3)
 REAL(dp), PARAMETER :: qMSA=1480. ! density of HIO3(kg/m^3)
 REAL(dp), PARAMETER :: qHIO3=4630. ! density of HIO3 in liquid form (assumed) (kg/m^3)
 REAL(dp), PARAMETER :: qHIO2=4630. ! density of HIO3 in liquid form (assumed) (kg/m^3)
 
 
 
 ! Species specific mass accommodation coefficients: 
 REAL(dp), PARAMETER :: aHC=1.
 REAL(dp), PARAMETER :: aH2SO4=1.
 REAL(dp), PARAMETER :: aDMA=1.
 REAL(dp), PARAMETER :: aHNO3=1.
 REAL(dp), PARAMETER :: aHCl=1.
 REAL(dp), PARAMETER :: aSO2=1.
 REAL(dp), PARAMETER :: aH2O2=1.
 REAL(dp), PARAMETER :: aorg=1.
 
 ! Species specific molar masses:
 REAL(dp), PARAMETER :: Mair=28.98D-3 !kg/mol
 REAL(dp), PARAMETER :: MNO2=46D-3  ! Mol/kg
 REAL(dp), PARAMETER :: MNO=30D-3  ! Mol/kg
 REAL(dp), PARAMETER :: MCO=28D-3  ! Mol/kg
 REAL(dp), PARAMETER :: MO3=48D-3  ! Mol/kg
 REAL(dp), PARAMETER :: MSO2=64D-3  ! Mol/kg
 REAL(dp), PARAMETER :: MSO4=96D-3  ! Mol/kg
 REAL(dp), PARAMETER :: MMSA=96.1D-3 ! Molar mass (kg/mol) of MSA
 REAL(dp), PARAMETER :: MHIO3=175.9D-3 ! Molar mass (kg/mol) of MSA
 REAL(dp), PARAMETER :: MHIO2=159.11D-3 ! Molar mass (kg/mol) of MSA
 REAL(dp), PARAMETER :: MDMA=46D-3 ! Molar mass (kg/mol) of DMA
 REAL(dp), PARAMETER :: MNH4=18D-3  ! Mol/kg
 REAL(dp), PARAMETER :: MNO3=62D-3  ! Mol/kg
 REAL(dp), PARAMETER :: MH2O2=34D-3  ! Mol/kg
 REAL(dp), PARAMETER :: MCl=35.5D-3  ! Mol/kg
 REAL(dp), PARAMETER :: MNa=23D-3  ! Mol/kg
 REAL(dp), PARAMETER :: MH2O=18D-3  ! Mol/kg
 REAL(dp), PARAMETER :: MEC=200D-3  ! Mol/kg
 REAL(dp), PARAMETER :: MMetal=200D-3  ! Mol/kg
 
 
 REAL(dp), PARAMETER :: delta_surf = 5D-10; ! Particle surface (mono-) layer width
 
 ! Logical index with values of 1 or 0 depending on if a process is or is not included 
 INTEGER, PARAMETER :: index_kinetic_multilayer = 0 ! Include kinetic multilayer (3-L) model 
 INTEGER, PARAMETER :: index_layer_diff = 0 ! Consider diffusion transport between the particle layers 
 
 
 ! Coefficients for different nucleation theories
 REAL(dp), PARAMETER :: A_nucl = 1E-8 ! (s^-1) Activation nucleation coefficient J=A_nucl*[H2SO4]
 REAL(dp), PARAMETER :: K_nucl = 1E-14 ! (cm^-3 s^-1) Kinetic nucleation coefficient J=K_nucl*[H2SO4]^2
 REAL(dp), PARAMETER :: K_ELVOC_H2SO4 = 1E-13 ! (cm^-3 s^-1) J=K_ELVOC_H2SO4*[H2SO4]*[ELVOC]

 ! Some other particle phase chemistry or state related properties:
 REAL(dp), PARAMETER :: Ka = 2.5119D-5  ! COOH (p) <-> H+(p) COO-(p)   % (pKa=4.6 for PINONIC and PINIC ACID)
 REAL(dp), PARAMETER :: pCO2 = 410D-6 ! partial pressure CO2 used for pH calculations

 INTEGER :: jjj
 INTEGER, DIMENSION(NSPEC_P-2), PARAMETER :: index_dry= (/1, 2, 3, 4, 5, 6, (jjj, jjj = 9,NSPEC_P) /) ! Particle phase compound index of non-water compounds

 
  
END MODULE constants
