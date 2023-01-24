    MODULE MEGAN_version_2

!***********************************************************************
!
!   CALL:
!      MODULE GAMMA_ETC
!         GAMMA_LAI
!         GAMMA_P
!         GAMMA_TISOP
!         GAMMA_TNISP
!         GAMMA_A
!         GAMMA_S
!
!   Created by Jack Chen 11/04
!   Modified by Tan 11/21/06 for MEGAN v2.0
!
!   History:
!
!       Jun, 2010  modifications to I/O and general structure - Anton Rusanen
!
!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Scientific algorithm
!
!             Emission = [EF] [GAMMA] [RHO]
!                        where 

!                        [EF]    = emission factor (ug/m2h)
!                        [GAMMA] = emission activity factor (non-dimension)
!                        [RHO]   = production and loss within plant canopies (non-dimensino)
!                        Assumption: [RHO] = 1  (11/27/06)   (See PDT_LOT_CP.EXT)
!
!             GAMMA    = [GAMMA_CE] [GAMMA_age] [GAMMA_SM]
!                        where 

!                        [GAMMA_CE]  = canopy correction factor
!                        [GAMMA_age] = leaf age correction factor
!                        [GAMMA_SM]  = soil moisture correction factor
!                        Assumption: [GAMMA_SM]  = 1  (11/27/06)
!             

!             GAMMA_CE = [GAMMA_LAI] [GAMMA_P] [GAMMA_T]
!                        where 

!                        [GAMMA_LAI] = leaf area index factor
!                        [GAMMA_P]   = PPFD emission activity factor
!                        [GAMMA_T]   = temperature response factor
!
!             Emission = [EF] [GAMMA_LAI] [GAMMA_P] [GAMMA_T] [GAMMA_age]
!                        Derivation:
!                        Emission = [EF] [GAMMA_etc] (1-LDF) + [EF] [GAMMA_etc] [LDF] [GAMMA_P]
!                        Emission = [EF] [GAMMA_etc] {(1-LDF) + [LDF] [GAMMA_P]}
!                        Emission = [EF] [GAMMA_etc] {(1-LDF) + [LDF] [GAMMA_P]}
!                        where LDF = light dependent function (non-dimension)
!                                    (See LD_FCT.EXT)
!
!             Final Equation
!             

!             Emission = [EF] [GAMMA_LAI] [GAMMA_T] [GAMMA_age] * { (1-LDF) + [LDF] [GAMMA_P] }
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      USE M2_GAMMA_ETC           ! Module containing gamma functions
      USE M2_RHO_SPC
      USE M2_LD_FCT
      USE M2_SPC_MGN
      
      !How does this thing use the functions from canopy?

      IMPLICIT NONE

      !INCLUDE 'PDT_LOS_CP.EXT'      !  Production and loss within canopy (RHO)
      !INCLUDE 'LD_FCT.EXT'          !  Light dependent factor (LDF)
      !INCLUDE 'SPC_MGN.EXT'         !  MEGAN species
     

      !..  Constants
      ! Length of the time step (days) 

	  ! if the model produces odd values the reason is probably this and/or the removal of the time loop + the layer things I did ( will have done, in the future) in canopy.F

      INTEGER, PARAMETER   :: TSTLEN = 1



      ! parameter for unit conversion
      REAL, PARAMETER    :: ug2tonne = 1E-12  ! convert microgram to metric tonne
      REAL, PARAMETER    :: hr2sec = 3600     ! convert hr to second
      REAL, PARAMETER    :: ug2g = 1E-6       ! convert microgram to gram

      INTEGER, PARAMETER :: NCOLS =1          ! Number of columns
      INTEGER, PARAMETER :: NROWS =1          ! Number of rows

      INTEGER, PARAMETER :: NCOL =1           ! Number of columns for wind
      INTEGER, PARAMETER :: NROW =1           ! Number of rows for wind

      INTEGER, PARAMETER :: NLAYS=1           ! Number of vertival layers
      INTEGER, PARAMETER :: NLAY= 1           ! Number of vertival layers for wind

      LOGICAL, PARAMETER :: TONPHR = .FALSE.  ! output in tons/hr flag

      INTEGER, PARAMETER :: NUM = N_MGN_SPC


      !..  Program I/O files

      CONTAINS 

      SUBROUTINE EMISSION_M2(Day, Layersin, LATin, LONGin, DATEin, TIMEin, PPFDin, D_PPFDin, &
                             TEMPin, DTEMPin, PRESin, HUMin, WINDin, SMOISTin, LADpin, LAIcin, LAIpin, Cantype, ERout, &
                             VARout, ER_Grass_out, ER_SB_out, ER_NT_out, ER_BT_out, GAM_PHO_out, GAM_CE_out, GAM_T_out, &
                             GAM_TMP_out, GAM_OTHER_out, Zin, kzin, EMI, sun_par, Sunleaftk, Shadeleaftk, sunfrac,EFin)
		


      INTEGER, INTENT(IN) ::  Layersin ,kzin        ! Number of Layers with vegetation in them, presumed from ground 0 to Layersin.
	  INTEGER, INTENT(IN) ::  DATEin, TIMEin        ! Date (yyyddd) and time (hhmmss)



      REAL, PARAMETER     ::  Avog = 6.0221E23      ! Avogadro-number [molecules cm3]



!      REAL                ::  EFin(366,22)          ! emission factors, 1-dimensional array, 20 slots by default
      REAL, INTENT(in)    ::  EFin(22)               ! Standard emission factor from (http://lar.wsu.edu/megan/guides.html) 

	  REAL, INTENT(IN)    ::  LATin, LONGin         ! LATitutde and LONGitude
	  REAL, INTENT(IN)    ::  DTEMPin               ! TEMPerature (K) and Daily average TEMPerature (K)

	  REAL, INTENT(IN)    ::  PPFDin, D_PPFDin      ! PAR and Daily average of PAR in (µmol/m2/s)
	  REAL, INTENT(IN)    ::  PRESin                ! Pressure (Pa)

	  REAL, INTENT(IN)    ::  LAIcin , LAIpin       ! LAIc is LAI for current month, LAIp for previous

	  REAL, INTENT(IN)    ::  SMOISTin              ! Wind speed (m/s) and Soil moisture, (a number belonging to [0,1])

	  !REAL, INTENT(IN)    ::  SRADin, DSRADin       ! Solar RADiation (W/mÂ²) and Daily average Solar RADiation (W/mÂ²)


	  REAL, INTENT(IN), DIMENSION(kzin) ::  LADpin  ! Leaf area density as [0,1] in the canopy

	  REAL, INTENT(IN), DIMENSION(:)    ::  TEMPin
	  REAL, INTENT(IN), DIMENSION(:)    ::  HUMin   ! Humidity (%)
	  REAL, INTENT(IN), DIMENSION(:)    ::  Zin     ! Layer top array

	   
	  REAL, DIMENSION(:), INTENT(IN) :: WINDin 
	   
	  REAL, INTENT(INOUT) :: ERout(N_MGN_SPC, kzin) ! Emission rates output (normally g/s, can be set with a boolean to ton/hr)


	  CHARACTER*16       VARout(N_MGN_SPC)       !VOC names in the same order as ERout


      REAL, DIMENSION(N_MGN_SPC, kzin),    INTENT(OUT) ::  ER_Grass_out, ER_SB_out, ER_NT_out, ER_BT_out  !emissions by plant functional type
      REAL, DIMENSION(4, kzin), INTENT(OUT)            ::  GAM_PHO_out, GAM_CE_out, GAM_T_out          ! Gamma factors for photons, CE?, and temperature
      REAL, DIMENSION(N_MGN_SPC, kzin),    INTENT(OUT) ::  GAM_TMP_out                                 ! Gamma factor for what?
      REAL, DIMENSION(N_MGN_SPC, 3),       INTENT(OUT) ::  GAM_OTHER_out                               !gamma factors soil moisture, leaf age and lai correction
      



      REAL,DIMENSION(kzin) ::  Sunleaftk, Shadeleaftk, sunfrac, sun_par



	  REAL, DIMENSION(Layersin) :: midpoints !Layer centerpoints for gamme_CE
      
      REAL, DIMENSION(N_MGN_SPC, kzin) :: ERtemporary

      REAL :: ISO_EMI, MBO_EMI, MYR_EMI, SAB_EMI, LIM_EMI, CAR_EMI, OCI_EMI, BPI_EMI, API_EMI, FAR_EMI, BCA_EMI, & 

	          MET_EMI, ACT_EMI, ACA_EMI, FOR_EMI, CH4_EMI, NO_EMI,  OMT_EMI, OSQ_EMI, CO_EMI,  LIN_EMI, CIN_EMI, &

			  emi_factor, time



      REAL, DIMENSION(kzin,22) :: EMI


      ! Program name
      CHARACTER*16  :: PROGNAME = 'MEGAN'
       
      !..  model parameters
      INTEGER       SDATE                            ! Start date YYYYDDD
      INTEGER       STIME                            ! Start time HHMMSS
      ! INTEGER       RLENG                          ! Run length HHMMSS
      INTEGER       TSTEP                            ! time step
      LOGICAL  ::   TONPHR = .FALSE.                 ! output in tons/hr flag
      ! LOGICAL  ::   ONLN_DTEMP = .TRUE.            ! online daily average
                                                     ! temperature calculation
      ! LOGICAL  ::   ONLN_DSRAD = .TRUE.            ! online daily average
                                                     ! solar radiation calculation
      !I/O API file parameters
      ! INTEGER       JDATE                          ! Date YYYYDDD from ECMAP
      ! INTEGER       JTIME                          ! Time HHMMSS from ECMAP
      ! INTEGER       MXREC                          ! total number of timesteps

      !..  Internal parameters
      ! internal paramters (status and buffer)
      INTEGER       IOS                              ! i/o status
      CHARACTER*256 MESG                             ! message buffer
      
      ! parameter for output species
      INTEGER, PARAMETER :: NEMIS = N_MGN_SPC
                                                     ! number of output emission variable
                                                     ! number of MEGAN species

      ! local variables and their descriptions:
      REAL          GAREA                            ! Area in one grid (meter^2)
      ! INTEGER       AVEBY                          ! Divider for daily average

      REAL          LDF                              ! Light dependent factor
      REAL          RHO                              ! Production and loss within canopy

      REAL, ALLOCATABLE    :: ER_BT( :,: )           ! output emission buffer
      REAL, ALLOCATABLE    :: ER_NT( :,: )           ! output emission buffer
      REAL, ALLOCATABLE    :: ER_SB( :,: )           ! output emission buffer
      REAL, ALLOCATABLE    :: ER_Grass( :,: )           ! output emission buffer  
      REAL, ALLOCATABLE    :: EF( :,: )              ! input annual emission factor
      REAL, ALLOCATABLE    :: LAT( :,: )             ! input latitude of grid cell
      REAL, ALLOCATABLE    :: LONG( :,: )            ! input longitude of grid cell

      REAL, ALLOCATABLE    :: LAIp( :,: )            ! previous monthly LAI
      REAL, ALLOCATABLE    :: LAIc( :,: )            ! current monthly LAI

      REAL, ALLOCATABLE    :: TEMP(:)                ! input hourly temperature (K)
      REAL, ALLOCATABLE    :: dummyTEMP(:, :, :)     ! for gamma_TNSIP
      REAL, ALLOCATABLE    :: SRAD( :,: )            ! input hourly radiation
      REAL, ALLOCATABLE    :: PPFD( :,: )            ! calculated PAR (umol/m2.s)

      REAL, ALLOCATABLE    :: D_SRAD( :,: )          ! daily solar radiation (Watt/m2)
      REAL, ALLOCATABLE    :: D_PPFD( :,: )          ! daily PAR (umol/m2.s)
      REAL, ALLOCATABLE    :: D_TEMP( :,: )          ! input daily temperature (K)

      REAL, ALLOCATABLE    :: GAM_PHO(:, :,: )       ! light correction factor
      REAL, ALLOCATABLE    :: GAM_TMP( :,:,: )       ! temperature correction factor
      REAL, ALLOCATABLE    :: GAM_LHT( :,: )         ! LAI correction factor
      REAL, ALLOCATABLE    :: GAM_AGE( :,: )         ! leaf age correction factor
      REAL, ALLOCATABLE    :: GAM_SMT( :,: )         ! Soil moilture correction factor
      REAL, ALLOCATABLE    :: GAM_T(:, :,: )
      REAL, ALLOCATABLE    :: GAM_CE (:, :,: )
      REAL, ALLOCATABLE    :: Wind(:) 
      REAL, ALLOCATABLE    :: Pres( :,: ) 
      REAL, ALLOCATABLE    :: Humidity(:) 
      REAL, ALLOCATABLE    :: DI( :,: ) 
      REAL, ALLOCATABLE    :: SMOIS( :,: )

      REAL, ALLOCATABLE    :: GAMMAfactors( :,:, :)  ! gamma factors (  GAMMA_T, GAMMA_PHO, GAMMA_CE), In layers, for different plant types 

      REAL, ALLOCATABLE    :: LAdp(:)



     ! INTEGER, ALLOCATABLE :: Cantype( :,: ) 
      INTEGER, INTENT(in) :: Cantype 
      INTEGER, ALLOCATABLE :: ISLTYP( :,: )
 
      INTEGER          ILINE                         ! current line
      CHARACTER(LEN=1000) LINE                       ! input line buffer
      INTEGER       CID, INX, INY                    ! Input grid x and y
      INTEGER, PARAMETER :: MXTCOL = 9               ! Columns in an input line
      CHARACTER*30     SEGMENT( MXTCOL )             ! Input line fields

      INTEGER, PARAMETER :: NVARS = MXTCOL - 3       ! Number of LAT, LONG, and
                                                     ! PFT factor variables
      CHARACTER*16 VNAME( NVARS )                    ! Variable names
      INTEGER, PARAMETER :: NPFT = 4                 ! guess, number of plant functional types
      REAL, ALLOCATABLE :: PFTFBUF( :, :, : )        ! PFT factor array+Lat and LONG
      REAL, ALLOCATABLE :: PFTF( :, :, : )           ! PFT factor array

      CHARACTER*16   VARIABLENAMES( 30 ) 

      INTEGER, PARAMETER:: NrTyp = 4, NrCha = 16
      REAL,DIMENSION(NrCha,NrTyp ):: Canopychar
        
      ! REAL, DIMENSION(:,:) :: Wind, Pres,Humidity,DI
      ! REAL,DIMENSION(NrCha,NrTyp) :: Canopychar 
      ! INTEGER,DIMENSION(:,:) :: Cantype
      

	  ! loop indices
      INTEGER       T, S, I, J ,K, VAR, VAR2, lay     ! Counters
      INTEGER       NMAP            ! Index
      INTEGER       IDATE           ! Looping date YYYYDDD
      INTEGER       ITIME           ! Looping time HHMMSS

      ! times
      INTEGER       MON             ! Month from YYYYDDD
      INTEGER       DAY             ! Day from YYYYDDD

      INTEGER :: test1,test2,test3

      ! months
      CHARACTER*3   MONTHS( 12 )
      DATA          MONTHS &
      &  / 'JAN' , 'FEB' , 'MAR' , 'APR' , 'MAY' , 'JUN' , &
      &    'JUL' , 'AUG' , 'SEP' , 'OCT' , 'NOV' , 'DEC'   /
      CHARACTER*2   MONNUM( 12 )
      DATA          MONNUM &
      &  /  '1 ' ,  '2 ' ,  '3 ' ,  '4 ' ,  '5 ' ,  '6 ' , &
      &     '7 ' ,  '8 ' ,  '9 ' ,  '10' ,  '11' ,  '12'   /


!#ifdef PARALLEL
!
!    CHARACTER(LEN=*), PARAMETER :: &
!
!         filename1 = 'sosa_in', &
!
!        filename2 = 'sosa_out'

!#elif LINUX

!    CHARACTER(LEN=*), PARAMETER :: &

!         filename1 = 'sosa_in', &

!         filename2 = 'sosa_out'

!    CHARACTER(LEN=*), PARAMETER :: &
!         filename1 = 'c://HY-data/BOY/documents/sosa_in', &
!         filename2 = 'c://HY-data/BOY/documents/sosa_out'
!#else
!
!    CHARACTER(LEN=*), PARAMETER :: &
!
!         filename1 = 'c://michael/sosa_in', &
!
!         filename2 = 'c://michael/sosa_out'

!#endif




!**********************************************************************

!======================================================================
!..  Begin program
!======================================================================
!----------------------------------------------------------------------
!....1) File set up and assign I/O parameters
!----------------------------------------------------------------------

      GAREA = 1 !1.296E9 !mÂ² 

      !..  Get input parameters

      SDATE= DATEin !2001187 !'Model start date (YYYYDDD)'
      STIME= TIMEin !'Model start time (HHMMSS)'
      !TSTEP= 1000 !'Model time step (HHMMSS)' ! not used

      DO S = 1, NEMIS
         VARIABLENAMES( S ) = TRIM( MGN_SPC( S ) )
      ENDDO


      VARIABLENAMES(NEMIS + 1) = 'D_TEMP'

      VARIABLENAMES(NEMIS + 2) = 'D_PPFD'

!----------------------------------------------------------------------
!....2) Process emission rates
!----------------------------------------------------------------------
!..  Allocate memory

      ALLOCATE ( ER_BT( NCOLS, NROWS ), STAT = IOS )
      ALLOCATE ( ER_NT( NCOLS, NROWS ), STAT = IOS )
      ALLOCATE ( ER_SB( NCOLS, NROWS ), STAT = IOS )
      ALLOCATE ( ER_Grass( NCOLS, NROWS ), STAT = IOS )
      ALLOCATE ( EF( NCOLS, NROWS ), STAT = IOS )
      ALLOCATE ( LAT( NCOLS, NROWS ), STAT = IOS ) 
      ALLOCATE ( LONG( NCOLS, NROWS ), STAT = IOS ) 
      ALLOCATE ( LAIp( NCOLS, NROWS ), STAT = IOS )
      ALLOCATE ( LAIc( NCOLS, NROWS ), STAT = IOS )
      ALLOCATE ( D_SRAD( NCOLS, NROWS ), STAT = IOS )
      ALLOCATE ( D_PPFD( NCOLS, NROWS ), STAT = IOS )
      ALLOCATE ( D_TEMP( NCOLS, NROWS ), STAT = IOS )
      ALLOCATE ( SRAD( NCOLS, NROWS ), STAT = IOS )
      ALLOCATE ( PPFD( NCOLS, NROWS ), STAT = IOS )
      ALLOCATE ( TEMP((Layersin+1) ), STAT = IOS ) 
      ALLOCATE (  dummyTEMP(NCOLS, NROWS, Layersin) )
      ALLOCATE ( GAM_PHO(NPFT, NCOLS, NROWS ), STAT = IOS )
      ALLOCATE ( GAM_TMP( NCOLS, NROWS, Layersin ), STAT = IOS )
      ALLOCATE ( GAM_LHT( NCOLS, NROWS ), STAT = IOS )
      ALLOCATE ( GAM_AGE( NCOLS, NROWS ), STAT = IOS )
      ALLOCATE ( GAM_SMT( NCOLS, NROWS ), STAT = IOS )
      ALLOCATE ( GAM_T(NPFT,NCOLS, NROWS ), STAT = IOS )      
      ALLOCATE ( Wind((Layersin+1)), STAT = IOS )
      ALLOCATE ( Pres( NCOLS, NROWS ), STAT = IOS ) 
      ALLOCATE ( Humidity((Layersin+1) ), STAT = IOS )
      ALLOCATE ( DI( NCOLS, NROWS ), STAT = IOS )
      !ALLOCATE ( Cantype ( NCOLS, NROWS ), STAT = IOS )
      ALLOCATE ( GAM_CE (NPFT, NCOLS, NROWS ), STAT = IOS )
      ALLOCATE ( ISLTYP ( NCOLS, NROWS ), STAT = IOS )
      ALLOCATE ( SMOIS ( NCOLS, NROWS ), STAT = IOS )
      

	  ! Allocate memory for PFTF
      ALLOCATE ( PFTF( NPFT, NCOLS, NROWS ), STAT = IOS )
      ALLOCATE (GAMMAfactors(3,4, Layersin))
      ALLOCATE (LADp(Layersin))





      !...  Get input file unit
	
      ! 1 = LAT, 2 = LONG,
      ! 3 = PFTF_BT, 4 = PFTF_NT,
      ! 5 = PFTF_SB, 6 = PFTF_HB
      !!! Order does matter !!!
      !well, not the same order anymore

       
      !.. Partial fractions of different plant types 
      !.. (doesn't need to sum to 100%, because parts of the grid can be barren)
      !.. You should hardcode appropriate values for your model here, or this can be added to input parameters
      PFTF(Cantype, 1, 1) = 100.
      PFTF = 100.

      !..  Get D_PPFD, D_TEMP, LAT and LONG (constant for each month)
      
      LAT    = LATin   ! 0 !latitude
      LONG   = LONGin  ! 0 !longitude
      D_TEMP = DTEMPin ! 240 ! Daily average temperature K
      D_PPFD = D_PPFDin ! 100.0 ! Daily average PAR
      !D_SRAD = DSRADin ! 100.0 ! Daily average shortwave radiation W/m2

      !.. Calculate daily PPDF
      !   PPFD: SRAD - short wave from sun (W/m2)
      !   assuming 4.766 (umol m-2 s-1) per (W m-2)
      !   assume 1/2 of SRAD is in 400-700nm band
      !D_PPFD = D_SRAD * 4.766 * 0.5

!--------INPUT OTHER VARIABLES----------------
        
      DI     = 0.0       ! something related to DIstomata, value 0 came with model
      SMOIS  = SMOISTin  ! Soil moisture
 	  ISLTYP = 6         ! Soiltype (I don't know if 6 is right, infact I don't know what type of soil it is, value came with model)
	
	  midpoints= centralizer(Zin, Layersin)

	  !rescaling to [0,1]
	  midpoints= midpoints/(Zin(Layersin))


      !--------INPUT standard emission potential
      ! 1 =  Isoprene
      ! 2 =  MBO (2methyl-3buten-2ol)
      ! 3 =  MYRC
      ! 4 =  Sabinene
      ! 5 =  Limonen
      ! 6 =  3-Carene
      ! 7 =  Ocimene
      ! 8 =  Beta-pinene
      ! 9 =  Alpha-pinene
      ! 10 = FARN
      ! 11 = Betacarophylene
      ! 12 = Methanol
      ! 13 = Aceton
      ! 14 = Acetaldehyde
      ! 15 = Formaldehyde
      ! 16 = Methan
      ! 17 = NO
      ! 18 = Other monoterpene
      ! 19 = Other sesquiterpenes
      ! 20 = CO
      ! 21 = Cineole
      ! 22 = Linalool

    !! IF (time .EQ. 0.0) THEN
    
    !     OPEN(unit=5,file=''//filename1//'/General/Hyytiala/EF_day.txt')
    !     DO J = 1,366
    !        READ(5,*) (EFin(J,I), I = 1,22)
    !     ENDDO
    !     CLOSE(5)
    !! ENDIF


      !--------INPUT canopy----------------'


      ! Input of the Canopy characteristics for NrTyp with:
     
      OPEN(44,file='input/megan/canopy.txt')

         DO J = 1,NrCha

           READ(44,*) (Canopychar(J,I), I = 1, NrTyp)

	     ENDDO

      CLOSE(44)


      IDATE = SDATE
      ITIME = STIME


      !..Initialize hourly variables
      ! TEMP(:) = TEMPin(1:(Layersin+1))!280 ! hourly temperature K

      !TEMP(1:Layersin)  = TEMpin(Layersin : 1 : -1)
      TEMP(1:Layersin)  = TEMpin(1:Layersin)
      TEMP( Layersin+1) = TEMpin(Layersin+1)

      !  dummyTEMP(1,1,:) = TEMPin(1:Layersin) !because calculation of Gamma_tmp requires a dimension(:,:,:) array
         
      dummyTEMP(1, 1 , 1 : Layersin) = TEMPin( Layersin : 1 : -1)
       !SRAD     =  SRADin                   ! hourly shortwave radiation W/m2
       PPFD     =  PPFDin                   ! hourly PAR
       LAIp     =  LAIpin                   ! Previous month LAI
       LAIc     =  LAIcin                   ! Current month LAI
       Wind(:)  =  WINDin( 1: (Layersin+1)) ! m/s


       Wind(1:Layersin) = Windin( Layersin : 1 : -1)
       Wind(Layersin+1) = Windin(Layersin+1)
         
       Pres =  PRESin !1013E3 ! Pa? probably
       Humidity(:) = HUMin(1: (Layersin+1)) ! %?
         
       !Humidity(1:Layersin) = HUMin( Layersin : 1 : -1)
       Humidity(1:Layersin) = HUMin(1:Layersin)
       Humidity(Layersin+1) = HUMin(Layersin+1)

       ! LAdp = LAdpin(1:Layersin)
       Ladp= Ladpin( LAyersin : 1 : -1 )
	
	    
       !..Cal Hourly variables
       !  PPFD: SRAD - short wave from sun (W/m2)
       !  assuming 4.766 (umol m-2 s-1) per (W m-2)
       !  assume 1/2 of SRAD is in 400-700nm band
       !PPFD = SRAD * 4.766 * 0.5
         
       ! We could of course give PPFD, instead of SRAD

       !..Go over all the chemical species
       
    
       DO S = 1, NEMIS


       !..Initialize variables
          EF = 0.
          ER_BT = 0.
          ER_NT = 0.
          ER_SB = 0. 
          ER_Grass = 0.
          GAM_PHO = 0.
          GAM_TMP = 0.
          GAM_LHT = 0.
          GAM_AGE = 0.
          GAM_SMT = 0.
          GAM_T = 0.
          GAM_CE = 0.

          !..Get EF
          ! EFBUF = 'EF_'//TRIM(VARIABLENAMES(S))
		  EF = 1. !EFin(S)
		
		  Gammafactors = -1.0

          !..Select algorithms for differet chemical species
          !  Due to the differences in the calculation for gamma,
          !  this logical condition will check for species and choose
          !  (call) the right algorithms to calculate the gamma and
          !  return back the main.


          CALL  GAMME_CE(IDATE,ITIME,NCOLS,NROWS,LAT, &
                         LONG,PPFD,TEMP,LAIc,Wind,Pres,Humidity, &
                         DI, Cantype,PFTF,Canopychar,NrCha,NrTyp,NPFT, &
                         GAM_T,GAM_PHO,GAM_CE, LAdp, midpoints, &
                         Layersin, GAMMAfactors, sun_par, &
                         Sunleaftk, Shadeleaftk, sunfrac)
     

          SELECT CASE ( VARIABLENAMES(S) )


          CASE ( 'ISOP' )

              CALL GAMMA_TISOP( NCOLS, NROWS, TEMP, D_TEMP, GAM_T )
              DO I=1, NCOLS
                DO J=1,NROWS
                     GAM_T  = GAM_T
                 ENDDO
              ENDDO

          CASE ('MBO','MYRC', 'SABI', 'LIMO', '3CAR', 'OCIM', 'BPIN', 'APIN', 'FARN', &
          'BCAR', 'MEOH', 'ACTO','ACTA','FORM', 'CH4', 'NO', 'OMTP','OSQT','CO')

	         CALL GAMMA_TNISP(NCOLS, NROWS, VARIABLENAMES(S), TEMP, Layersin, GAM_TMP)

          CASE DEFAULT
         

		     MESG = 'Error: Chemical species, invalid variable: ' //TRIM(VARIABLENAMES(S))

	         WRITE(*,*) MESG
	         STOP
          ENDSELECT

          CALL GAMMA_LAI(NCOLS, NROWS, LAIc, GAM_LHT)

          ! CALL GAMMA_P(IDATE, ITIME, LAT, LONG, NCOLS, NROWS, PPFD, D_PPFD, GAM_PHO)

          CALL GAMMA_A(IDATE, ITIME, NCOLS, NROWS, VARIABLENAMES(S), LAIp, LAIc, TSTLEN, D_TEMP, GAM_AGE)

          ! CALL GAMMA_S( NCOLS, NROWS, GAM_SMT )


          CALL  GAMMA_S(NCOLS, NROWS, SMOIS, ISLTYP, GAM_SMT)


          DO VAR = 1, NEMIS
		     IF((TRIM(VARIABLENAMES(S))) .EQ. (TRIM(LDF_SPC(VAR)))) THEN
		         LDF= LDF_FCT(VAR)
		     END IF

		     IF((TRIM(VARIABLENAMES(S))) .EQ. (TRIM(RHO_SPC(VAR)))) THEN
		         RHO= RHO_FCT(VAR)
		     END IF
	      END DO	     


          !Same thing but separately for every layer.

          SELECT CASE ( VARIABLENAMES(S) )
            
            CASE ('ISOP')
            
               DO I=1,NCOLS
                  DO J=1,NROWS
                     DO lay=1, Layersin
                        IF (TONPHR) THEN
                           ER_BT_out(S,lay) = (GAMMAfactors(3,1,lay) * GAM_AGE(I,J) * GAM_SMT(I,J) * RHO) * GAREA * ug2tonne
                           ER_NT_out(S,lay) = (GAMMAfactors(3,2,lay) * GAM_AGE(I,J) * GAM_SMT(I,J) * RHO) * GAREA * ug2tonne
                           ER_Grass_out(S,lay) = (GAMMAfactors(3,3,lay) * GAM_AGE(I,J) * GAM_SMT(I,J) * RHO) * GAREA * ug2tonne
                           ER_SB_out(S,lay) = (GAMMAfactors(3,4,lay) * GAM_AGE(I,J) * GAM_SMT(I,J) * RHO) * GAREA * ug2tonne
                        ELSE
                           ER_BT_out(S,lay) = (GAMMAfactors(3,1,lay) * GAM_AGE(I,J) * GAM_SMT(I,J) * RHO)
                           ER_NT_out(S,lay) = (GAMMAfactors(3,2,lay) * GAM_AGE(I,J) * GAM_SMT(I,J) * RHO)
                           ER_Grass_out(S,lay) = (GAMMAfactors(3,3,lay) * GAM_AGE(I,J) * GAM_SMT(I,J) * RHO)
                           ER_SB_out(S,lay) = (GAMMAfactors(3,4,lay) * GAM_AGE(I,J) * GAM_SMT(I,J) * RHO)
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
        
		    CASE ('MBO','MYRC','SABI','LIMO','3CAR','OCIM','BPIN', 'APIN','FARN','BCAR','MEOH','ACTO','ACTA','FORM', &
                  'CH4',  'NO','OMTP','OSQT','CO')
          

		       DO I=1,NCOLS
                  DO J=1,NROWS
                     Do lay=1, Layersin
                        IF ( TONPHR ) THEN
                           ER_BT_out(S,lay) = (EF(I,J) * GAM_TMP(I,J, lay) * GAM_AGE(I,J) * GAM_LHT(I,J) * GAM_SMT(I,J) * RHO) * &

						                      ((1-LDF) + (GAMMAfactors(2,1,lay)*LDF)) * GAREA * ug2tonne
                           

						   ER_NT_out(S,lay) = (EF(I,J) * GAM_TMP(I,J, lay) * GAM_AGE(I,J) * GAM_LHT(I,J) * GAM_SMT(I,J) * RHO) * &
                                              ((1-LDF) + (GAMMAfactors(2,2,lay)*LDF)) * GAREA * ug2tonne
                           

						   ER_SB_out(S,lay) = (EF(I,J) * GAM_TMP(I,J, lay) * GAM_AGE(I,J) * GAM_LHT(I,J) * GAM_SMT(I,J) * RHO) * &
                                              ((1-LDF) + (GAMMAfactors(2,3,lay)*LDF)) * GAREA * ug2tonne
                           

						   ER_Grass_out(S,lay) = (EF(I,J) * GAM_TMP(I,J, lay) * GAM_AGE(I,J) * GAM_LHT(I,J) * GAM_SMT(I,J) * RHO) * &
                                              ((1-LDF) + (GAMMAfactors(2,4,lay)*LDF)) * GAREA * ug2tonne


                        ELSE
! Original MEGAN:                        
!                           ER_BT_out(S,lay) = (EF(I,J) * GAM_TMP(I,J,lay) * GAM_AGE(I,J) * GAM_LHT(I,J) * GAM_SMT(I,J) * RHO) * &
!                                              ((1-LDF) + (GAMMAfactors(2,1,lay)*LDF))
! Dittes suggestion (divide the emissions into pool and synthesised instead of correcting the pool emissions with a light-dependent factor)
! The GAMMAfactors(3,x,lay) describes the light-dependency of the isoprene emissions but should work also for MT. 
! See Ghirardo et al. Plant, Cell and Environment (2010) 33, 781Ð792 for details about the factors of light and temperature dependent emissions

						   ER_BT_out(S,lay) = 0.02*(EF(I,J) * GAM_TMP(I,J,lay) * GAM_AGE(I,J) * GAM_LHT(I,J) * GAM_SMT(I,J) * RHO) + &
                                              0.98*(GAMMAfactors(3,1,lay) * GAM_AGE(I,J) * GAM_SMT(I,J) * RHO) 

! Original MEGAN:
!						   ER_NT_out(S,lay) = (EF(I,J) * GAM_TMP(I,J,lay) * GAM_AGE(I,J) * GAM_LHT(I,J) * GAM_SMT(I,J) * RHO) * &
!                                              ((1-LDF) + (GAMMAfactors(2,2,lay)*LDF))
! Dittes suggestion (divide the emissions into pool and synthesised instead of correcting the pool emissions with a light-dependent factor)
						   ER_NT_out(S,lay) = 0.42*(EF(I,J) * GAM_TMP(I,J,lay) * GAM_AGE(I,J) * GAM_LHT(I,J) * GAM_SMT(I,J) * RHO) + & ! *0.42
                                              0.58*(GAMMAfactors(3,2,lay) * GAM_AGE(I,J) * GAM_SMT(I,J) * RHO) ! *0.58
                                              
                                              
! Original MEGAN:
!						   ER_Grass_out(S,lay) = (EF(I,J) * GAM_TMP(I,J,lay) * GAM_AGE(I,J) * GAM_LHT(I,J) * GAM_SMT(I,J) * RHO) * &
!                                              ((1-LDF) + (GAMMAfactors(2,3,lay)*LDF))
! Dittes suggestion (divide the emissions into pool and synthesised instead of correcting the pool emissions with a light-dependent factor)
						   ER_SB_out(S,lay) = 0.02*(EF(I,J) * GAM_TMP(I,J,lay) * GAM_AGE(I,J) * GAM_LHT(I,J) * GAM_SMT(I,J) * RHO) + &
                                              0.98*(GAMMAfactors(3,3,lay) * GAM_AGE(I,J) * GAM_SMT(I,J) * RHO)                                               

! Original MEGAN:
!                           ER_SB_out(S,lay) = (EF(I,J) * GAM_TMP(I,J,lay) * GAM_AGE(I,J) * GAM_LHT(I,J) * GAM_SMT(I,J) * RHO) * &
!                                              ((1-LDF) + (GAMMAfactors(2,4,lay)*LDF))
! Dittes suggestion (divide the emissions into pool and synthesised instead of correcting the pool emissions with a light-dependent factor)
						   ER_SB_out(S,lay) = 0.02*(EF(I,J) * GAM_TMP(I,J,lay) * GAM_AGE(I,J) * GAM_LHT(I,J) * GAM_SMT(I,J) * RHO) + &
                                              0.98*(GAMMAfactors(3,4,lay) * GAM_AGE(I,J) * GAM_SMT(I,J) * RHO) 
                           

                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
       
	        CASE DEFAULT
               WRITE(*,*) 'Error: Chemical species, invalid variable: '
               WRITE(*,*) TRIM(VARIABLENAMES(S))
	           STOP
            

          ENDSELECT

!----------------------------------------------------------------------
!....3) Write out the calculated ER and met data
!----------------------------------------------------------------------
          ! Write emission to file

          IF (S <= N_MGN_SPC) THEN
             VARout(S) = VARIABLENAMES( S)
          END IF



    
          IF (S .EQ. 9) THEN
             DO lay = 2, 2 !looping over plant functional types
                DO J = 1,layersin
                   GAM_T_out(lay, J)     = GAMMAfactors(1, lay, (layersin-J+1)) !Note GAM_T is not actually used anywhere
                   GAM_PHO_out(lay, J)   = GAMMAfactors(2, lay, (layersin-J+1))
                   GAM_CE_out(lay, J)    = GAMMAfactors(3, lay, (layersin-J+1))
                ENDDO
             END DO
          ENDIF

          GAM_OTHER_out(S, 1) = GAM_AGE(1,1)
          GAM_OTHER_out(S, 2) = GAM_LHT(1,1)
          GAM_OTHER_out(S, 3) = GAM_SMT(1,1)
    
          !GAM_TMP_out(S, 1 : Layersin) = GAM_TMP(1,1, Layersin : 1 : -1)
          GAM_TMP_out(S, 1 : Layersin) = GAM_TMP(1,1, 1:Layersin)

      ENDDO ! End loop for emission species (S)



IF  (Cantype == 1) THEN
      ERtemporary = ER_NT_out
END IF      

IF  (Cantype == 2) THEN
      ERtemporary = ER_BT_out
END IF
      
IF  (Cantype == 3) THEN
      ERtemporary = ER_Grass_out
END IF      

IF  (Cantype == 4) THEN
      ERtemporary = ER_SB_out
END IF


!      ERtemporary = ER_Grass_out + ER_NT_out + ER_SB_out + ER_BT_out

      !ERtemporary = ER_NT_out

    

      !Erout= Ertemporary



      ERout( : , 1 : Layersin) = Ertemporary( :, Layersin : 1 : -1)

      !Rout = Ertemporary


      ! Standard emission potential [ng/g(needledryweight)/h into g/cm2/s] for Hyytiala spring from Hakola et al. Biogeosc., 3, 2006

      ! (with g(needledryweight into cm2 by CarboEuroflux data: Standing leaf biomass Hyytiala 0.538 kg/m2 or 0.0538 g/cm2)

! Use the standard emission factor from http://lar.wsu.edu/megan/guides.html (ug/m^2/h) and 
! convert it to values given in g/cm^2/s


      ISO_EMI = EFin(1)*1D-10/3.6D3 !EFin(Day,1)  * 1E-9 / 3600 * 0.0538
      
      MBO_EMI =  EFin(2)*1D-10/3.6D3 !EFin(Day,2)  * 1E-9 / 3600 * 0.0538
      
      MYR_EMI = EFin(3)*1D-10/3.6D3 ! EFin(Day,3)  * 1E-9 / 3600 * 0.0538

      SAB_EMI = EFin(4)*1D-10/3.6D3 ! EFin(Day,4)  * 1E-9 / 3600 * 0.0538

      LIM_EMI = EFin(5)*1D-10/3.6D3 !EFin(Day,5)  * 1E-9 / 3600 * 0.0538

      CAR_EMI = EFin(6)*1D-10/3.6D3 !EFin(Day,6)  * 1E-9 / 3600 * 0.0538

      OCI_EMI =  EFin(7)*1D-10/3.6D3 !EFin(Day,7)  * 1E-9 / 3600 * 0.0538

      BPI_EMI =  EFin(8)*1D-10/3.6D3 !EFin(Day,8)  * 1E-9 / 3600 * 0.0538

      API_EMI =  EFin(9)*1D-10/3.6D3 !EFin(Day,9)  * 1E-9 / 3600 * 0.0538

      FAR_EMI =  EFin(10)*1D-10/3.6D3 !EFin(Day,10) * 1E-9 / 3600 * 0.0538

      BCA_EMI =  EFin(11)*1D-10/3.6D3 !EFin(Day,11) * 1E-9 / 3600 * 0.0538
     
      MET_EMI =  EFin(12)*1D-10/3.6D3 !EFin(Day,12) * 1E-9 / 3600 * 0.0538

      ACT_EMI =  EFin(13)*1D-10/3.6D3 !EFin(Day,13) * 1E-9 / 3600 * 0.0538

      ACA_EMI =  EFin(14)*1D-10/3.6D3 !EFin(Day,14) * 1E-9 / 3600 * 0.0538
      
      FOR_EMI =  EFin(15)*1D-10/3.6D3 !EFin(Day,15) * 1E-9 / 3600 * 0.0538
      
      CH4_EMI =  EFin(16)*1D-10/3.6D3 !EFin(Day,16) * 1E-9 / 3600 * 0.0538
      
      NO_EMI  =  EFin(17)*1D-10/3.6D3 !EFin(Day,17) * 1E-9 / 3600 * 0.0538

      OMT_EMI =  EFin(18)*1D-10/3.6D3 !EFin(Day,18) * 1E-9 / 3600 * 0.0538
      
      OSQ_EMI =  EFin(19)*1D-10/3.6D3 !EFin(Day,19) * 1E-9 / 3600 * 0.0538
      
      CO_EMI  =  EFin(20)*1D-10/3.6D3 !EFin(Day,20) * 1E-9 / 3600 * 0.0538

      CIN_EMI =  EFin(21)*1D-10/3.6D3 !EFin(Day,21) * 1E-9 / 3600 * 0.0538

      LIN_EMI =  EFin(22)*1D-10/3.6D3 !EFin(Day,22) * 1E-9 / 3600 * 0.0538



      ! Emission with EF in g/cm2/s with ER into molecules/cm3/s

      DO lay = 2,Layersin

	     EMI(lay,1)  = ERout(1,lay)  * ISO_EMI /  68 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)
	     
	     EMI(lay,2)  = ERout(2,lay)  * MBO_EMI /  86 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)

         EMI(lay,3)  = ERout(3,lay)  * MYR_EMI / 136 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)

         EMI(lay,4)  = ERout(4,lay)  * SAB_EMI / 136 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)

         EMI(lay,5)  = ERout(5,lay)  * LIM_EMI / 136 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)

         EMI(lay,6)  = ERout(6,lay)  * CAR_EMI / 136 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)

         EMI(lay,7)  = ERout(7,lay)  * OCI_EMI / 136 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)

         EMI(lay,8)  = ERout(8,lay)  * BPI_EMI / 136 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)

         EMI(lay,9)  = ERout(9,lay)  * API_EMI / 136 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)

         EMI(lay,10) = ERout(10,lay) * FAR_EMI / 204 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)

         EMI(lay,11) = ERout(11,lay) * BCA_EMI / 204 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)
        
         EMI(lay,12) = ERout(12,lay) * MET_EMI /  32 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)

         EMI(lay,13) = ERout(13,lay) * ACT_EMI /  58 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)
         
         EMI(lay,14) = ERout(14,lay) * ACA_EMI /  44 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)

         EMI(lay,15) = ERout(15,lay) * FOR_EMI /  30 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)

         EMI(lay,16) = ERout(16,lay) * CH4_EMI /  16 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)
         
         EMI(lay,17) = ERout(17,lay) * NO_EMI  /  30 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)
         
         EMI(lay,18)  = ERout(18,lay)  * OMT_EMI / 136 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)
         
         EMI(lay,19) = ERout(19,lay) * OSQ_EMI / 204 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)
        
         EMI(lay,20) = ERout(20,lay) * CO_EMI  /  28 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)

         EMI(lay,21) = 0.            * CIN_EMI / 154 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)

         EMI(lay,22) = 0.            * LIN_EMI / 154 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)

	  ENDDO

      EMI(1,1)  = ERout(1,1)  * ISO_EMI /  68 * Avog / Zin(1) / 100 * Ladpin(1)
      EMI(1,2) =  ERout(2,1) * MBO_EMI /  86 * Avog / Zin(1) / 100 * Ladpin(1)
      EMI(1,3)  = ERout(3,1)  * MYR_EMI / 136 * Avog / Zin(1) / 100 * Ladpin(1)
      EMI(1,4)  = ERout(4,1)  * SAB_EMI / 136 * Avog / Zin(1) / 100 * Ladpin(1)
      EMI(1,5)  = ERout(5,1)  * LIM_EMI / 136 * Avog / Zin(1) / 100 * Ladpin(1)
      EMI(1,6)  = ERout(6,1)  * CAR_EMI / 136 * Avog / Zin(1) / 100 * Ladpin(1)
      EMI(1,7)  = ERout(7,1)  * OCI_EMI / 136 * Avog / Zin(1) / 100 * Ladpin(1)
      EMI(1,8)  = ERout(8,1)  * BPI_EMI / 136 * Avog / Zin(1) / 100 * Ladpin(1)
      EMI(1,9)  = ERout(9,1)  * API_EMI / 136 * Avog / Zin(1) / 100 * Ladpin(1)
      EMI(1,10) = ERout(10,1) * FAR_EMI / 204 * Avog / Zin(1) / 100 * Ladpin(1)
      EMI(1,11) = ERout(11,1) * BCA_EMI / 204 * Avog / Zin(1) / 100 * Ladpin(1)
      EMI(1,12) = ERout(12,1) * MET_EMI /  32 * Avog / Zin(1) / 100 * Ladpin(1)
      EMI(1,13) = ERout(13,1) * ACT_EMI /  58 * Avog / Zin(1) / 100 * Ladpin(1)
      EMI(1,14) = ERout(14,1) * ACA_EMI /  44 * Avog / Zin(1) / 100 * Ladpin(1)
      EMI(1,15) = ERout(15,1) * FOR_EMI /  30 * Avog / Zin(1) / 100 * Ladpin(1)
      EMI(1,16) = ERout(16,1) * CH4_EMI /  16 * Avog / Zin(1) / 100 * Ladpin(1)
      EMI(1,17) = ERout(17,1) * NO_EMI  /  30 * Avog / Zin(1) / 100 * Ladpin(1)
      EMI(1,18) = ERout(18,1)  * OMT_EMI / 136 * Avog / Zin(1) / 100 * Ladpin(1)
      EMI(1,19) = ERout(19,1) * OSQ_EMI / 204 * Avog / Zin(1) / 100 * Ladpin(1)
      EMI(1,20) = ERout(20,1) * CO_EMI  /  28 * Avog / Zin(1) / 100 * Ladpin(1)
      EMI(1,21) = 0.          * CIN_EMI / 154 * Avog / Zin(1) / 100 * Ladpin(1)
      EMI(1,22) = 0.          * LIN_EMI / 154 * Avog / Zin(1) / 100 * Ladpin(1)
    

	  DO lay = Layersin+1, kzin

         sun_par(lay) = 1.

         DO J = 1,22

            EMI(lay,J) = 0.

         ENDDO

      ENDDO


      DEALLOCATE ( ER_BT      )   ! output emission buffer
      DEALLOCATE ( ER_NT      )   ! output emission buffer
      DEALLOCATE ( ER_SB      )   ! output emission buffer
      DEALLOCATE ( ER_Grass      )   ! output emission buffer

      DEALLOCATE ( EF      )   ! input annual emission factor
      DEALLOCATE ( LAT     )   ! input latitude of grid cell
      DEALLOCATE ( LONG    )   ! input longitude of grid cell

      DEALLOCATE ( LAIp    )   ! previous monthly LAI
      DEALLOCATE ( LAIc    )   ! current monthly LAI

      DEALLOCATE ( TEMP    )   ! input hourly temperature (K)
      DEALLOCATE ( SRAD    )   ! input hourly radiation
      DEALLOCATE ( PPFD    )   ! calculated PAR (umol/m2.s)

      DEALLOCATE ( D_SRAD  )   ! daily solar radiation (Watt/m2)
      DEALLOCATE ( D_PPFD  )   ! daily PAR (umol/m2.s)
      DEALLOCATE ( D_TEMP  )   ! input daily temperature (K)

      DEALLOCATE ( GAM_PHO )   ! light correction factor
      DEALLOCATE ( GAM_TMP )   ! temperature correction factor
      DEALLOCATE ( GAM_LHT )   ! LAI correction factor
      DEALLOCATE ( GAM_AGE )   ! leaf age correction factor
      DEALLOCATE ( GAM_SMT )   ! Soil moilture correction factor
      DEALLOCATE ( GAM_T )

      DEALLOCATE ( Wind )   
      DEALLOCATE ( Pres )   
      DEALLOCATE ( Humidity )
      DEALLOCATE ( DI )   
      !DEALLOCATE ( Cantype )
      DEALLOCATE ( GAM_CE )
      DEALLOCATE ( ISLTYP  )
      DEALLOCATE ( SMOIS )
      DEALLOCATE ( PFTF)
      DEALLOCATE ( GAMMAfactors)


!======================================================================
!..  FORMAT
!======================================================================
1000  FORMAT( A )
1010  FORMAT( 40( A, :, I8, :, 1X ) )

!======================================================================

!..  End program

!======================================================================



      END SUBROUTINE EMISSION_M2




!======================================================================

!..  Some added functionality, can be moved to a new module

!======================================================================

     
      Pure function Timeformatter(input) result(output)
      
         IMPLICIT NONE
      
         REAL, INTENT(IN) :: input
         INTEGER :: output
      
         REAL :: tmp
         INTEGER :: h, m ,s
      
         tmp= input
         output=0
      
         h   = INT(tmp/3600) ! truncation of values should work as intented
         tmp = tmp - (h*3600)
         m   = INT(tmp/60)
         tmp = tmp - (m*60)
         s   = INT(tmp)
      
         output = (10000*h) + (100*m) + s
      
      END function timeformatter


      
      Pure function centralizer(Inputheightarray, Inputlayers) result(output)
      
         IMPLICIT NONE
      
         REAL, DIMENSION(:), INTENT(IN) :: Inputheightarray
         INTEGER, INTENT(IN) :: Inputlayers
         REAL, DIMENSION(Inputlayers) :: output
            
         INTEGER :: ind
         REAL :: temp
           
           
         output= Inputheightarray(1:Inputlayers)

           
         DO ind=1, Inputlayers
            temp=0
            
	         IF(ind == 1) THEN
		         temp= Inputheightarray(ind)
	         ELSE
		         temp = Inputheightarray(ind) - Inputheightarray(ind-1)
	         END IF
	
	         output(ind) = output(ind) - (temp*0.5)
                  
        END DO
         
      END function centralizer      
      
END MODULE MEGAN_version_2
