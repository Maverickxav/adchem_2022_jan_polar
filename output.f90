MODULE output

   USE second_Precision,  ONLY : dp    ! KPP Numerical type
   USE constants ! Set case specifics (module also includes physical constants)
   USE second_Parameters, ONLY : NSPEC, ind_H2SO4, ind_HNO3

   IMPLICIT NONE

CONTAINS

   SUBROUTINE write_output(time,conc,index_cond,N_bins,V_bins,d_p,d_p_dry,c_p,v_dep,v_wet1,&
   dens_p,index_st,molar_mass,psat,Egases,kz,vd_gas,tr,PN_marine,&
   Jnucl_N,Jnucl_D,CS_H2SO4,RH,TEMP,pH)

    REAL(dp), INTENT(in) :: time
    REAL(dp), DIMENSION(Nz,NSPEC), INTENT(in) :: conc
    REAL(dp), DIMENSION(6,NSPEC), INTENT(in) :: Egases
    REAL(dp), DIMENSION(NCOND), INTENT(in) :: molar_mass,psat
    INTEGER, DIMENSION(NCOND-34), INTENT(in) :: index_cond
    REAL(dp), DIMENSION(Nz,nr_bins), INTENT(in) :: N_bins,V_bins,d_p,dens_p,d_p_dry,pH
    REAL(dp), DIMENSION(Nz,NSPEC_P,nr_bins), INTENT(in) :: c_p
    REAL(dp), DIMENSION(nr_bins), INTENT(in) :: v_dep, v_wet1,PN_marine
    REAL(dp), DIMENSION(Nz+1), INTENT(in) :: kz
    REAL(dp), DIMENSION(Nz), INTENT(in) :: Jnucl_N,Jnucl_D, CS_H2SO4,RH,TEMP
    REAL(dp), DIMENSION(12), INTENT(in) :: vd_gas
    INTEGER, INTENT(in) :: index_st
    INTEGER, INTENT(in) :: tr

    REAL(dp), DIMENSION(Nz,NSPEC_P) :: c_p_tot
    
    INTEGER  :: time_in_ms,i ! time in milliseconds
    INTEGER, PARAMETER :: one_hour = 60*60 ! in seconds
    INTEGER, PARAMETER :: ten_min = 10*60 ! in seconds
	LOGICAL, SAVE :: first_time = .TRUE.
	LOGICAL, SAVE :: first_time2 = .TRUE.
	

    c_p_tot = SUM(c_p,DIM=3) ! conc of each cond. comp at each height

    time_in_ms = FLOOR(1000*time)
    !IF (tr >= index_st-2880) THEN ! Start saving 2 hours upwind the station
    !IF (MODULO(time_in_ms, 1000*ten_min) == 0) THEN ! what a hack
    !WRITE(201,*) N_bins
    !END IF
    !END IF

    ! write output every full hour    
    ! IF (MODULO(time_in_ms, 1000*one_hour) == 0) THEN ! what a hack
       ! WRITE(*,'(A8,F6.3,A6)') 'time = ', time/(24*one_hour), '  days'   
	   ! !write(*,*) c_p_tot_tot(1,:)
	   ! write(*,*)sum(N_bins(1,:))*1D-6
	   ! write(*,*)sum(V_bins(1,:))*1D12
	   
	   ! !write(*,*)sum(N_bins((pi*d_p_dry**3D0)/6D0,DIM=2)*1D12
	   ! IF (tr >= index_st-288) THEN ! Start saving 2 day upwind the station
         ! WRITE(200,*) conc
         ! WRITE(202,*) V_bins
         ! WRITE(203,*) c_p_tot
         ! WRITE(204,*) c_p(1,:,:)
         ! WRITE(205,*) dens_p
         ! WRITE(206,*) v_dep
         ! WRITE(207,*) v_wet1
         ! WRITE(208,*) psat
         ! WRITE(209,*) SUM(Egases,DIM=1)
         ! WRITE(210,*) kz
         ! !WRITE(211,*) cNH3
         ! WRITE(219,*) d_p
         ! WRITE(221,*) d_p_dry
         ! WRITE(224,*) vd_gas
         ! !WRITE(225,*) ra_gas
         ! !WRITE(226,*) rb_gas
         ! !WRITE(227,*) rc_gas
         ! WRITE(228,*) PN_marine 
         ! !WRITE(229,*) elvoc_yield 
         ! WRITE(230,*) CS_H2SO4
         ! WRITE(231,*) Jnucl*1D-6
         ! !WRITE(231,*) Jnucl_org		 
         ! !WRITE(233,*) dimer_O
         ! !WRITE(234,*) dimer_N
         ! !WRITE(235,*) dimer_H
      ! END IF   
  
     ! END IF   
	
     ! write output when trajectory passes station
     !DO i = 1,nr_st
        IF (tr == index_st) THEN
		IF (first_time) THEN ! do only once
		WRITE(204,*) c_p
        WRITE(212,*) conc
        WRITE(213,*) N_bins
        !WRITE(214,*) V_bins
        WRITE(215,*) c_p_tot
        WRITE(216,*) dens_p
        WRITE(217,*) psat
        WRITE(222,*) d_p_dry
		WRITE(230,*) CS_H2SO4
		WRITE(231,*) Jnucl_N*1D-6
        WRITE(300,*) Jnucl_D*1D-6
		WRITE(232,*) RH
        WRITE(233,*) TEMP 
        WRITE(234,*) pH 		
		END IF
		first_time = .FALSE.
        END IF
     !END DO

     ! write only once
     IF (tr == 1) THEN
	 IF (first_time2) THEN ! do only once
        WRITE(218,*) index_cond
        WRITE(220,*) molar_mass ! kg/mol cond species
	 END IF
        first_time2 = .FALSE.	 
     END IF

     END SUBROUTINE write_output

END MODULE output     
