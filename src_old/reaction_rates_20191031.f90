MODULE reaction_rates_20191031

    USE second_Precision,  ONLY : dp ! KPP Numerical type
    USE second_Global, ONLY : NPHOT, NKVALUES
    USE constants

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: getKVALUES, getJVALUES

CONTAINS

    SUBROUTINE getKVALUES(KVALUES,TEMP,M,O2,H2O)
        ! Reaction rate coefficients from MCMv3.2 (http://mcm.leeds.ac.uk/MCM/roots.htt), 31 (the last are added
        ! in MALTE-BOX and right now put to 0.00001

        IMPLICIT NONE

        REAL(dp), INTENT(in) :: TEMP, M, O2, H2O ! Temperature, conc of intert molec, O2 and H2O in [molec/cm^3]
        REAL(dp), INTENT(out) :: KVALUES(NKVALUES)

        REAL(dp) :: F1, F10, F12, F13, F14, F15, F16, F17, F2, F3, F4, F7, &
            F8, F9, FC, FC1, FC10, FC12, FC13, FC14, FC15, FC16, FC17, FC2, FC3, &
            FC4, FC7, FC8, FC9, FCC, FCD, FD, K1, K10, K100, K10i, K120, K12i, &
            K130, K13i, K140, K14i, K150, K15i, K160, K16i, K170, K17i, K1i, K2, &
            K20, K2i, K3, K30, K3i, K4, K40, K4i, K70, K7i, K80, K8i, K90, K9i, &
            KC0, KCi, KD0, KDi, KR1, KR10, KR12, KR13, KR14, KR15, KR16, KR17, &
            KR2, KR3, KR4, KR7, KR8, KR9, KRC, KRD, NC, NC1, NC10, NC12, NC13, &
            NC14, NC15, NC16, NC17, NC2, NC3, NC4, NC7, NC8, NC9, NCD

        INTEGER :: i

        K10 =   1.0D-31*M*(TEMP/300)**(-1.6)
        K1I =   3.00D-11*(TEMP/300)**0.3
        KR1 =   K10/K1I
        FC1 =   0.85
        NC1 =   0.75-1.27*(log10(FC1))
        F1  =   10**(log10(FC1)/(1+(log10(KR1)/NC1)**2))
        KVALUES(1) = (K10*K1I)*F1/(K10+K1I) ! KMT01

        K20 =   1.3D-31*M*(TEMP/300)**(-1.5)
        K2I =   2.3D-11*(TEMP/300)**0.24
        KR2 =   K20/K2I
        FC2 =   0.6
        NC2 =   0.75-1.27*(log10(FC2))
        F2  =   10**(log10(FC2)/(1+(log10(KR2)/NC2)**2))
        KVALUES(2)   =   (K20*K2I)*F2/(K20+K2I) ! KMT02

        K30 =   3.6D-30*M*(TEMP/300)**(-4.1)
        K3I =   1.9D-12*(TEMP/300)**0.2
        KR3 =   K30/K3I
        FC3 =   0.35
        NC3 =   0.75-1.27*(log10(FC3))
        F3  =   10**(log10(FC3)/(1+(log10(KR3)/NC3)**2))
        KVALUES(3)   =   (K30*K3I)*F3/(K30+K3I) ! KMT03

        K40 =   1.3D-3*M*(TEMP/300)**(-3.5)*exp(-11000/TEMP)
        K4I =   9.7D+14*(TEMP/300)**0.1*exp(-11080/TEMP)
        KR4 =   K40/K4I
        FC4 =   0.35
        NC4 =   0.75-1.27*(log10(FC4))
        F4  =   10**(log10(FC4)/(1+(log10(KR4)/NC4)**2))
        KVALUES(4)   =   (K40*K4I)*F4/(K40+K4I) ! KMT04
        KVALUES(5)   =   1.44D-13*(1+(M/4.2D+19)) ! KMT05
        KVALUES(6)   =   1   +   (1.40D-21*exp(2200/TEMP)*H2O) ! KMT06

        K70 =   7.4D-31*M*(TEMP/300)**(-2.4)
        K7I =   3.3D-11*(TEMP/300)**(-0.3)
        KR7 =   K70/K7I
        FC7 =   exp(-TEMP/1420)
        NC7 =   0.75-1.27*(log10(FC7))
        F7  =   10**(log10(FC7)/(1+(log10(KR7)/NC7)**2))
        KVALUES(7)   =   (K70*K7I)*F7/(K70+K7I) ! KMT07

        K80 =   3.3D-30*M*(TEMP/300)**(-3.0)
        K8I =   4.1D-11
        KR8 =   K80/K8I
        FC8 =   0.4
        NC8 =   0.75-1.27*(log10(FC8))
        F8  =   10**(log10(FC8)/(1+(log10(KR8)/NC8)**2))
        KVALUES(8)   =   (K80*K8I)*F8/(K80+K8I) ! KMT08

        K90 =   1.8D-31*M*(TEMP/300)**(-3.2)
        K9I =   4.7D-12
        KR9 =   K90/K9I
        FC9 =   0.6
        NC9 =   0.75-1.27*(log10(FC9))
        F9  =   10**(log10(FC9)/(1+(log10(KR9)/NC9)**2))
        KVALUES(9)   =   (K90*K9I)*F9/(K90+K9I) ! KMT09

        K100    =   4.10D-05*M*exp(-10650/TEMP)
        K10I    =   4.8D+15*exp(-11170/TEMP)
        KR10    =   K100/K10I
        FC10    =   0.6
        NC10    =   0.75-1.27*(log10(FC10))
        F10 =   10**(log10(FC10)/(1+(log10(KR10)/NC10)**2))
        KVALUES(10)   =   (K100*K10I)*F10/(K100+K10I) ! KMT10

        K1  =   2.40D-14*exp(460/TEMP)
        K3  =   6.50D-34*exp(1335/TEMP)
        K4  =   2.70D-17*exp(2199/TEMP)
        K2  =   (K3*M)/(1+(K3*M/K4))
        KVALUES(11)   =   K1  +   K2 ! KMT11

        K120    =   4.5D-31*M*(TEMP/300)**(-3.9)
        K12I    =   1.3D-12*(TEMP/300)**(-0.7)
        KR12    =   K120/K12I
        FC12    =   0.525
        NC12    =   0.75-1.27*(log10(FC12))
        F12 =   10**(log10(FC12)/(1.0+(log10(KR12)/NC12)**2))
        KVALUES(12)   =   (K120*K12I*F12)/(K120+K12I) ! KMT12

        K130    =   2.5D-30*M*(TEMP/300)**(-5.5)
        K13I    =   1.8D-11
        KR13    =   K130/K13I
        FC13    =   0.36
        NC13    =   0.75-1.27*(log10(FC13))
        F13 =   10**(log10(FC13)/(1+(log10(KR13)/NC13)**2))
        KVALUES(13)   =   (K130*K13I)*F13/(K130+K13I) ! KMT13

        K140    =   9.0D-5*exp(-9690/TEMP)*M
        K14I    =   1.1D+16*exp(-10560/TEMP)
        KR14    =   K140/K14I
        FC14    =   0.4
        NC14    =   0.75-1.27*(log10(FC14))
        F14 =   10**(log10(FC14)/(1+(log10(KR14)/NC14)**2))
        KVALUES(14)   =   (K140*K14I)*F14/(K140+K14I) ! KMT14

        K150    =   8.6D-29*M*(TEMP/300)**(-3.1)
        K15I    =   9.0D-12*(TEMP/300)**(-0.85)
        KR15    =   K150/K15I
        FC15    =   0.48
        NC15    =   0.75-1.27*(log10(FC15))
        F15 =   10**(log10(FC15)/(1+(log10(KR15)/NC15)**2))
        KVALUES(15)   =   (K150*K15I)*F15/(K150+K15I) ! KMT15

        K160    =   8D-27*M*(TEMP/300)**(-3.5)
        K16I    =   3.0D-11*(TEMP/300)**(-1)
        KR16    =   K160/K16I
        FC16    =   0.5
        NC16    =   0.75-1.27*(log10(FC16))
        F16 =   10**(log10(FC16)/(1+(log10(KR16)/NC16)**2))
        KVALUES(16)   =   (K160*K16I)*F16/(K160+K16I) ! KMT16

        K170    =   5.0D-30*M*(TEMP/300)**(-1.5)
        K17I    =   1.0D-12
        KR17    =   K170/K17I
        FC17    =   0.17*exp(-51/TEMP)+exp(-TEMP/204)
        NC17    =   0.75-1.27*(log10(FC17))
        F17 =   10**(log10(FC17)/(1.0+(log10(KR17)/NC17)**2))
        KVALUES(17)   =   (K170*K17I*F17)/(K170+K17I) ! KMT17
		
        KVALUES(18)   =   9.5D-39*O2*exp(5270/TEMP)/(1+7.5D-29*O2*exp(5610/TEMP)) ! KMT18
		
		
        !KC0 =   2.7D-28*M*(TEMP/300)**(-7.1)
        !KCI =   1.2D-11*(TEMP/300)**(-0.9)
        !KRC =   KC0/KCI
        !FCC =   0.30
        !NC  =   0.75-1.27*(log10(FCC))
        !FC  =   10**(log10(FCC)/(1+(log10(KRC)/NC)**2))
        !KVALUES(19)   =   (KC0*KCI)*FC/(KC0+KCI) ! KFPAN
        
        KC0	= 3.28D-28*M*(TEMP/300)**(-6.87)
        KCI = 1.125D-11*(TEMP/300)**(-1.105)
        KRC = KC0/KCI
        FCC = 0.30
        NC  = 0.75-1.27*(LOG10(FCC))
        FC  = 10**(LOG10(FCC)/(1+(LOG10(KRC)/NC)**2))	
        KVALUES(19)   =   (KC0*KCI)*FC/(KC0+KCI) ! KFPAN


        KD0 =   4.90D-3*exp(-12100/TEMP)*M
        KDI =   5.4D+16*exp(-13830/TEMP)
        KRD =   KD0/KDI
        FCD =   0.30
        NCD =   0.75-1.27*(log10(FCD))
        FD  =   10**(log10(FCD)/(1+(log10(KRD)/NCD)**2))
        KVALUES(20)   =   (KD0*KDI)*FD/(KD0+KDI) ! KBPAN

        KVALUES(21)  =   2.7D-12*exp(360/TEMP) ! KRO2NO
        KVALUES(22) =   2.91D-13*exp(1300/TEMP) ! KRO2HO2
        KVALUES(23)  =   5.2D-13*exp(980/TEMP) ! KAPHO2
        KVALUES(24)   =   7.5D-12*exp(290/TEMP) ! KAPNO
        KVALUES(25) =   2.3D-12 ! KRO2NO3
        KVALUES(26)  =   1.4D-12*exp(-1860/TEMP) ! KNO3AL
        KVALUES(27)    =   1.00D+06 ! KDEC
        KVALUES(28) =   2.50D-14*exp(-300/TEMP) ! KROPRIM
        KVALUES(29)  =   2.50D-14*exp(-300/TEMP) ! KROSEC
        KVALUES(30)  =   1.03D-13*exp(365/TEMP) ! KCH3O2
        KVALUES(31)   =   3.5D-13 ! K298CH3O2

        ! KVALUES(32:42) not in MCMv3.2
        DO i = 32, NKVALUES
            KVALUES(i) = 0.00001
        END DO

        ! KVALUES new in MCMv3.3
        KVALUES(43) = 4.6D9*exp(-8380./TEMP)*exp(1D8/TEMP**3.)
        KVALUES(44) = 1.5D11*exp(-9750./TEMP)
        KVALUES(45) = 3D7*exp(-5300./TEMP)
		
		! CLO + CLO +M -> CL2O2
		K10 =   2D-32*(M-O2)*(TEMP/300)**(-4.0)
        K1I =   1.00D-11
        FC1 =   0.45
        !KVALUES(46) = K10/(1.0+K10/K1I)*FC1**((1.0+log10(K10/K1I))**(-2.0)) 
		KVALUES(46) = K1I*((K10/K1I)/(1.0+K10/K1I))*(10**(log10(FC1)/(1+(log10(K10/K1I)/(0.75-1.27*log10(FC1)))**2.0))) !K10/(1.0+K10/K1I)*FC1**((1.0+log10(K10/K1I))**(-2.0)) 
		
		! CL + O2 -> CLO2
		K10 =   2.2D-33*M*(TEMP/300)**(-3.1)
        K1I =   1.8D-10
        FC1 =   0.6
        !KVALUES(47) = K10/(1.0+K10/K1I)*FC1**((1.0+log10(K10/K1I))**(-2.0))
        KVALUES(47) = K1I*((K10/K1I)/(1.0+K10/K1I))*(10**(log10(FC1)/(1+(log10(K10/K1I)/(0.75-1.27*log10(FC1)))**2.0))) !K10/(1.0+K10/K1I)*FC1**((1.0+log10(K10/K1I))**(-2.0)) 
				

		! CL2O2 -> CLO + CLO
		K10 =   3.7D-7*(M-O2)*EXP(-7690/TEMP)
        K1I =   7.9D15*EXP(-8820/TEMP)
        FC1 =   0.45
        !KVALUES(48) = K10/(1.0+K10/K1I)*FC1**((1.0+log10(K10/K1I))**(-2.0)) 
		KVALUES(48) = K1I*((K10/K1I)/(1.0+K10/K1I))*(10**(log10(FC1)/(1+(log10(K10/K1I)/(0.75-1.27*log10(FC1)))**2.0))) !K10/(1.0+K10/K1I)*FC1**((1.0+log10(K10/K1I))**(-2.0)) 
		
		! CLO + OCLO -> CL2O3
		K10 =   6.2D-32*(M-O2)*(TEMP/300)**(-4.7)
        K1I =   2.4D-11
        FC1 =   0.6
        !KVALUES(49) = K10/(1.0+K10/K1I)*FC1**((1.0+log10(K10/K1I))**(-2.0)) 
		KVALUES(49) = K1I*((K10/K1I)/(1.0+K10/K1I))*(10**(log10(FC1)/(1+(log10(K10/K1I)/(0.75-1.27*log10(FC1)))**2.0))) !K10/(1.0+K10/K1I)*FC1**((1.0+log10(K10/K1I))**(-2.0)) 
		
		! CL2O3 -> CLO + OCLO
		K10 =   1.4D-10*(M-O2)*EXP(-3810/TEMP)
        K1I =   2.5D12*EXP(-4940/TEMP)
        FC1 =   0.6
        !KVALUES(50) = K10/(1.0+K10/K1I)*FC1**((1.0+log10(K10/K1I))**(-2.0)) 
		KVALUES(50) = K1I*((K10/K1I)/(1.0+K10/K1I))*(10**(log10(FC1)/(1+(log10(K10/K1I)/(0.75-1.27*log10(FC1)))**2.0))) !K10/(1.0+K10/K1I)*FC1**((1.0+log10(K10/K1I))**(-2.0)) 
		
        ! CL +NO -> CLNO
		KVALUES(51)=7.6D-32*M*(TEMP/300)**(-1.8)
		
		! CL +NO2 -> CLNO2
		K10 =   1.8D-31*M*(TEMP/300)**(-2.0)
        K1I =   1D-10*(TEMP/300)**(-1.0)
        FC1 =   0.6
        !KVALUES(52) = K10/(1.0+K10/K1I)*FC1**((1.0+log10(K10/K1I))**(-2.0)) 
		KVALUES(52) = K1I*((K10/K1I)/(1.0+K10/K1I))*(10**(log10(FC1)/(1+(log10(K10/K1I)/(0.75-1.27*log10(FC1)))**2.0))) !K10/(1.0+K10/K1I)*FC1**((1.0+log10(K10/K1I))**(-2.0)) 
		
		! CLO + NO2 -> CLNO3
		K10 =   1.6D-31*(M-O2)*(TEMP/300)**(-3.4)
        K1I =   7D-11
        FC1 =   0.4
        KVALUES(53) = K1I*((K10/K1I)/(1.0+K10/K1I))*(10**(log10(FC1)/(1+(log10(K10/K1I)/(0.75-1.27*log10(FC1)))**2.0))) !K10/(1.0+K10/K1I)*FC1**((1.0+log10(K10/K1I))**(-2.0)) 
		
!		DMS + Br -> CH3SCH3Br
		K10 =   3.7D-29*(M)*(TEMP/300)**(-5.3)
        K1I =   1.5D-10*(TEMP/300)**(-2.0)
        FC1 =   0.6
        KVALUES(54) = K1I*((K10/K1I)/(1.0+K10/K1I))*(10**(log10(FC1)/(1+(log10(K10/K1I)/(0.75-1.27*log10(FC1)))**2.0))) !K10/(1.0+K10/K1I)*FC1**((1.0+log10(K10/K1I))**(-2.0)) 
		
!		Br + NO2 -> BrNO2 : ; // G149 3-body reaction Atkinson et al., 2007
        K10 =   4.2D-31*(M-O2)*(TEMP/300)**(-2.4)
        K1I =   2.8D-11
        FC1 =   0.55
        KVALUES(55) = K1I*((K10/K1I)/(1.0+K10/K1I))*(10**(log10(FC1)/(1+(log10(K10/K1I)/(0.75-1.27*log10(FC1)))**2.0))) !K10/(1.0+K10/K1I)*FC1**((1.0+log10(K10/K1I))**(-2.0)) 
		
!       BrO + NO2 -> BrNO3 : ; // G152 3-body reaction Atkinson et al., 2007 
        K10 =   4.7D-31*(M-O2)*(TEMP/300)**(-3.1)
        K1I =   1.8D-11
        FC1 =   0.4
        KVALUES(56) = K1I*((K10/K1I)/(1.0+K10/K1I))*(10**(log10(FC1)/(1+(log10(K10/K1I)/(0.75-1.27*log10(FC1)))**2.0))) !K10/(1.0+K10/K1I)*FC1**((1.0+log10(K10/K1I))**(-2.0)) 


!       I + NO -> INO : KMT57; // G239 Atkinson et al., 2007
        K10 =   1.8D-32*(M-O2)*(TEMP/300)**(-1.0)
        K1I =   1.7D-11
        FC1 =   0.6
        KVALUES(57) = K1I*((K10/K1I)/(1.0+K10/K1I))*(10**(log10(FC1)/(1+(log10(K10/K1I)/(0.75-1.27*log10(FC1)))**2.0))) !K10/(1.0+K10/K1I)*FC1**((1.0+log10(K10/K1I))**(-2.0)) 

!       I + NO2 -> INO2 : KMT58; // G240 Atkinson et al., 2007
        K10 =   3D-31*(M-O2)*(TEMP/300)**(-1.0)
        K1I =   6.6D-11
        FC1 =   0.63
        KVALUES(58) = K1I*((K10/K1I)/(1.0+K10/K1I))*(10**(log10(FC1)/(1+(log10(K10/K1I)/(0.75-1.27*log10(FC1)))**2.0))) !K10/(1.0+K10/K1I)*FC1**((1.0+log10(K10/K1I))**(-2.0)) 

! IO + NO2 -> INO3 : KMT59; // G244 Atkinson et al., 2007
        K10 =   7.7D-31*(M-O2)*(TEMP/300)**(-5.0)
        K1I =   1.6D-11
        FC1 =   0.4
        KVALUES(59) = K1I*((K10/K1I)/(1.0+K10/K1I))*(10**(log10(FC1)/(1+(log10(K10/K1I)/(0.75-1.27*log10(FC1)))**2.0))) !K10/(1.0+K10/K1I)*FC1**((1.0+log10(K10/K1I))**(-2.0)) 

! CCL3O2 + NO2 -> CCL3O2NO2 : KMT60; // G203 Atkinson et al., 2008
        K10 =   9.2D-29*(M-O2)*(TEMP/298)**(-6.0)
        K1I =   1.5D-12*(TEMP/298)**(-0.7)
        FC1 =   0.32
       KVALUES(60) = K1I*((K10/K1I)/(1.0+K10/K1I))*(10**(log10(FC1)/(1+(log10(K10/K1I)/(0.75-1.27*log10(FC1)))**2.0))) !K10/(1.0+K10/K1I)*FC1**((1.0+log10(K10/K1I))**(-2.0)) 

! CCL3O2NO2 -> CCL3O2 + NO2 : KMT61; // G204 Atkinson et al., 2008
        K10 =   4.3D-3*(M-O2)*(-10235/TEMP)
        K1I =   4.8D16*(-11820/TEMP)
        FC1 =   0.32
       KVALUES(61) = K1I*((K10/K1I)/(1.0+K10/K1I))*(10**(log10(FC1)/(1+(log10(K10/K1I)/(0.75-1.27*log10(FC1)))**2.0))) !K10/(1.0+K10/K1I)*FC1**((1.0+log10(K10/K1I))**(-2.0)) 
		
		!write(*,*) KVALUES(46:53)
		
		!stop

    END SUBROUTINE getKVALUES

    SUBROUTINE getJVALUES(JVALUES,TEMP,i_surface)
        ! Photolysis rate coefficients from MCMv3.2 (http://mcm.leeds.ac.uk/MCM/parameters/photolysis.htt)

        IMPLICIT NONE

        REAL(dp), INTENT(out)      :: JVALUES(NPHOT)
        REAL(dp), INTENT(in)       :: TEMP
        REAL(dp), DIMENSION(510), INTENT(in) :: i_surface

        ! Absorption cross section and quantum yield data used to calculate the photolysis rates in MCMv3.2
        REAL(dp), DIMENSION(510), SAVE :: wave_len, cs_1, cs_2, cs_3, cs_4, cs_5, cs_6, cs_7, cs_8, cs_11, &
            cs_12, cs_13, cs_14, cs_15, cs_16, cs_17, cs_18, cs_19, cs_21, cs_22, cs_23, cs_24, cs_31, cs_32, &
            cs_33, cs_34, cs_35, cs_41, cs_51, cs_52, cs_53, cs_54, cs_55, cs_56, cs_57, cs_61, cs_62, cs_63,&
			cs_64, cs_65, cs_66, cs_67, cs_68, cs_69, cs_70, cs_71, cs_72, cs_73, cs_74, cs_75, cs_76, cs_77, &
			cs_78, cs_79, cs_80, cs_81, cs_82, cs_83, cs_84, cs_85, cs_86, cs_87, cs_88, cs_89, cs_90, cs_91, &
			cs_92, cs_93, cs_94, cs_95, cs_96, &
			qy_1, qy_2, qy_3, qy_4, qy_5, qy_6, qy_7, qy_8, qy_11, qy_12, qy_13, qy_14, qy_15, qy_16, qy_17, &
			qy_18, qy_19, qy_21, qy_22, qy_23, qy_24, qy_31, qy_32, qy_33, qy_34, qy_35, qy_41, qy_51, qy_52, &
			qy_53, qy_54, qy_55, qy_56, qy_57, qy_61, qy_62, qy_63, qy_64, qy_65, qy_66, qy_67, qy_68, qy_69,&
			qy_70, qy_71, qy_72, qy_73, qy_74, qy_75, qy_76, qy_77, qy_78, qy_79, qy_80, qy_81, qy_82, qy_83, qy_84, &
			qy_85, qy_86, qy_87, qy_88, qy_89, qy_90, qy_91, qy_92, qy_93, qy_94, qy_95, qy_96, &
			cs_int_1, cs_int_2, cs_int_8, cs_int_11, cs_int_12, cs_int_21, cs_int_51, cs_int_52

        REAL(dp) :: lambda1(451), lambda(511), lambda2(20), lambda3(29), lambda4(11), dlambda(510), &
            cs_rawT_1(60), &
            cs_rawT_2(60), cs_rawT_8(83), cs_rawT_11(150), cs_rawT_12(150), cs_rawT_21(135), cs_rawT_51(21), &
            cs_rawT_52(22), cs_rawT_54(37), cs_int_54(110)
            
        REAL(dp), SAVE :: wave_len_1(60), wave_len_2(60), wave_len_8(83), wave_len_11(150), &
            wave_len_12(150), wave_len_21(135), wave_len_51(21), wave_len_52(22), wave_len_54(37),cs_raw_1(60), &
            cs_raw_2(60), cs_raw_8(83), cs_raw_11(150), &
            cs_raw_12(150), cs_raw_21(135), cs_raw_51(21), cs_raw_52(22), cs_raw_54(37), cs_T_1(60), &
            cs_T_2(60), cs_T_8(83), cs_T_11(150), cs_T_12(150), cs_T1_21(135), cs_T2_21(135), cs_T3_21(135), &
            cs_T_51(21), cs_T_52(22), cs_T_54(37)    
            
       

        LOGICAL, SAVE :: first_call = .TRUE.
        INTEGER       :: i, j, k


        ! The following is run only once, only on the first time chemistry is called
        ! Input of absorption cross spectrum and quantum yields for different molecules,
        ! temperature independent data are interpolated to 1 nm resolution in matlab-script
        ! numbers represent the respective photolysis rate in MCM
        IF (first_call) THEN ! do only once
            first_call = .FALSE.
            ! Temperature independent absorption cross sections
            OPEN(12, FILE='input/photolysis/cs_int_3.dat',STATUS='OLD')
            OPEN(13, FILE='input/photolysis/cs_int_4.dat',STATUS='OLD')
            OPEN(14, FILE='input/photolysis/cs_int_5.dat',STATUS='OLD')
            OPEN(15, FILE='input/photolysis/cs_int_6.dat',STATUS='OLD')
            OPEN(16, FILE='input/photolysis/cs_int_7.dat',STATUS='OLD')
            OPEN(17, FILE='input/photolysis/cs_int_13.dat',STATUS='OLD')
            OPEN(18, FILE='input/photolysis/cs_int_14.dat',STATUS='OLD')
            OPEN(19, FILE='input/photolysis/cs_int_15.dat',STATUS='OLD')
            OPEN(20, FILE='input/photolysis/cs_int_16.dat',STATUS='OLD')
            OPEN(21, FILE='input/photolysis/cs_int_17.dat',STATUS='OLD')
            OPEN(22, FILE='input/photolysis/cs_int_18.dat',STATUS='OLD')
            OPEN(23, FILE='input/photolysis/cs_int_19.dat',STATUS='OLD')
            OPEN(24, FILE='input/photolysis/cs_int_22.dat',STATUS='OLD')
            OPEN(25, FILE='input/photolysis/cs_int_23.dat',STATUS='OLD')
            OPEN(26, FILE='input/photolysis/cs_int_24.dat',STATUS='OLD')
            OPEN(27, FILE='input/photolysis/cs_int_31.dat',STATUS='OLD')
            OPEN(28, FILE='input/photolysis/cs_int_32.dat',STATUS='OLD')
            OPEN(29, FILE='input/photolysis/cs_int_33.dat',STATUS='OLD')
            OPEN(30, FILE='input/photolysis/cs_int_34.dat',STATUS='OLD')
            OPEN(31, FILE='input/photolysis/cs_int_35.dat',STATUS='OLD')
            OPEN(32, FILE='input/photolysis/cs_int_41.dat',STATUS='OLD')
            OPEN(33, FILE='input/photolysis/cs_int_53.dat',STATUS='OLD')
            OPEN(34, FILE='input/photolysis/cs_int_55.dat',STATUS='OLD')
            OPEN(35, FILE='input/photolysis/cs_int_56.dat',STATUS='OLD')
            OPEN(36, FILE='input/photolysis/cs_int_57.dat',STATUS='OLD')
			
			OPEN(91, FILE='input/photolysis/cs_int_61.dat',STATUS='OLD')
			OPEN(92, FILE='input/photolysis/cs_int_62.dat',STATUS='OLD')
			OPEN(93, FILE='input/photolysis/cs_int_63.dat',STATUS='OLD')
			OPEN(94, FILE='input/photolysis/cs_int_64.dat',STATUS='OLD')
			OPEN(95, FILE='input/photolysis/cs_int_65.dat',STATUS='OLD')
			OPEN(96, FILE='input/photolysis/cs_int_66.dat',STATUS='OLD')
			OPEN(97, FILE='input/photolysis/cs_int_67.dat',STATUS='OLD')
			OPEN(98, FILE='input/photolysis/cs_int_68.dat',STATUS='OLD')
			OPEN(99, FILE='input/photolysis/cs_int_69.dat',STATUS='OLD')
			OPEN(100, FILE='input/photolysis/cs_int_70.dat',STATUS='OLD')
			OPEN(101, FILE='input/photolysis/cs_int_71.dat',STATUS='OLD')
			OPEN(102, FILE='input/photolysis/cs_int_72.dat',STATUS='OLD')
			OPEN(103, FILE='input/photolysis/cs_int_73.dat',STATUS='OLD')
			OPEN(104, FILE='input/photolysis/cs_int_74.dat',STATUS='OLD')
			OPEN(105, FILE='input/photolysis/cs_int_75.dat',STATUS='OLD')
			OPEN(106, FILE='input/photolysis/cs_int_76.dat',STATUS='OLD')
			OPEN(107, FILE='input/photolysis/cs_int_77.dat',STATUS='OLD')
			OPEN(108, FILE='input/photolysis/cs_int_78.dat',STATUS='OLD')
			OPEN(109, FILE='input/photolysis/cs_int_79.dat',STATUS='OLD')
			
			OPEN(118, FILE='input/photolysis/cs_int_80.dat',STATUS='OLD')
			OPEN(119, FILE='input/photolysis/cs_int_81.dat',STATUS='OLD')
			OPEN(120, FILE='input/photolysis/cs_int_82.dat',STATUS='OLD')
			OPEN(121, FILE='input/photolysis/cs_int_83.dat',STATUS='OLD')
			OPEN(122, FILE='input/photolysis/cs_int_84.dat',STATUS='OLD')
			OPEN(123, FILE='input/photolysis/cs_int_85.dat',STATUS='OLD')
			OPEN(124, FILE='input/photolysis/cs_int_86.dat',STATUS='OLD')
			OPEN(125, FILE='input/photolysis/cs_int_87.dat',STATUS='OLD')
			OPEN(126, FILE='input/photolysis/cs_int_88.dat',STATUS='OLD')
			OPEN(127, FILE='input/photolysis/cs_int_89.dat',STATUS='OLD')
			OPEN(128, FILE='input/photolysis/cs_int_90.dat',STATUS='OLD')
			OPEN(129, FILE='input/photolysis/cs_int_91.dat',STATUS='OLD')
			OPEN(130, FILE='input/photolysis/cs_int_92.dat',STATUS='OLD')
			
			OPEN(144, FILE='input/photolysis/cs_int_93.dat',STATUS='OLD')
			OPEN(145, FILE='input/photolysis/cs_int_94.dat',STATUS='OLD')
			OPEN(146, FILE='input/photolysis/cs_int_95.dat',STATUS='OLD')
			OPEN(147, FILE='input/photolysis/cs_int_96.dat',STATUS='OLD')
			
            ! Temperature independent quantum yields:
            OPEN(37, FILE='input/photolysis/qy_int_1.dat',STATUS='OLD')
            OPEN(38, FILE='input/photolysis/qy_int_2.dat',STATUS='OLD')
            OPEN(39, FILE='input/photolysis/qy_int_3.dat',STATUS='OLD')
            OPEN(40, FILE='input/photolysis/qy_int_4.dat',STATUS='OLD')
            OPEN(41, FILE='input/photolysis/qy_int_5.dat',STATUS='OLD')
            OPEN(42, FILE='input/photolysis/qy_int_6.dat',STATUS='OLD')
            OPEN(43, FILE='input/photolysis/qy_int_7.dat',STATUS='OLD')
            OPEN(44, FILE='input/photolysis/qy_int_8.dat',STATUS='OLD')
            OPEN(45, FILE='input/photolysis/qy_int_11.dat',STATUS='OLD')
            OPEN(46, FILE='input/photolysis/qy_int_12.dat',STATUS='OLD')
            OPEN(47, FILE='input/photolysis/qy_int_13.dat',STATUS='OLD')
            OPEN(48, FILE='input/photolysis/qy_int_14.dat',STATUS='OLD')
            OPEN(49, FILE='input/photolysis/qy_int_15.dat',STATUS='OLD')
            OPEN(50, FILE='input/photolysis/qy_int_16.dat',STATUS='OLD')
            OPEN(51, FILE='input/photolysis/qy_int_17.dat',STATUS='OLD')
            OPEN(52, FILE='input/photolysis/qy_int_18.dat',STATUS='OLD')
            OPEN(53, FILE='input/photolysis/qy_int_19.dat',STATUS='OLD')
            OPEN(54, FILE='input/photolysis/qy_int_21.dat',STATUS='OLD')
            OPEN(55, FILE='input/photolysis/qy_int_22.dat',STATUS='OLD')
            OPEN(56, FILE='input/photolysis/qy_int_23.dat',STATUS='OLD')
            OPEN(57, FILE='input/photolysis/qy_int_24.dat',STATUS='OLD')
            OPEN(58, FILE='input/photolysis/qy_int_31.dat',STATUS='OLD')
            OPEN(59, FILE='input/photolysis/qy_int_32.dat',STATUS='OLD')
            OPEN(60, FILE='input/photolysis/qy_int_33.dat',STATUS='OLD')
            OPEN(61, FILE='input/photolysis/qy_int_34.dat',STATUS='OLD')
            OPEN(62, FILE='input/photolysis/qy_int_35.dat',STATUS='OLD')
            OPEN(63, FILE='input/photolysis/qy_int_41.dat',STATUS='OLD')
            OPEN(64, FILE='input/photolysis/qy_int_51.dat',STATUS='OLD')
            OPEN(65, FILE='input/photolysis/qy_int_52.dat',STATUS='OLD')
            OPEN(66, FILE='input/photolysis/qy_int_53.dat',STATUS='OLD')
            OPEN(67, FILE='input/photolysis/qy_int_54.dat',STATUS='OLD')
            OPEN(68, FILE='input/photolysis/qy_int_55.dat',STATUS='OLD')
            OPEN(69, FILE='input/photolysis/qy_int_56.dat',STATUS='OLD')
            OPEN(70, FILE='input/photolysis/qy_int_57.dat',STATUS='OLD')
			
			OPEN(80, FILE='input/photolysis/qy_int_61.dat',STATUS='OLD')
			OPEN(81, FILE='input/photolysis/qy_int_62.dat',STATUS='OLD')
			OPEN(82, FILE='input/photolysis/qy_int_63.dat',STATUS='OLD')
			OPEN(83, FILE='input/photolysis/qy_int_64.dat',STATUS='OLD')
			OPEN(84, FILE='input/photolysis/qy_int_65.dat',STATUS='OLD')
			OPEN(85, FILE='input/photolysis/qy_int_66.dat',STATUS='OLD')
			OPEN(86, FILE='input/photolysis/qy_int_67.dat',STATUS='OLD')
			OPEN(87, FILE='input/photolysis/qy_int_68.dat',STATUS='OLD')
			OPEN(88, FILE='input/photolysis/qy_int_69.dat',STATUS='OLD')
			OPEN(89, FILE='input/photolysis/qy_int_70.dat',STATUS='OLD')
			OPEN(90, FILE='input/photolysis/qy_int_71.dat',STATUS='OLD')
			OPEN(110, FILE='input/photolysis/qy_int_72.dat',STATUS='OLD')
			OPEN(111, FILE='input/photolysis/qy_int_73.dat',STATUS='OLD')
			OPEN(112, FILE='input/photolysis/qy_int_74.dat',STATUS='OLD')
			OPEN(113, FILE='input/photolysis/qy_int_75.dat',STATUS='OLD')
			OPEN(114, FILE='input/photolysis/qy_int_76.dat',STATUS='OLD')
			OPEN(115, FILE='input/photolysis/qy_int_77.dat',STATUS='OLD')
			OPEN(116, FILE='input/photolysis/qy_int_78.dat',STATUS='OLD')
			OPEN(117, FILE='input/photolysis/qy_int_79.dat',STATUS='OLD')

            OPEN(131, FILE='input/photolysis/qy_int_80.dat',STATUS='OLD')
            OPEN(132, FILE='input/photolysis/qy_int_81.dat',STATUS='OLD')
			OPEN(133, FILE='input/photolysis/qy_int_82.dat',STATUS='OLD')
			OPEN(134, FILE='input/photolysis/qy_int_83.dat',STATUS='OLD')
			OPEN(135, FILE='input/photolysis/qy_int_84.dat',STATUS='OLD')
			OPEN(136, FILE='input/photolysis/qy_int_85.dat',STATUS='OLD')
			OPEN(137, FILE='input/photolysis/qy_int_86.dat',STATUS='OLD')
			OPEN(138, FILE='input/photolysis/qy_int_87.dat',STATUS='OLD')
			OPEN(139, FILE='input/photolysis/qy_int_88.dat',STATUS='OLD')
			OPEN(140, FILE='input/photolysis/qy_int_89.dat',STATUS='OLD')
			OPEN(141, FILE='input/photolysis/qy_int_90.dat',STATUS='OLD')
			OPEN(142, FILE='input/photolysis/qy_int_91.dat',STATUS='OLD')
			OPEN(143, FILE='input/photolysis/qy_int_92.dat',STATUS='OLD')
			
			OPEN(148, FILE='input/photolysis/qy_int_93.dat',STATUS='OLD')
			OPEN(149, FILE='input/photolysis/qy_int_94.dat',STATUS='OLD')
			OPEN(150, FILE='input/photolysis/qy_int_95.dat',STATUS='OLD')
			OPEN(151, FILE='input/photolysis/qy_int_96.dat',STATUS='OLD')
			
            ! Temperature dependent absorption cross sections
            OPEN(71, FILE='input/photolysis/cs_1.dat',STATUS='OLD')
            OPEN(72, FILE='input/photolysis/cs_2.dat',STATUS='OLD')
            OPEN(73, FILE='input/photolysis/cs_8.dat',STATUS='OLD')
            OPEN(74, FILE='input/photolysis/cs_11.dat',STATUS='OLD')
            OPEN(75, FILE='input/photolysis/cs_12.dat',STATUS='OLD')
            OPEN(76, FILE='input/photolysis/cs_21.dat',STATUS='OLD')
            OPEN(77, FILE='input/photolysis/cs_51.dat',STATUS='OLD')
            OPEN(78, FILE='input/photolysis/cs_52.dat',STATUS='OLD')
            OPEN(79, FILE='input/photolysis/cs_54.dat',STATUS='OLD')

            DO i = 1,510 ! length of wavelength-array
                ! wavelength [nm], Absorption cross section
                READ(12,*) wave_len(i), cs_3(i)
                READ(13,*) wave_len(i), cs_4(i)
                READ(14,*) wave_len(i), cs_5(i)
                READ(15,*) wave_len(i), cs_6(i)
                READ(16,*) wave_len(i), cs_7(i)
                READ(17,*) wave_len(i), cs_13(i)
                READ(18,*) wave_len(i), cs_14(i)
                READ(19,*) wave_len(i), cs_15(i)
                READ(20,*) wave_len(i), cs_16(i)
                READ(21,*) wave_len(i), cs_17(i)
                READ(22,*) wave_len(i), cs_18(i)
                READ(23,*) wave_len(i), cs_19(i)
                READ(24,*) wave_len(i), cs_22(i)
                READ(25,*) wave_len(i), cs_23(i)
                READ(26,*) wave_len(i), cs_24(i)
                READ(27,*) wave_len(i), cs_31(i)
                READ(28,*) wave_len(i), cs_32(i)
                READ(29,*) wave_len(i), cs_33(i)
                READ(30,*) wave_len(i), cs_34(i)
                READ(31,*) wave_len(i), cs_35(i)
                READ(32,*) wave_len(i), cs_41(i)
                READ(33,*) wave_len(i), cs_53(i)
                READ(34,*) wave_len(i), cs_55(i)
                READ(35,*) wave_len(i), cs_56(i)
                READ(36,*) wave_len(i), cs_57(i)
				
				READ(91,*) wave_len(i), cs_61(i)
				READ(92,*) wave_len(i), cs_62(i)
				READ(93,*) wave_len(i), cs_63(i)
				READ(94,*) wave_len(i), cs_64(i)
				READ(95,*) wave_len(i), cs_65(i)
				READ(96,*) wave_len(i), cs_66(i)
				READ(97,*) wave_len(i), cs_67(i)
				READ(98,*) wave_len(i), cs_68(i)
				READ(99,*) wave_len(i), cs_69(i)
				READ(100,*) wave_len(i), cs_70(i)
				READ(101,*) wave_len(i), cs_71(i)
				
				READ(102,*) wave_len(i), cs_72(i)
				READ(103,*) wave_len(i), cs_73(i)
				READ(104,*) wave_len(i), cs_74(i)
				READ(105,*) wave_len(i), cs_75(i)
				READ(106,*) wave_len(i), cs_76(i)
				READ(107,*) wave_len(i), cs_77(i)
				READ(108,*) wave_len(i), cs_78(i)
				READ(109,*) wave_len(i), cs_79(i)
				
				READ(118,*) wave_len(i), cs_80(i)
				READ(119,*) wave_len(i), cs_81(i)
				READ(120,*) wave_len(i), cs_82(i)
				READ(121,*) wave_len(i), cs_83(i)
				READ(122,*) wave_len(i), cs_84(i)
				READ(123,*) wave_len(i), cs_85(i)
				READ(124,*) wave_len(i), cs_86(i)
				READ(125,*) wave_len(i), cs_87(i)
				READ(126,*) wave_len(i), cs_88(i)
				READ(127,*) wave_len(i), cs_89(i)
				READ(128,*) wave_len(i), cs_90(i)
				READ(129,*) wave_len(i), cs_91(i)
				READ(130,*) wave_len(i), cs_92(i)
				
				READ(144,*) wave_len(i), cs_93(i)
				READ(145,*) wave_len(i), cs_94(i)
				READ(146,*) wave_len(i), cs_95(i)
				READ(147,*) wave_len(i), cs_96(i)
				
                ! Wavelength [nm], Quantum yield
                READ(37,*) wave_len(i), qy_1(i)
                READ(38,*) wave_len(i), qy_2(i)
                READ(39,*) wave_len(i), qy_3(i)
                READ(40,*) wave_len(i), qy_4(i)
                READ(41,*) wave_len(i), qy_5(i)
                READ(42,*) wave_len(i), qy_6(i)
                READ(43,*) wave_len(i), qy_7(i)
                READ(44,*) wave_len(i), qy_8(i)
                READ(45,*) wave_len(i), qy_11(i)
                READ(46,*) wave_len(i), qy_12(i)
                READ(47,*) wave_len(i), qy_13(i)
                READ(48,*) wave_len(i), qy_14(i)
                READ(49,*) wave_len(i), qy_15(i)
                READ(50,*) wave_len(i), qy_16(i)
                READ(51,*) wave_len(i), qy_17(i)
                READ(52,*) wave_len(i), qy_18(i)
                READ(53,*) wave_len(i), qy_19(i)
                READ(54,*) wave_len(i), qy_21(i)
                READ(55,*) wave_len(i), qy_22(i)
                READ(56,*) wave_len(i), qy_23(i)
                READ(57,*) wave_len(i), qy_24(i)
                READ(58,*) wave_len(i), qy_31(i)
                READ(59,*) wave_len(i), qy_32(i)
                READ(60,*) wave_len(i), qy_33(i)
                READ(61,*) wave_len(i), qy_34(i)
                READ(62,*) wave_len(i), qy_35(i)
                READ(63,*) wave_len(i), qy_41(i)
                READ(64,*) wave_len(i), qy_51(i)
                READ(65,*) wave_len(i), qy_52(i)
                READ(66,*) wave_len(i), qy_53(i)
                READ(67,*) wave_len(i), qy_54(i)
                READ(68,*) wave_len(i), qy_55(i)
                READ(69,*) wave_len(i), qy_56(i)
                READ(70,*) wave_len(i), qy_57(i)
				
				READ(80,*) wave_len(i), qy_61(i)
				READ(81,*) wave_len(i), qy_62(i)
				READ(82,*) wave_len(i), qy_63(i)
				READ(83,*) wave_len(i), qy_64(i)
				READ(84,*) wave_len(i), qy_65(i)
				READ(85,*) wave_len(i), qy_66(i)
				READ(86,*) wave_len(i), qy_67(i)
				READ(87,*) wave_len(i), qy_68(i)
				READ(88,*) wave_len(i), qy_69(i)
				READ(89,*) wave_len(i), qy_70(i)
				READ(90,*) wave_len(i), qy_71(i)
				
				READ(110,*) wave_len(i), qy_72(i)
				READ(111,*) wave_len(i), qy_73(i)
				READ(112,*) wave_len(i), qy_74(i)
				READ(113,*) wave_len(i), qy_75(i)
				READ(114,*) wave_len(i), qy_76(i)
				READ(115,*) wave_len(i), qy_77(i)
				READ(116,*) wave_len(i), qy_78(i)
				READ(117,*) wave_len(i), qy_79(i)
				
				READ(131,*) wave_len(i), qy_80(i)
				READ(132,*) wave_len(i), qy_81(i)
				READ(133,*) wave_len(i), qy_82(i)
				READ(134,*) wave_len(i), qy_83(i)
				READ(135,*) wave_len(i), qy_84(i)
				READ(136,*) wave_len(i), qy_85(i)
				READ(137,*) wave_len(i), qy_86(i)
				READ(138,*) wave_len(i), qy_87(i)
				READ(139,*) wave_len(i), qy_88(i)
				READ(140,*) wave_len(i), qy_89(i)
				READ(141,*) wave_len(i), qy_90(i)
				READ(142,*) wave_len(i), qy_91(i)
				READ(143,*) wave_len(i), qy_92(i)
				
				READ(148,*) wave_len(i), qy_93(i)
				READ(149,*) wave_len(i), qy_94(i)
				READ(150,*) wave_len(i), qy_95(i)
				READ(151,*) wave_len(i), qy_96(i)
				
            END DO

            ! Extract data from un-interpolated temperature dependent cs and qy
            DO i = 1,60
                READ(71,*) wave_len_1(i), cs_raw_1(i), cs_T_1(i)
                READ(72,*) wave_len_2(i), cs_raw_2(i), cs_T_2(i)
            END DO
            DO i = 1,83
                READ(73,*) wave_len_8(i), cs_raw_8(i), cs_T_8(i)
            END DO
            DO i = 1,150
                READ(74,*) wave_len_11(i), cs_raw_11(i), cs_T_11(i)
                READ(75,*) wave_len_12(i), cs_raw_12(i), cs_T_12(i)
            END DO
            DO i = 1,135
                READ(76,*) wave_len_21(i), cs_raw_21(i), cs_T1_21(i), cs_T2_21(i), cs_T3_21(i)
            END DO
            DO i = 1,21
                READ(77,*) wave_len_51(i), cs_raw_51(i), cs_T_51(i)
            END DO
            DO i = 1,22
                READ(78,*) wave_len_52(i), cs_raw_52(i), cs_T_52(i)
            END DO
            DO i = 1,37
                READ(79,*) wave_len_54(i), cs_raw_54(i), cs_T_54(i)
            END DO

            DO i = 12,79
                CLOSE(i)
            END DO
            

        END IF

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


        ! Calculate temperature dependent absorption cross sections
        cs_rawT_1 = cs_raw_1*exp(cs_T_1/TEMP)
        cs_rawT_2 = cs_raw_2*exp(cs_T_2/TEMP)
        cs_rawT_8 = cs_raw_8*exp(cs_T_8*1.d-3*(TEMP-298))
        cs_rawT_11 = cs_raw_11+cs_T_11*TEMP
        cs_rawT_12 = cs_raw_12+cs_T_12*TEMP
        cs_rawT_21 = cs_raw_21*(1+cs_T1_21*TEMP+cs_T2_21*TEMP**2+cs_T3_21*TEMP**3)
        cs_rawT_51 = cs_raw_51*exp(cs_T_51*(TEMP-298))
        cs_rawT_52 = cs_raw_52*exp(cs_T_52*(TEMP-298))
        cs_rawT_54 = cs_raw_54*exp(cs_T_54*(TEMP-298))

        ! Interpolate to 1 nm wavelength resolution
        ! O3
        DO i = 1,60-1
            cs_int_1(i) = cs_rawT_1(i) + (cs_rawT_1(i+1)-cs_rawT_1(i)) * &
                ((wave_len_1(i)+(wave_len_1(i+1)-wave_len_1(i))/2)-wave_len_1(i)) &
                / (wave_len_1(i+1)-wave_len_1(i))
            cs_int_2(i) = cs_rawT_2(i) + (cs_rawT_2(i+1)-cs_rawT_2(i)) * &
                ((wave_len_2(i)+(wave_len_2(i+1)-wave_len_2(i))/2)-wave_len_2(i)) &
                / (wave_len_2(i+1)-wave_len_2(i))
        END DO
        cs_1(1:40) = 0
        cs_1(41:99) = cs_int_1(1:59)
        cs_1(100:510) = 0
        cs_2(1:40) = 0
        cs_2(41:99) = cs_int_2(1:59)
        cs_2(100:510) = 0
        ! HNO3
        k = 0
        DO i = 1,83-1
            DO j = 0,2-1
                cs_int_8(1+k) = cs_rawT_8(i) + (cs_rawT_8(i+1)-cs_rawT_8(i)) * &
                    ((wave_len_8(i)+(1+j*2)*(wave_len_8(i+1)-wave_len_8(i))/4)-wave_len_8(i)) &
                    / (wave_len_8(i+1)-wave_len_8(i))
            END DO
        END DO
        cs_8(1:100) = cs_int_8(66:165)
        cs_8(101:510) = 0
        ! HCHO
        DO i = 1,150-1
            cs_int_11(i) = cs_rawT_11(i) + (cs_rawT_11(i+1)-cs_rawT_11(i)) * &
                ((wave_len_11(i)+(wave_len_11(i+1)-wave_len_11(i))/2)-wave_len_11(i)) &
                / (wave_len_11(i+1)-wave_len_11(i))
            cs_int_12(i) = cs_rawT_12(i) + (cs_rawT_12(i+1)-cs_rawT_12(i)) * &
                ((wave_len_12(i)+(wave_len_12(i+1)-wave_len_12(i))/2)-wave_len_12(i)) &
                / (wave_len_12(i+1)-wave_len_12(i))
        END DO
        cs_11(1:125) = cs_int_11(25:149)
        cs_11(126:510) = 0
        cs_12(1:125) = cs_int_12(25:149)
        cs_12(126:510) = 0
        ! CH3COCH3
        DO i = 1,135-1
            cs_int_21(i) = cs_rawT_21(i) + (cs_rawT_21(i+1)-cs_rawT_21(i)) * &
                ((wave_len_21(i)+(wave_len_21(i+1)-wave_len_21(i))/2)-wave_len_21(i)) &
                / (wave_len_21(i+1)-wave_len_21(i))
        END DO
        cs_21(1:99) = cs_int_21(36:134)
        cs_21(100:510) = 0
        ! CH3NO3
        k = 0
        DO i = 1,21-1
            DO j = 0,5-1
                cs_int_51(1+k) = cs_rawT_51(i) + (cs_rawT_51(i+1)-cs_rawT_51(i)) * &
                    ((wave_len_51(i)+(1+j*2)*(wave_len_51(i+1)-wave_len_51(i))/10)-wave_len_51(i)) &
                    / (wave_len_51(i+1)-wave_len_51(i))
                k = k + 1
            END DO
        END DO
        cs_51(1:90) = cs_int_51(11:100)
        cs_51(91:510) = 0
        ! C2H5NO3
        k = 0
        DO i = 1,22-1
            DO j = 0,5-1
                cs_int_52(1+k) = cs_rawT_52(i) + (cs_rawT_52(i+1)-cs_rawT_52(i)) * &
                    ((wave_len_52(i)+(1+j*2)*(wave_len_52(i+1)-wave_len_52(i))/10)-wave_len_52(i)) &
                    / (wave_len_52(i+1)-wave_len_52(i))
                k = k + 1
            END DO
        END DO
        cs_52(1:90) = cs_int_52(16:105)
        cs_52(91:510) = 0
        ! i-C3H7NO3
        k = 0
        DO i = 15,37-1
            DO j = 0,5-1
                cs_int_54(1+k) = cs_rawT_54(i) + (cs_rawT_54(i+1)-cs_rawT_54(i)) * &
                    ((wave_len_54(i)+(1+j*2)*(wave_len_54(i+1)-wave_len_54(i))/10)-wave_len_54(i)) &
                    / (wave_len_54(i+1)-wave_len_54(i))
                k = k + 1
            END DO
        END DO
        cs_54(1:110) = cs_int_54(1:110)
        cs_54(111:510) = 0


        ! Calculate the photolysis rates
        DO i = 1, NPHOT
            JVALUES(i) = 0
        END DO
        JVALUES(1)=sum(cs_1*qy_1*i_surface*dlambda*4*pi) ! O3 = O1D
        JVALUES(2)=sum(cs_2*qy_2*i_surface*dlambda*4*pi) ! O3 = O
        JVALUES(3)=sum(cs_3*qy_3*i_surface*dlambda*4*pi) ! H2O2 = OH + OH
        JVALUES(4)=sum(cs_4*qy_4*i_surface*dlambda*4*pi) ! NO2 = NO + O
        JVALUES(5)=sum(cs_5*qy_5*i_surface*dlambda*4*pi) ! NO3 = NO
        JVALUES(6)=sum(cs_6*qy_6*i_surface*dlambda*4*pi) ! NO3 = NO2 + O
        JVALUES(7)=sum(cs_7*qy_7*i_surface*dlambda*4*pi) ! HONO = OH + NO
        JVALUES(8)=sum(cs_8*qy_8*i_surface*dlambda*4*pi) ! HNO3 = OH + NO2

        JVALUES(11)=sum(cs_11*qy_11*i_surface*dlambda*4*pi) ! HCHO = CO + HO2 + HO2
        JVALUES(12)=sum(cs_12*qy_12*i_surface*dlambda*4*pi) ! HCHO = H2 + CO
        JVALUES(13)=sum(cs_13*qy_13*i_surface*dlambda*4*pi) ! CH3CHO = CH3O2 + HO2 + CO
        JVALUES(14)=sum(cs_14*qy_14*i_surface*dlambda*4*pi) ! C2H5CHO = C2H5O2 + HO2 + CO
        JVALUES(15)=sum(cs_15*qy_15*i_surface*dlambda*4*pi) ! C3H7CHO = n-C3H7 + HCO
        JVALUES(16)=sum(cs_16*qy_16*i_surface*dlambda*4*pi) ! C3H7CHO = C2H4 + CH3CHO
        JVALUES(17)=sum(cs_17*qy_17*i_surface*dlambda*4*pi) ! IPRCHO = IC3H7O2 + HO2 + CO
        JVALUES(18)=sum(cs_18*qy_18*i_surface*dlambda*4*pi) ! MACR = MACO3 + HO2
        JVALUES(19)=sum(cs_19*qy_19*i_surface*dlambda*4*pi) ! MACR = CH3CO3 + HCHO + CO + HO2

        JVALUES(21)=sum(cs_21*qy_21*i_surface*dlambda*4*pi) ! CH3COCH3 = CH3CO3 + CH3O2
        JVALUES(22)=sum(cs_22*qy_22*i_surface*dlambda*4*pi) ! MEK = CH3CO3 + C2H5O2
        JVALUES(23)=sum(cs_23*qy_23*i_surface*dlambda*4*pi) ! MVK = C3H6 + CO
        JVALUES(24)=sum(cs_24*qy_24*i_surface*dlambda*4*pi) ! MVK = CH3CO3 + HCHO + CO + HO2

        JVALUES(31)=sum(cs_31*qy_31*i_surface*dlambda*4*pi) !  GLYOX = CO + CO + H2
        JVALUES(32)=sum(cs_32*qy_32*i_surface*dlambda*4*pi) ! GLYOX = HCHO + CO
        JVALUES(33)=sum(cs_33*qy_33*i_surface*dlambda*4*pi) ! GLYOX = CO + CO + HO2 + HO2
        JVALUES(34)=sum(cs_34*qy_34*i_surface*dlambda*4*pi) ! MGLYOX = CH3CO3 + CO + HO2
        JVALUES(35)=sum(cs_35*qy_35*i_surface*dlambda*4*pi) ! MGLYOX = CH3CO3 + CO + HO2

        JVALUES(41)=sum(cs_41*qy_41*i_surface*dlambda*4*pi) ! e.g. CH3OOH->CH3O+OH or C2H5OOH -> C2H5O + OH

        JVALUES(51)=sum(cs_51*qy_51*i_surface*dlambda*4*pi) ! CH3NO3 -> CH3O + NO2
        JVALUES(52)=sum(cs_52*qy_52*i_surface*dlambda*4*pi) ! C2H5NO3 -> C2H5O + NO2
        JVALUES(53)=sum(cs_53*qy_54*i_surface*dlambda*4*pi) ! NC3H7NO3 -> n-C3H7O + NO2
        JVALUES(54)=sum(cs_54*qy_54*i_surface*dlambda*4*pi) ! IC3H7NO3 -> CH3C(O)CH3 + NO2
        JVALUES(55)=sum(cs_55*qy_55*i_surface*dlambda*4*pi) ! TC4H9NO3 -> t-C4H9O + NO2
        JVALUES(56)=sum(cs_56*qy_56*i_surface*dlambda*4*pi) ! NOA -> CH3C(O)CH2(O) + NO2
        JVALUES(57)=sum(cs_57*qy_57*i_surface*dlambda*4*pi) ! NOA -> CH3CO + HCHO + NO2
        
		JVALUES(61)=sum(cs_61*qy_61*i_surface*dlambda*4*pi) ! CL2+hv->2CL Pg1 Photolysis rate Table S14 Bräuer et al. (2013) Ref: Sander et al. (2006)
        JVALUES(62)=sum(cs_62*qy_62*i_surface*dlambda*4*pi) ! CLO+hv->CL+O Pg2 Photolysis rate Table S14 Bräuer et al. (2013) Ref: Sander et al. (2006)
        JVALUES(63)=sum(cs_63*qy_63*i_surface*dlambda*4*pi) ! OCLO+hv->CLO+O Pg3 Photolysis rate Table S14 Bräuer et al. (2013) Ref: Sander et al. (2006)
        JVALUES(64)=sum(cs_64*qy_64*i_surface*dlambda*4*pi) ! CL2O2+hv->CL+CLO2 Pg4 Photolysis rate Table S14 Bräuer et al. (2013) Ref Atkinson et al. (2007)
        JVALUES(65)=sum(cs_65*qy_65*i_surface*dlambda*4*pi) ! CL2O3+hv->CLO+OCLO Pg5 Photolysis rate Table S14 Bräuer et al. (2013) Ref Atkinson et al. (2007)
        JVALUES(66)=sum(cs_66*qy_66*i_surface*dlambda*4*pi) ! HOCL+hv->CL+OH Pg6 Photolysis rate Table S14 Bräuer et al. (2013) Ref Atkinson et al. (2007)
        JVALUES(67)=sum(cs_67*qy_67*i_surface*dlambda*4*pi) ! CLNO+hv->CL+NO Pg7 Photolysis rate Table S14 Bräuer et al. (2013) Ref Atkinson et al. (2007)
        JVALUES(68)=sum(cs_68*qy_68*i_surface*dlambda*4*pi) ! CLNO2+hv->CL+NO2 Pg8 Photolysis rate Table S14 Bräuer et al. (2013) Ref Atkinson et al. (2007)
        JVALUES(69)=sum(cs_69*qy_69*i_surface*dlambda*4*pi) ! CLNO3+hv->CL+NO3 Pg9 Photolysis rate Table S14 Bräuer et al. (2013) Ref Atkinson et al. (2007)
        JVALUES(70)=sum(cs_70*qy_70*i_surface*dlambda*4*pi) ! CLNO3+hv->CLO+NO2 Pg10 Photolysis rate Table S14 Bräuer et al. (2013) Ref: Sander et al. (2006)
        JVALUES(71)=sum(cs_71*qy_71*i_surface*dlambda*4*pi) ! CH3SCH2CL+hv->CH3S+CH2CLO2 gPD07 Photolysis rate Table S7 Hoffmann et al. (2016) Ref: Sander et al. (2006)
        
		JVALUES(72)=sum(cs_72*qy_72*i_surface*dlambda*4*pi) ! Br2 = Br + Br ; // Pg19  Photolysis rate Table S14 Ref: Atkinson et al., 2007
		JVALUES(73)=sum(cs_73*qy_73*i_surface*dlambda*4*pi) ! BrO = Br + O : ; Pg20 Photolysis rate Table S14 Ref: Atkinson et al., 2007
		JVALUES(74)=sum(cs_74*qy_74*i_surface*dlambda*4*pi) ! OBrO = BrO + O : ; Pg21 Photolysis rate Table S14 Ref: Atkinson et al., 2007
		JVALUES(75)=sum(cs_75*qy_75*i_surface*dlambda*4*pi) ! HOBr = Br + OH : ; Pg22 Photolysis rate Table S14 Ref: Atkinson et al., 2007
		JVALUES(76)=sum(cs_76*qy_76*i_surface*dlambda*4*pi) ! BrNO2 = Br + NO2 : ; Pg23 Photolysis rate Table S14 Ref: Atkinson et al., 2007
		JVALUES(77)=sum(cs_77*qy_77*i_surface*dlambda*4*pi) ! BrNO3 = Br + NO3 : ; Pg24 Photolysis rate Table S14 Ref: Atkinson et al., 2007
		JVALUES(78)=sum(cs_78*qy_78*i_surface*dlambda*4*pi) ! BrNO3 = BrO + NO2 : ; Pg25 Photolysis rate Table S14 Ref: Atkinson et al., 2007
		JVALUES(79)=sum(cs_79*qy_79*i_surface*dlambda*4*pi) ! BrCL = Br + CL : ; Pg26   Photolysis rate Table S14 Ref: Atkinson et al., 2007
		
		
		JVALUES(80)=sum(cs_80*qy_80*i_surface*dlambda*4*pi) ! I2 = I + I : J(80); // Pg36 Atkinson et al., 2007
		JVALUES(81)=sum(cs_81*qy_81*i_surface*dlambda*4*pi) ! IO = I + O : J(81); // Pg37 Atkinson et al., 2007
		JVALUES(82)=sum(cs_82*qy_82*i_surface*dlambda*4*pi) ! OIO = I : J(82); // Pg38 
		JVALUES(83)=sum(cs_83*qy_83*i_surface*dlambda*4*pi) ! OIO = IO + O : J(83); // Pg39 
		JVALUES(84)=sum(cs_84*qy_84*i_surface*dlambda*4*pi) ! I2O2 = I + I : J(84); // Pg40
		JVALUES(85)=sum(cs_85*qy_85*i_surface*dlambda*4*pi) ! HI = I + HO2 : J(85); // Pg41 Atkinson et al., 2007
		JVALUES(86)=sum(cs_86*qy_86*i_surface*dlambda*4*pi) ! HOI = I + OH : J(86); // Pg42 Atkinson et al., 2007
		JVALUES(87)=sum(cs_87*qy_87*i_surface*dlambda*4*pi) ! INO = I + NO : J(87); // Pg43 Sander et al., 2006
		JVALUES(88)=sum(cs_88*qy_88*i_surface*dlambda*4*pi) ! INO2 = I + NO2 : J(88); // Pg44 Sander et al., 2006
		JVALUES(89)=sum(cs_89*qy_89*i_surface*dlambda*4*pi) ! INO3 = I + NO3 : J(89); // Pg45 Sander et al., 2006
		JVALUES(90)=sum(cs_90*qy_90*i_surface*dlambda*4*pi) ! INO3 = IO + NO2 : J(90); // Pg46 Sander et al., 2006
		JVALUES(91)=sum(cs_91*qy_91*i_surface*dlambda*4*pi) ! ICL = I + CL : J(91); // Pg47 Atkinson et al., 2007
		JVALUES(92)=sum(cs_92*qy_92*i_surface*dlambda*4*pi) ! IBr = I + Br : J(92); // Pg48 Atkionson et al., 2007		
		
		JVALUES(93)=sum(cs_93*qy_93*i_surface*dlambda*4*pi) ! CHBr3 = Br + CHBr2O2 : J(93); // Pg31 Sander et al., 2015		
		JVALUES(94)=sum(cs_94*qy_94*i_surface*dlambda*4*pi) ! CH2Br2 = Br + CH2BrO2 : J(94); // Pg32 Atkionson et al., 2007		
		JVALUES(95)=sum(cs_95*qy_95*i_surface*dlambda*4*pi) ! COBr2 = Br + Br + CO : J(95); // Pg33 Sander et al., 2015		
		JVALUES(96)=sum(cs_96*qy_96*i_surface*dlambda*4*pi) ! CH3I = I + CH3O2 : J(96); // Pg54 Sander et al., 2015		
		
        !Some safety check
        DO i = 1,NPHOT
            IF (JVALUES(i) < -1e9 .OR. JVALUES(i) > 1e9 ) then
                WRITE(*,*) 'Note: JVALUES(i) has bad value.'
                WRITE(*,*) 'i, JVALUES(i) = ', i, JVALUES(i)
                STOP
            ENDIF
        ENDDO

    END SUBROUTINE getJVALUES

END MODULE reaction_rates_20191031
