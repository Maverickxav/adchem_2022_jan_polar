MODULE thermodynamicsDMS

USE constants
USE acidityDMS
USE second_Precision, ONLY : dp    ! KPP Numerical type
IMPLICIT NONE

PRIVATE

PUBLIC :: thermodyn_AIOMFAC_inorg, thermodyn_AIOMFAC_inorg_bulk, UNIFAC_activitycoefficients_org, UNIFAC_parameters

CONTAINS

SUBROUTINE thermodyn_AIOMFAC_inorg(T,c_p_L,N_bins,y_L,&
    cNH3,cHNO3,cHCl,aw,pCO2,pH,Kprim_HNO3,Kprim_HCl,&
	Kprim_CH3SO3H,Kprim_HIO3,Hprim_NH3,Kprim_NH3,fHSO4,fSO4,fNO3,fCl,&
	fCH3SO3,fHIO3,mHCO3,mCO3,mOH,W)

REAL(dp), INTENT(in)  :: T,cNH3,cHNO3,cHCl,pCO2
REAL(dp), DIMENSION(nr_bins), INTENT(in)  :: N_bins, aw

REAL(dp), DIMENSION(NSPEC_P,nr_bins), INTENT(inout) :: c_p_L
REAL(dp), DIMENSION(NSPEC_P+1,nr_bins), INTENT(inout) :: y_L
REAL(dp), DIMENSION(NSPEC_P+1,nr_bins) :: y_start
REAL(dp), DIMENSION(nr_bins), INTENT(out) :: pH, Kprim_HNO3, Kprim_HCl, Kprim_CH3SO3H, Kprim_HIO3, &
Hprim_NH3, Kprim_NH3,fHSO4, fSO4, fNO3, fCl, fCH3SO3,fHIO3, mHCO3, mCO3, mOH, W  
REAL(dp), PARAMETER :: xw_limit=5D-1  ! Lowest allowed water activity for inorganic salt aerosol particles activity coeff calculation
!REAL(dp), PARAMETER :: aw_limit=2D-1  ! Lowest allowed water activity for inorganic salt aerosol particles activity coeff calculation
REAL(dp), DIMENSION(nr_bins) :: a_w,n_SVI,n_NV,n_CLI,n_CH3SO3I,n_HIO3,n_NIII,n_Na,n_DMA, &
n_H2O_inorg,mH,fNH4,mNIII,mNV,mClI,mCH3SO3I,mNa,mcations,mSVI,mHIO3,pNH3pHNO3_liqid,&
pNH3pHCl_liqid,CNH3CHNO3_liqid,CNH3CHCl_liqid,mCl,mCH3SO3,mNH4,mNO3,cw,cH,Hprim_HCl,Hprim_CH3SO3H,Hprim_HIO3,Hprim_HNO3,yHIO3,yCH3SO3
REAL(dp), DIMENSION(7) :: ni
REAL(dp), DIMENSION(10) :: n_inorg
REAL(dp) :: error_tot, Rprim,KNH4NO3s,KNH4NO3l,&
KNH4NO3_solid,KNH4Cll,KNH4Cl_solid,K_3,K_4,K_5,K_6,Nr_iterations, &
pHold, xw_inorg, CNH3CHNO3_solid, CNH3CHCl_solid, C0, CNH3m_old, CHNO3m_old, &
CHClm_old, CNH3s, KHNO3, KHCl,KCH3SO3H,KHIO3, K_2, H_NH3, K_NH3, H_HCl,H_CH3SO3, &
H_HNO3, H_CH3SO3H, H_HIO3,xw_AIOMFAC
REAL(dp), DIMENSION(8) :: yi,y_old
REAL(dp), DIMENSION(NSPEC_P,nr_bins) :: n_X,c_p_old
REAL(dp) :: lnyHIO3,mtot,A1,A2,A3,A4,A5,A6,A7,T_factor
REAL(dp) :: kappa_org,vp_dry_org,vp_water_org
INTEGER :: j
REAL(dp), DIMENSION(8) :: Rti, Qti
REAL(dp), DIMENSION(3,12) :: bcai
REAL(dp), DIMENSION(2,12) :: ccai


! Combnined KHNO3*H_HNO3 from Table B7 in Jacobson (2005), Temperature dependence assigned to the Henry's law coefficient
KHNO3=12D0! mol/kg HNO3(aq)<->H+ + NO3-
H_HNO3=2.1D5*EXP(29.17*(298.15/T-1D0)+16.83*(1D0+LOG(298.15/T)-298.15/T))  !*EXP(8700D0*(1D0/T-1D0/298D0)) ! mol/(kg atm) Henry's law coeff for HNO3(g)<->HNO3(aq)
K_3=H_HNO3*KHNO3 ! equilibrium constant HNO3(g) <-> H+ + NO3-

!KHCl and KHCl*H_HCl from Table B7 in Jacobson (2005)
KHCl=1.72D6*EXP(23.15*(298.15/T-1D0)) ! 5D11 !*exp(6890.0*(1D0/T-1D0/298D0))/2.9D5 !1.72D6*EXP(23.15*(298./T-1D0)) ! mol/kg HCl(aq)<->H+ + Cl-
K_4=1.97D6*EXP(30.19*(298.15/T-1D0)+19.91*(1D0+LOG(298.15/T)-298.15/T)) !HCl(g) <-> H+ + Cl- mol^2 kg^-2 atm^-1 
H_HCl=K_4/KHCl !HCl(g)<->HCl(aq)

! KCH3SO3H and H_CH3SO3H from COSMOTherm
KCH3SO3H=855.0667*EXP(1980.5*(1D0/T-1D0/298.15D0))! COSMOTherm (Noora)
H_CH3SO3H=1D1*7.7716275D7*EXP(9.9034D3*(1D0/T-1D0/298.15D0)) ! mol/(kg atm) Henry's law coeff for CH3SO2H(g)<->CH3SO2H(aq), COSMOtherm
K_5=H_CH3SO3H*KCH3SO3H; ! equilibrium constant CH3SO3H(g) <-> H+ + CH3SO3- mol^2 kg^-2 atm^-1 COSMOTherm

! KHIO3 and H_HIO3 from COSMOTherm
KHIO3=155.2387*EXP(1729.8*(1D0/T-1D0/298.15D0))! COSMOTherm (Noora)
H_HIO3=5.02572D8*EXP(1.0695D4*(1D0/T-1D0/298.15D0)) ! Henry's law coeff for HIO3(g)<->HIO3(aq), COSMOTherm
K_6=H_HIO3*KHIO3; ! equilibrium constant HIO3(g) <-> H+ + IO3- mol^2 kg^-2 atm^-1 

! Middle range contribution:
! binary cation-anion MR interaction parameters, Table 5 from Zuend et
! al., 2008 and 2011
bcai(:,1)=(/0.182003, 0.243340, 0.8/)         ! H+ Cl-
bcai(:,2)=(/0.210638, 0.122694, 0.8/)        ! H+ NO3-
bcai(:,3)=(/0.286343, -5.99615, 1.36861/)     ! H+ SO42-
bcai(:,4)=(/0.0215532, 0.562966, 0.142442/)   ! H+ HSO4-
bcai(:,5)=(/0.053741, 0.079771, 0.8/)         ! Na+ Cl-
bcai(:,6)=(/0.001164, -0.102546, 0.410453/)   ! Na+ NO3-
bcai(:,7)=(/0.001891, -0.424184, 0.8/)        ! Na+ SO42-
bcai(:,8)=(/0.0153214, 0.4, 0.423635/)        ! Na+ HSO4-
bcai(:,9)=(/0.001520, 0.049074, 0.116801/)    ! NH4+ Cl-
bcai(:,10)=(/-0.000057, -0.171746, 0.260000/) ! NH4+ NO3-
bcai(:,11)=(/0.000373, -0.906075, 0.545109/)  ! NH4+ SO42-
bcai(:,12)=(/0.00759735, 0.143012, 0.203954/) ! NH4+ HSO4-

ccai(:,1)=(/0.033319, 0.504672/)  ! H+ Cl-
ccai(:,2)=(/-0.101736, 1.676420/) ! H+ NO3-
ccai(:,3)=(/-0.535977, 0.9072/)   ! H+ SO42-
ccai(:,4)=(/0.0703842, 0.714194/) ! H+ HSO4-
ccai(:,5)=(/0.024553, 0.562981/)  ! Na+ Cl-
ccai(:,6)=(/0.002535, 0.512657/)  ! Na+ NO3- 
ccai(:,7)=(/-0.223851, 1.053620/) ! Na+ SO42-
ccai(:,8)=(/0.00350072, 0.40/)    ! Na+ HSO4-
ccai(:,9)=(/0.011112, 0.653256/)  ! NH4+ Cl-
ccai(:,10)=(/0.005510, 0.529762/)  ! NH4+ NO3-
ccai(:,11)=(/-0.000379, 0.354206/) ! NH4+ SO42-
ccai(:,12)=(/0.00631184, 0.825386/)! NH4+ HSO4-


! UNIFAC short range (SR) interactions of organics and water:
! Rt and Qt are the relative van der Waals subgroup volume and surface area parameters, respectively.
Rti=(/1.78,&  ! H+  Zuend et al., 2008 table 1, considering dynamic hydration
0.38,&       ! Na+  Zuend et al., 2008 table 1, considering dynamic hydration
0.69,&       ! NH4+   Zuend et al., 2008 table 1, considering dynamic hydration
0.99,&       ! Cl-  Zuend et al., 2008 table 1, considering dynamic hydration 
0.95,&       ! NO3-  Zuend et al., 2008 table 1, considering dynamic hydration
3.34,&       ! SO42-  Zuend et al., 2008 table 1, considering dynamic hydration
1.65,&       ! HSO4-  Zuend et al., 2008 table 1, considering dynamic hydration 
0.92/)       ! H2O, Hansen et al., 1991, EAIM
 

Qti=(/2.7,&    ! H+  Zuend et al., 2008 table 1, considering dynamic hydration
0.62,&       ! Na+  Zuend et al., 2008 table 1, considering dynamic hydration
0.78,&       ! NH4+   Zuend et al., 2008 table 1, considering dynamic hydration
0.99,&       ! Cl-   Zuend et al., 2008 table 1, considering dynamic hydration
0.97,&       ! NO3-   Zuend et al., 2008 table 1, considering dynamic hydration
3.96,&       ! SO42-   Zuend et al., 2008 table 1, considering dynamic hydration
1.4,&        ! HSO4-   Zuend et al., 2008 table 1, considering dynamic hydration
1.4/)        !H2O, Hansen et al., 1991, EAIM

a_w=aw
!where (a_w < aw_limit) a_w = aw_limit


DO j = 1,NSPEC_P
n_X(j,:)=c_p_L(j,:)/(N_bins*1D-6)/Na ! mol of compounds in each singe particle 
END DO

! moles of water soluble inorganic comp. and solvents (water and organic comp.):
n_SVI=n_X(1,:)
n_NV=n_X(2,:)
n_ClI=n_X(3,:)
n_NIII=n_X(4,:)
n_Na=n_X(5,:)
n_CH3SO3I=n_X(9,:)
n_HIO3=n_X(10,:)
n_DMA=n_X(11,:)
n_H2O_inorg=n_X(7,:)
W=n_H2O_inorg*MH2O ! water content in kg of each single inorganic particle fraction
c_p_old=c_p_L
pH=0D0

y_start=y_L

DO j = 1,nr_bins
yi=y_L(1:8,j)

error_tot=1D0
Nr_iterations=0D0

DO WHILE (error_tot>1D-2 .AND. Nr_iterations<=1D3) ! Model can be sensitive to this Error tolerance. 
y_old=yi
pHold=pH(j)

!yCH3SO3(j)=yi(6) ! Assume the same activity coefficient for MSA as H2SO4 yCH3SO3=1D0
yCH3SO3=1D0 ! Assumed
!yHIO3=1D0 ! Assumed

! Molality of ions in potentialy supersaturated solutions in the iorganic + water phase:
mNIII(j)=n_NIII(j)/W(j)     ! molality of NH4+NH3 (mol/kg water)
mNV(j)=n_NV(j)/W(j)         ! molality of NO3- (mol/kg water)
mClI(j)=n_ClI(j)/W(j)       ! molality of Cl- (mol /kg water)
mcations(j)=n_Na(j)/W(j)+n_DMA(j)/W(j)         ! molality of Na+ and DMA+ (mol/kg water)
mSVI(j)=n_SVI(j)/W(j)       ! molality of S(VI)
mCH3SO3I(j)=n_CH3SO3I(j)/W(j) ! molality of CH3SO3- (mol /kg water)
mHIO3(j)=n_HIO3(j)/W(j) ! molality of IO3- (mol /kg water)

!mCOOHtot(j)=0D0
mtot=mHIO3(j)!mNIII(j)+mNV(j)+mClI(j)+mNa(j)+mSVI(j)+mCH3SO3I(j)+mHIO3(j)

! Empirical IO3- activity coefficients Goldman et al Journal of Solution Chemistry, Vol. 3, No. 8, 1974
A1=-0.687624D0; A2=0.200921D0; A3=-1.36944D-2; A4=1.49411D-3
A5=-1.02754D-4;A6=3.53394D-6; A7=-4.63371D-8
lnyHIO3 = 3D0*A1*sqrt(mtot) + 2D0*A2*mtot + (3D0/2D0)*A3*(mtot**2) + (4D0/3D0)*A4*(mtot**3) +&
    (5D0/4D0)*A5*(mtot**4)+(6D0/5D0)*A6*(mtot**5)+(7D0/6D0)*A7*(mtot**6)
yHIO3(j)=EXP(lnyHIO3);
IF (yHIO3(j)<1D-2) THEN
yHIO3(j)=1D-2
END IF

CALL ACIDITY_DMS(mNV(j),mNIII(j),mSVI(j),mClI(j),mcations(j),mCH3SO3I(j),T,&
yi(1),yi(7),yi(6),yi(3),yi(5),yi(4),yCH3SO3(j),pCO2,&
mH(j),fHSO4(j),fSO4(j),fNO3(j),fCl(j),fCH3SO3(j),fNH4(j),mHCO3(j),mCO3(j),mOH(j))


fHIO3(j)=KHIO3/(mH(j)*yi(1)*yHIO3(j)+KHIO3) ! [IO3-]/([HIO3]+[IO3-]), Assuming yIO3=1

pH(j)=-LOG10(mH(j))

ni=(/mH(j)*W(j), n_Na(j), n_NIII(j)*fNH4(j), n_ClI(j)*fCl(j), &
n_NV(j)*fNO3(j), n_SVI(j)*fSO4(j), n_SVI(j)*fHSO4(j)/) ! Update ion composition (moles) for ADCHAM

n_inorg=(/mH(j)*W(j), n_Na(j), n_NIII(j), n_ClI(j), &
n_NV(j), n_SVI(j)*fSO4(j), n_SVI(j)*fHSO4(j), n_CH3SO3I(j), 2D0*n_HIO3(j), n_DMA(j)/) ! Update ion composition (moles), assume HIO3 => 1 H+ and 1 IO3-

xw_inorg=a_w(j)/yi(8) ! updated mole fraction of water in the inorg particle fraction in size bin j
! For safety, keep xw_inorg<1D0 (Added by Carlton, after debugging)
IF (xw_inorg >= 1D0) THEN
         xw_inorg = 1D0-1D-10
END IF

xw_AIOMFAC=xw_inorg
! For AIOMFAC activity coeff calculations the water mole fraction has to be >0.3 
if (xw_AIOMFAC<xw_limit) xw_AIOMFAC=xw_limit

n_H2O_inorg(j)=xw_AIOMFAC*SUM(n_inorg)/(1D0-xw_AIOMFAC) ! updated moles of water in inorganic  particle phase at equilibrium

CALL AIOMFAC_activitycoefficients_H2O_inorg(MH2O,ni,&
n_H2O_inorg(j),Rti,Qti,bcai,ccai,T,yi) ! activity coefficients inorganic particle phase

error_tot=sum(abs(yi-y_old))+abs(pHold-pH(j))
yi=y_old+(yi-y_old)/(2D0+Nr_iterations/1D2)

! Update the particle water content:
xw_inorg=a_w(j)/yi(8) ! updated mole fraction of water in the inorg particle fraction in size bin j
! For safety, keep xw_inorg<1D0 (Added by Carlton, after debugging)
IF (xw_inorg >= 1D0) THEN
         xw_inorg = 1D0-1D-10
END IF
n_H2O_inorg(j)=xw_inorg*SUM(n_inorg)/(1D0-xw_inorg) ! updated moles of water in inorganic  particle phase at equilibrium
W(j)=n_H2O_inorg(j)*MH2O ! water content in kg of each single inorganic particle fraction

Nr_iterations=Nr_iterations+1D0

IF (Nr_iterations>999.0) THEN
write(*,*) 'Problem with error convergence in thermodynamics code' 
write(*,*) error_tot
END IF
END DO
y_L(1:8,j)=yi

n_X(7,j)=n_H2O_inorg(j)
c_p_L(:,j)=n_X(:,j)*N_bins(j)*1D-6*Na

! Water content in organic phase, Kappa Köhler theory
!kappa_org=0.1
!vp_dry_org=SUM(c_p_L(12:NSPEC_P,j)*VX(12:NSPEC_P))
!vp_water_org=kappa_org*vp_dry_org/(1D0/aw(j)-1D0)
!c_p_L(8,j)=vp_water_org/VX(8)
END DO

!y_L=(y_start*y_L)**(1D0/2D0) ! Prevent too large steep changes in y_L within one time step

! molality of NH4+, NO3-, Cl- and CH3SO3-
mNH4=mNIII*fNH4
mNO3=mNV*fNO3
mCl=mCLI*fCl

KNH4NO3s=14.9*EXP(-10.4*(298./T-1D0)+17.56*(1D0+LOG(298./T)-298./T)) ! mol^2/kg^2  NH4NO3(s)<-> NH4+ + NO3-
KNH4NO3l=2.58D17*EXP(64.02*(298./T-1D0)+11.44*(1D0+LOG(298./T)-298./T)) ! (mol^2 kg^-2 atm^-2) [NH3(g)]+[HNO3(g)]<->NH4+ + NO3-
KNH4NO3_solid=KNH4NO3s/KNH4NO3l ! NH4NO3(s) <-> NH3(g) + HNO3(g)

KNH4Cll=2.12D17*EXP(65.08*(298./T-1D0)+14.51*(1D0+LOG(298./T)-298./T)) ! (mol^2 kg^-2 atm^-2) [NH3(g)]+[HCl(g)]<->NH4+ + Cl-
KNH4Cl_solid=1.039D-16*EXP(-71.04*(298./T-1D0)+2.4*(1D0+LOG(298./T)-298./T)) ! NH4Cl(s) <-> NH3(g) + HCl(g)


pNH3pHNO3_liqid=y_L(5,:)*y_L(3,:)*mNO3*mNH4/KNH4NO3l ! pHNO3*pNH3 over particle surface in equilibrium (atm^2)

Rprim=82.06D-6 ! m^3 atm mol^-1 K^-1

CNH3CHNO3_liqid=pNH3pHNO3_liqid*(Rprim*T)**(-2D0)
CNH3CHNO3_solid=KNH4NO3_solid*(Rprim*T)**(-2D0)


pNH3pHCl_liqid=y_L(4,:)*y_L(3,:)*mCl*mNH4/KNH4Cll ! pHCl*pNH3 over particle surface in equilibrium (atm^2)

CNH3CHCl_liqid=pNH3pHCl_liqid*(Rprim*T)**(-2D0)
CNH3CHCl_solid=KNH4Cl_solid*(Rprim*T)**(-2D0)

CNH3m_old=cNH3/Na*1D6        ! vapor mole concentration (mol m^-3)
CHNO3m_old=cHNO3/Na*1D6      ! vapor mole concentration (mol m^-3)
CHClm_old=cHCl/Na*1D6        ! vapor mole concentration (mol m^-3)

! If NH4NO3 but not NH4Cl is present:
!CHNO3s=0D0
!C0=CNH3m_old-CHNO3m_old
!CNH3s=0.5*C0+0.5*SQRT(C0**2D0+4D0*(KNH4NO3_solid)*(Rprim*T)**(-2D0))
!where (CNH3CHNO3_liqid>CNH3CHNO3_solid) CHNO3s=CNH3CHNO3_solid/CNH3s ! HNO3 saturation concentration above a NH4NO3 solid salt surface

! If NH4Cl but not NH4NO3 is present:
!CHCls=0D0
!C0=CNH3m_old-CHClm_old
!CNH3s=0.5*C0+0.5*SQRT(C0**2D0+4D0*(KNH4Cl_solid)*(Rprim*T)**(-2D0))
!where (CNH3CHCl_liqid>CNH3CHCl_solid) CHCls=CNH3CHCl_solid/CNH3s ! HCl saturation concentration above a NH4Cl solid salt surface


! If both NH4NO3 and NH4Cl is present:
!C0=CNH3m_old-CHNO3m_old-CHClm_old
!CNH3s=0.5*C0+0.5*sqrt(C0**2D0+4D0*(KNH4NO3_solid+KNH4Cl_solid)*(Rprim*T)**(-2D0))
!where (CNH3CHNO3_liqid>CNH3CHNO3_solid .AND. CNH3CHCl_liqid>CNH3CHCl_solid) & 
!CHNO3s=CNH3CHNO3_solid/CNH3s ! HNO3 saturation concentration above a NH4NO3 solid salt surface
!where (CNH3CHNO3_liqid>CNH3CHNO3_solid .AND. CNH3CHCl_liqid>CNH3CHCl_solid) & 
!CHCls=CNH3CHCl_solid/CNH3s ! HCl saturation concentration above a NH4Cl solid salt surface

!CHCls=0D0
! Adjusted equilibrium coefficient for HNO3, HCl and CH3SO3
! dissolution ( Jacobson, 2005):
!KCH3SO3H=10**(-1.9) ! mol/kg CH3SO3H(aq)<->H+ + CH3SO3-



cw=n_H2O_inorg*N_bins ! mol water /m^3 air in each size bin

cH=mH*W*N_bins ! mol H+ /m^3 air;
Hprim_HNO3=MH2O*cw*Rprim*T*K_3/(mH*y_L(1,:)*y_L(5,:))  
Hprim_HCl=MH2O*cw*Rprim*T*K_4/(mH*y_L(1,:)*y_L(4,:))
Hprim_CH3SO3H=MH2O*cw*Rprim*T*K_5/(mH*y_L(1,:)*yCH3SO3)  ! Assume yCH3SO3=1
Hprim_HIO3=MH2O*cw*Rprim*T*K_6/(mH*y_L(1,:)*yHIO3)  ! Assume yHIO3=1


Kprim_HNO3=Hprim_HNO3+H_HNO3*KHNO3*(MH2O*cw)**2D0*Rprim*T/(cH*y_L(1,:)*y_L(5,:))
Kprim_HCl=Hprim_HCl+H_HCl*KHCl*(MH2O*cw)**2D0*Rprim*T/(cH*y_L(1,:)*y_L(4,:)) 
Kprim_CH3SO3H=Hprim_CH3SO3H+H_CH3SO3H*KCH3SO3H*(MH2O*cw)**2D0*Rprim*T/(cH*y_L(1,:)*yCH3SO3) ! Assume yCH3SO3=1
Kprim_HIO3=Hprim_HIO3+H_HIO3*KHIO3*(MH2O*cw)**2D0*Rprim*T/(cH*y_L(1,:)*yHIO3) ! Assume yHIO3=1


K_2=1.03D11*EXP(34.81*(298./T-1D0)-5.39*(1D0+LOG(298./T)-298./T)) ! equilibrium constant NH3(g)+H+ <->NH4+
H_NH3=57.6*EXP(13.79*(298./T-1D0)-5.39*(1D0+LOG(298./T)-298./T)) ! mol/(kg atm) NH3(g)<->NH3(aq)
K_NH3=K_2/H_NH3 ! kg/mol NH3(aq)+H+<->NH4+
Hprim_NH3=H_NH3*Rprim*T*MH2O*cw ! mol/mol
Kprim_NH3=K_NH3/(MH2O*cw)*y_L(1,:)/y_L(3,:) ! m^3/mol
END SUBROUTINE thermodyn_AIOMFAC_inorg

SUBROUTINE thermodyn_AIOMFAC_inorg_bulk(T,c_p_bulk,pCO2,pH_bulk,yi_bulk)
		
REAL(dp), INTENT(in)  :: T,pCO2
REAL(dp), DIMENSION(NSPEC_P) :: VX
REAL(dp), DIMENSION(NSPEC_P), INTENT(in) :: c_p_bulk
REAL(dp), INTENT(inout) :: pH_bulk
REAL(dp), INTENT(inout), DIMENSION(8) :: yi_bulk
REAL(dp), DIMENSION(8) :: y_old
REAL(dp) :: pHold,W  
REAL(dp) :: n_SVI,n_NV,n_CLI,n_CH3SO3I,n_HIO3,n_NIII,n_Na, n_DMA,&
n_H2O_inorg,mH,fNH4,fHSO4, fSO4, fNO3, fCl, fCH3SO3,fHIO3, mHCO3, mCO3, mOH,&
mNIII,mNV,mClI,mCH3SO3I,mNa,mcations,mSVI,mHIO3,mCl,mCH3SO3,mNH4,mNO3,yHIO3,yCH3SO3
REAL(dp), DIMENSION(7) :: ni
REAL(dp), DIMENSION(10) :: n_inorg
REAL(dp), DIMENSION(NSPEC_P) :: n_X
REAL(dp) :: lnyHIO3,mtot,A1,A2,A3,A4,A5,A6,A7,T_factor
INTEGER :: j

REAL(dp), DIMENSION(8) :: Rti, Qti
REAL(dp), DIMENSION(3,12) :: bcai
REAL(dp), DIMENSION(2,12) :: ccai
REAL(dp) :: error_tot,Nr_iterations
! Middle range contribution:
! binary cation-anion MR interaction parameters, Table 5 from Zuend et
! al., 2008 and 2011
bcai(:,1)=(/0.182003, 0.243340, 0.8/)         ! H+ Cl-
bcai(:,2)=(/0.210638, 0.122694, 0.8/)        ! H+ NO3-
bcai(:,3)=(/0.286343, -5.99615, 1.36861/)     ! H+ SO42-
bcai(:,4)=(/0.0215532, 0.562966, 0.142442/)   ! H+ HSO4-
bcai(:,5)=(/0.053741, 0.079771, 0.8/)         ! Na+ Cl-
bcai(:,6)=(/0.001164, -0.102546, 0.410453/)   ! Na+ NO3-
bcai(:,7)=(/0.001891, -0.424184, 0.8/)        ! Na+ SO42-
bcai(:,8)=(/0.0153214, 0.4, 0.423635/)        ! Na+ HSO4-
bcai(:,9)=(/0.001520, 0.049074, 0.116801/)    ! NH4+ Cl-
bcai(:,10)=(/-0.000057, -0.171746, 0.260000/) ! NH4+ NO3-
bcai(:,11)=(/0.000373, -0.906075, 0.545109/)  ! NH4+ SO42-
bcai(:,12)=(/0.00759735, 0.143012, 0.203954/) ! NH4+ HSO4-

ccai(:,1)=(/0.033319, 0.504672/)  ! H+ Cl-
ccai(:,2)=(/-0.101736, 1.676420/) ! H+ NO3-
ccai(:,3)=(/-0.535977, 0.9072/)   ! H+ SO42-
ccai(:,4)=(/0.0703842, 0.714194/) ! H+ HSO4-
ccai(:,5)=(/0.024553, 0.562981/)  ! Na+ Cl-
ccai(:,6)=(/0.002535, 0.512657/)  ! Na+ NO3- 
ccai(:,7)=(/-0.223851, 1.053620/) ! Na+ SO42-
ccai(:,8)=(/0.00350072, 0.40/)    ! Na+ HSO4-
ccai(:,9)=(/0.011112, 0.653256/)  ! NH4+ Cl-
ccai(:,10)=(/0.005510, 0.529762/)  ! NH4+ NO3-
ccai(:,11)=(/-0.000379, 0.354206/) ! NH4+ SO42-
ccai(:,12)=(/0.00631184, 0.825386/)! NH4+ HSO4-


! UNIFAC short range (SR) interactions of organics and water:
! Rt and Qt are the relative van der Waals subgroup volume and surface area parameters, respectively.
Rti=(/1.78,&  ! H+  Zuend et al., 2008 table 1, considering dynamic hydration
0.38,&       ! Na+  Zuend et al., 2008 table 1, considering dynamic hydration
0.69,&       ! NH4+   Zuend et al., 2008 table 1, considering dynamic hydration
0.99,&       ! Cl-  Zuend et al., 2008 table 1, considering dynamic hydration 
0.95,&       ! NO3-  Zuend et al., 2008 table 1, considering dynamic hydration
3.34,&       ! SO42-  Zuend et al., 2008 table 1, considering dynamic hydration
1.65,&       ! HSO4-  Zuend et al., 2008 table 1, considering dynamic hydration 
0.92/)       ! H2O, Hansen et al., 1991, EAIM
 

Qti=(/2.7,&    ! H+  Zuend et al., 2008 table 1, considering dynamic hydration
0.62,&       ! Na+  Zuend et al., 2008 table 1, considering dynamic hydration
0.78,&       ! NH4+   Zuend et al., 2008 table 1, considering dynamic hydration
0.99,&       ! Cl-   Zuend et al., 2008 table 1, considering dynamic hydration
0.97,&       ! NO3-   Zuend et al., 2008 table 1, considering dynamic hydration
3.96,&       ! SO42-   Zuend et al., 2008 table 1, considering dynamic hydration
1.4,&        ! HSO4-   Zuend et al., 2008 table 1, considering dynamic hydration
1.4/)        !H2O, Hansen et al., 1991, EAIM


n_X=c_p_bulk/Na ! mol of compounds cm^-3 

! moles of water soluble inorganic comp. and solvents (water and organic comp.):
n_SVI=n_X(1)
n_NV=n_X(2)
n_ClI=n_X(3)
n_NIII=n_X(4)
n_Na=n_X(5)
n_CH3SO3I=n_X(9)
n_HIO3=n_X(10)
n_DMA=n_X(11)
n_H2O_inorg=n_X(7)
W=n_H2O_inorg*MH2O ! water content in kg/cm^3

error_tot=1D0
Nr_iterations=0D0
mH=10**(-pH_bulk)

DO WHILE (error_tot>1D-2 .AND. Nr_iterations<=1D3) ! Model can be sensitive to this Error tolerance. 
y_old=yi_bulk
pHold=pH_bulk

yCH3SO3=1D0 ! Assumed
!yHIO3=1D0 ! Assumed


! Molality of ions in potentialy supersaturated solutions in the iorganic + water phase:
mNIII=n_NIII/W     ! molality of NH4+NH3 (mol/kg water)
mNV=n_NV/W         ! molality of NO3- (mol/kg water)
mClI=n_ClI/W       ! molality of Cl- (mol /kg water)
mcations=n_Na/W+n_DMA/W         ! molality of Na+ and DMA+ (mol/kg water)
mSVI=n_SVI/W       ! molality of S(VI)
mCH3SO3I=n_CH3SO3I/W ! molality of CH3SO3- (mol /kg water)
mHIO3=n_HIO3/W ! molality of IO3- (mol /kg water)

mtot=mHIO3!mNIII(j)+mNV(j)+mClI(j)+mNa(j)+mSVI(j)+mCH3SO3I(j)+mHIO3(j)
! Empirical IO3- activity coefficients Goldman et al Journal of Solution Chemistry, Vol. 3, No. 8, 1974
A1=-0.687624D0; A2=0.200921D0; A3=-1.36944D-2; A4=1.49411D-3
A5=-1.02754D-4;A6=3.53394D-6; A7=-4.63371D-8
lnyHIO3 = 3D0*A1*sqrt(mtot) + 2D0*A2*mtot + (3D0/2D0)*A3*(mtot**2) + (4D0/3D0)*A4*(mtot**3) +&
    (5D0/4D0)*A5*(mtot**4)+(6D0/5D0)*A6*(mtot**5)+(7D0/6D0)*A7*(mtot**6)
yHIO3=EXP(lnyHIO3);
IF (yHIO3<1D-2) THEN
yHIO3=1D-2
END IF

CALL ACIDITY_DMS(mNV,mNIII,mSVI,mClI,mcations,mCH3SO3I,T,&
yi_bulk(1),yi_bulk(7),yi_bulk(6),yi_bulk(3),yi_bulk(5),yi_bulk(4),yCH3SO3,pCO2,&
mH,fHSO4,fSO4,fNO3,fCl,fCH3SO3,fNH4,mHCO3,mCO3,mOH)


ni=(/mH*W, n_Na, n_NIII*fNH4, n_ClI*fCl, &
n_NV*fNO3, n_SVI*fSO4, n_SVI*fHSO4/) ! Update ion composition (moles) for ADCHAM

CALL AIOMFAC_activitycoefficients_H2O_inorg(MH2O,ni,&
n_H2O_inorg,Rti,Qti,bcai,ccai,T,yi_bulk) ! activity coefficients inorganic particle phase

pH_bulk=-LOG10(mH)
error_tot=sum(abs(yi_bulk-y_old))+abs(pHold-pH_bulk)
Nr_iterations=Nr_iterations+1D0

IF (Nr_iterations>999.0) THEN
write(*,*) 'Problem with error convergence in thermodynamics code' 
write(*,*) error_tot
END IF
END DO

END SUBROUTINE thermodyn_AIOMFAC_inorg_bulk

SUBROUTINE AIOMFAC_activitycoefficients_H2O_inorg(MH2O,ni,n_H2O,Rti,Qti,bcai,ccai,T,yi) 
INTEGER :: Nr_ions=7,j,i
REAL(dp), DIMENSION(8), INTENT(inout) :: yi
REAL(dp), DIMENSION(7), INTENT(in)  :: ni
REAL(dp), INTENT(in) :: MH2O, n_H2O, T
REAL(dp), DIMENSION(8) :: yi_old
REAL(dp) :: qw=0.998D3 ! Density of pure water
REAL(dp) :: Tref=298.15
REAL(dp) :: ew,A_DH,b_DH, Rcc, Qcca, Rest, I_strength
REAL(dp), DIMENSION(7) :: zi
REAL(dp), DIMENSION(3) :: mc, Alpha2, Alpha3, Beta2
REAL(dp), DIMENSION(4) :: ma
REAL(dp), DIMENSION(3,3) :: Beta3
REAL(dp), DIMENSION(4,3) :: Gamma

REAL(dp), DIMENSION(8), INTENT(in)    :: Rti, Qti
REAL(dp), DIMENSION(3,12), INTENT(in) :: bcai
REAL(dp), DIMENSION(2,12), INTENT(in) :: ccai

REAL(dp), DIMENSION(12) :: Bca, Cca, Bcaprim, Ccaprim
REAL(dp), DIMENSION(8) :: x,lnyLR,lnyMR,lnySR,lnyC,lnyR,Qm,O,Q,l,lny
REAL(dp), DIMENSION(7) :: mi,lnySRref

! ni(1)=H+  ni(2)=Na+ ni(3)=NH4+ ni(4)=Cl- ni(5)=NO3- ni(6)=SO42- ni(7)=HSO4-
! ni number of moles of the ionic components
! nw number of moles of the solvent components

yi_old=yi
zi=1D0;  zi(6)=2D0 ! number of elemental charges of ion i

x(1:Nr_ions)=ni; x(Nr_ions+1)=n_H2O; 
x=x/SUM(x) ! Mole fractions  
DO j=1,Nr_ions
IF (x(j)<1D-10) THEN
x(j)=1D-10 ! Prevent divided by 0
END IF
END DO

ew=78.051+31989.4*(1D0/T-1D0/Tref) ! relative static permittivity (dimensionless) for water

mi=ni/(n_H2O*MH2O) ! Molality (mol/kg) of ions in the aqueous solution

I_strength=0.5D0*SUM(mi*zi**2) ! Ion strength (mol/kg)

A_DH=1.327757E5*SQRT(qw)/((ew*T)**(3D0/2D0)) !  Debye-Hückle parameter (kg^1/2 mol^-1/2)
b_DH=6.359696*SQRT(qw)/SQRT(ew*T)


! Long range contribution:

! long range (LR) activity coefficient expression for the ions (inorganic compounds):
lnyLR(1:Nr_ions)=(-zi**2D0*A_DH*SQRT(I_strength))/(1D0+b_DH*SQRT(I_strength))

! long range (LR) activity coefficient expression for the solvens (organic compounds and water):
lnyLR(Nr_ions+1)=2D0*A_DH*MH2O/(b_DH**3D0)*(1D0+b_DH*SQRT(I_strength)-&
1D0/(1D0+b_DH*sqrt(I_strength))-2D0*log(1D0+b_DH*sqrt(I_strength))) 
 
 ! Mid range contribution:
Bca=bcai(1,:)+bcai(2,:)*EXP(-bcai(3,:)*SQRT(I_strength))
Cca=ccai(1,:)*EXP(-ccai(2,:)*SQRT(I_strength))
Bcaprim=-0.5D0*bcai(2,:)*EXP(-bcai(3,:)*SQRT(I_strength))*bcai(3,:)*I_strength**(-1D0/2D0) ! partial derivative dBca/dI
Ccaprim=-0.5D0*ccai(1,:)*EXP(-ccai(2,:)*SQRT(I_strength))*ccai(2,:)*I_strength**(-1D0/2D0) ! partial derivative  dCca/dI

mc=mi(1:3) ! molality of cations +
ma=mi(4:7) ! molality of anions -

Rcc=-1.54486D-1 ! Table 6 from Zuend et al., 2011 (kg mol^-1) aqueous electrolyte interaction parameter for NH4+ H+ interaction
Qcca=4.48354D-4 ! Table 6 from Zuend et al., 2011 (kg mol^-1) aqueous electrolyte interaction parameter for NH4+ H+ HSO4- interaction


DO j = 1,3
Alpha2(j)=SUM((Bca((j-1)*4+1:j*4)+I_strength*Bcaprim((j-1)*4+1:j*4))*mc(j)*ma) 
Alpha3(j)=SUM((2D0*Cca((j-1)*4+1:j*4)+I_strength*Ccaprim((j-1)*4+1:j*4))*mc(j)*ma)
Beta2(j)=SUM(Bcaprim((j-1)*4+1:j*4)*mc(j)*ma);

DO i = 1,3
Beta3(i,j)=SUM((Cca((j-1)*4+1:j*4)*ABS(zi(i))+&
    Ccaprim((j-1)*4+1:j*4)*(zi(i)**2D0)/2D0*SUM(mi*ABS(zi)))*mc(j)*ma)    
END DO

DO i = 1,4
Gamma(i,j)=SUM((Cca((j-1)*4+1:j*4)*abs(zi(3+i))+&
Ccaprim((j-1)*4+1:j*4)*(zi(3+i)**2D0)/2D0*SUM(mi*ABS(zi)))*mc(j)*ma)
END DO
END DO


DO j = 1,3
IF (j == 1) THEN
Rest=Rcc*mc(3)+Qcca*mc(3)*ma(4)
ELSE IF (j == 2) THEN
Rest = 0D0
ELSE
Rest=Rcc*mc(1)+Qcca*mc(1)*ma(4) 
END IF
! eq. 19 Zuend et al., 2008
lnyMR(j)=SUM(Bca((j-1)*4+1:j*4)*ma)+((zi(j)**2D0)/2D0)*SUM(Beta2)+&   
SUM(Cca((j-1)*4+1:j*4)*ma)*SUM(mi*abs(zi))+SUM(Beta3(j,:))+Rest
END DO

DO j = 1,4
IF (j==4) THEN
Rest=Qcca*mc(1)*mc(3)
ELSE
Rest=0D0
END IF
! eq. 20 Zuend et al., 2008

lnyMR(3+j)=SUM(Bca((/j,4+j,8+j/))*mc)+(zi(3+j)**2D0)/2D0*SUM(Beta2)+&
    SUM(Cca((/j,4+j,8+j/))*mc)*SUM(mi*abs(zi))+SUM(Gamma(j,:))+Rest
END DO

! eq. 17 Zuend et al., 2008
lnyMR(Nr_ions+1)=-MH2O*SUM(Alpha2)-&  
    MH2O*SUM(mi*ABS(zi))*SUM(Alpha3)-MH2O*Rcc*mc(3)*mc(1)-& 
    MH2O*2D0*Qcca*mc(3)*mc(1)*ma(4) 

l=1D1/2D0*(Rti-Qti)-(Rti-1D0) ! Eq. 25 Zuend et al., 2008
O=Rti*x/SUM(Rti*x)
Q=Qti*x/SUM(qti*x)

lnyC=LOG(O/x)+1D1/2D0*Qti*LOG(Q/O)+l-O/x*SUM(x*l) ! combinatory part of UNIFAC (Fredenslund et al., 1975)

Qm=Qti*x/SUM(Qti*x) ! Relative surface area fraction of subgroup m


lnyR=Qti*(1D0-LOG(SUM(Qm))-SUM(Qm/SUM(Qm))) ! Eq. 29, Zuend et al., 2008

lnySR=lnyC+lnyR
lnySRref=LOG(Rti(1:Nr_ions)/Rti(Nr_ions+1))+1D0-Rti(1:Nr_ions)/Rti(Nr_ions+1)+&
1D1/2D0*Qti(1:Nr_ions)*(LOG(Rti(Nr_ions+1)*Qti(1:Nr_ions)/(Rti(1:Nr_ions)*&
Qti(Nr_ions+1)))-1D0+Rti(1:Nr_ions)*Qti(Nr_ions+1)/(Rti(Nr_ions+1)*Qti(1:Nr_ions))) ! equation 32 from Zuend et al., 2008

lnySR(1:Nr_ions)=lnySR(1:Nr_ions)-lnySRref

lny(1:Nr_ions)=lnyLR(1:Nr_ions)+lnyMR(1:Nr_ions)+lnySR(1:Nr_ions)-LOG(1D0+MH2O*SUM(mi))  ! -log(1+MH2O*SUM(mi)) gives the activity coefficients on a molality basis!!!!
lny(Nr_ions+1)=lnyLR(Nr_ions+1)+lnyMR(Nr_ions+1)+lnySR(Nr_ions+1)

yi=EXP(lny) ! single compound activity coefficient

!yi=1D0
! Put lower and upper bounds on the activity coefficients. 
!where (yi<1D-2) yi=1D-2
!where (yi>1D2) yi=1D2

END SUBROUTINE AIOMFAC_activitycoefficients_H2O_inorg


SUBROUTINE UNIFAC_activitycoefficients_org(ns,Qt,Rt,v,lnTj,Ymn,y_org)

REAL(dp), DIMENSION(NCOND+1), INTENT(in) :: ns
REAL(dp), DIMENSION(52), INTENT(in) :: Qt,Rt 
REAL(dp), DIMENSION(52,NCOND+1), INTENT(in) :: v, lnTj
REAL(dp), DIMENSION(52,52), INTENT(in) :: Ymn
REAL(dp), DIMENSION(NCOND+1) :: r,qq,l,x,O,Q,lnyC,lnyR
REAL(dp) :: qw, sumXm
REAL(dp), DIMENSION(52) :: Qmj,Xm,Qm,lnT, sumQmj,sumQm,part2,part3
REAL(dp), DIMENSION(NCOND+1), INTENT(out) :: y_org
INTEGER :: i,j

qw=1D3 ! density of water (kg/m^3)

sumXm=0D0 
DO j=1,52
sumXm=sumXm+SUM(v(j,:)*ns)
END DO

DO j=1,52
Xm(j)=SUM(v(j,:)*ns)/sumXm ! mole fraction subgroup m in the mixture, intotal 52 different subgroups
END DO

x=ns/SUM(ns) ! mole fraction of compound j

! UNIFAC
r=0D0 
qq=0D0 

DO j=1,NCOND+1
r(j)=SUM(v(:,j)*Rt) 
qq(j)=SUM(v(:,j)*Qt) 
END DO

l=1D1/2D0*(r-qq)-(r-1D0)
O=r*x/SUM(r*x)
Q=qq*x/SUM(qq*x)
lnyC=LOG(O/x)+1D1/2D0*qq*LOG(Q/O)+l-O/x*SUM(x*l) ! combinatory part of UNIFAC (Fredenslund et al., 1975)

Qm=Qt*Xm/SUM(Qt*Xm) ! Relative surface area fraction of subgroup m

lnT=0D0 
sumQm=0D0  
DO i=1,52
   DO j=1,52
    sumQm(j)=Qm(j)*Ymn(i,j)/SUM(Qm*Ymn(:,j))
   END DO
lnT(i)=Qt(i)*(1D0-LOG(SUM(Qm*Ymn(:,i)))-SUM(sumQm)) ! Eq. 29, Zuend et al., 2008
END DO


lnyR=0D0 
DO j=1,NCOND+1
lnyR(j)=SUM(v(:,j)*(lnT-lnTj(:,j)))
END DO

y_org=EXP(lnyC+lnyR) ! single compound activity coefficient
END SUBROUTINE UNIFAC_activitycoefficients_org

SUBROUTINE UNIFAC_parameters(T,UNIFAC_groups,Qt,Rt,v,lnTj,Ymn)
REAL(dp), INTENT(in) :: T
REAL(dp), DIMENSION(52,NCOND), INTENT(in) :: UNIFAC_groups
REAL(dp), DIMENSION(52), INTENT(out) :: Qt,Rt 
REAL(dp), DIMENSION(52,NCOND+1), INTENT(out) :: v, lnTj
REAL(dp), DIMENSION(52,52), INTENT(out) :: Ymn
INTEGER, DIMENSION(52) :: main_subgroup_index
REAL(dp), DIMENSION(52,NCOND+1) :: Xmj
REAL(dp), DIMENSION(52) :: sumQmj,Qmj
REAL(dp), DIMENSION(23,23) :: anm
INTEGER :: i,j,jj

v=0D0 ! number of subgroups of type t in component j (solvents)
v(28,1)=1; ! stoichiometric number subgroups for water
v(1:52,2:NCOND+1)=UNIFAC_groups! number of UNIFAC subgroups in organic compounds

! UNIFAC short range (SR) interactions of organics and water:
! Rt and Qt are the relative van der Waals subgroup volume and surface area parameters, respectively.
Rt=(/0.9011,& ! CH3, hexane, Hansen et al., 1991, EAIM
0.6744,& ! CH2, octane, Hansen et al., 1991, EAIM
0.4469,& ! CH, octane, Hansen et al., 1991, EAIM
0.2195,& ! C, octane, Hansen et al., 1991, EAIM
0.9011,& ! CH3[OH], [bonded to OH-group] Marcolli and Peter, 2005
0.6744,& ! CH2[OH], [bonded to OH-group] Marcolli and Peter, 2005
0.4469,& ! CH[OH], [bonded to OH-group] Marcolli and Peter, 2005
0.2195,& ! C[OH], [bonded to OH-group] Marcolli and Peter, 2005
0.9011,& ! CH3[alc-tail] Marcolli and Peter, 2005
0.6744,& ! CH2[alc-tail], [bonded to OH-group] Marcolli and Peter, 2005
0.4469,& ! CH[alc-tail], [bonded to OH-group] Marcolli and Peter, 2005
0.2195,& ! C[alc-tail], [bonded to OH-group] Marcolli and Peter, 2005
0.9011,& ! CH3[alc] Marcolli and Peter, 2005
0.6744,& ! CH2[alc], [bonded to OH-group] Marcolli and Peter, 2005
0.4469,& ! CH[alc], [bonded to OH-group] Marcolli and Peter, 2005
0.2195,& ! C[alc], [bonded to OH-group] Marcolli and Peter, 2005
1.3454,& ! CH2=CH, Hansen et al., 1991, EAIM
1.1167,& ! CH=CH, Hansen et al., 1991, EAIM
1.1173,& ! CH2=C, Hansen et al., 1991, EAIM
0.8886,& ! CH=C, Hansen et al., 1991, EAIM
0.6605,& ! C=C, Hansen et al., 1991, EAIM
0.5313,& ! ACH, Hansen et al., 1991, EAIM
0.3652,& ! AC, Hansen et al., 1991, EAIM
1.2663,& ! ACCH3, Hansen et al., 1991, EAIM
1.0396,& ! ACCH2, Hansen et al., 1991, EAIM
0.8121,& ! ACCH, Hansen et al., 1991, EAIM
1.00,& ! OH, Hansen et al., 1991, EAIM
0.92,& ! H2O, Hansen et al., 1991, EAIM
0.8952,& ! ACOH, Hansen et al., 1991, EAIM
1.6724,& ! CH3CO, Hansen et al., 1991, EAIM
1.4457,& ! CH2CO, Hansen et al., 1991, EAIM
0.998,& ! HCO, Hansen et al., 1991, EAIM
1.9031,& ! CH3COO, Hansen et al., 1991, EAIM
1.6764,& ! CH2COO, Hansen et al., 1991, EAIM
1.242,& ! HCOO, Hansen et al., 1991, EAIM
1.145,& ! CH3O, Hansen et al., 1991, EAIM
0.9183,& ! CH2O, Hansen et al., 1991, EAIM
0.6908,& ! CHO, Hansen et al., 1991, EAIM
1.3013,& ! COOH, Hansen et al., 1991, EAIM
2.0086,& ! CH3NO2, Hansen et al., 1991, EAIM
1.7818,& ! CH2NO2, Hansen et al., 1991, EAIM
1.5544,& ! CHNO2, Hansen et al., 1991, EAIM
1.4199,& ! ACNO2, Hansen et al., 1991, EAIM
2.1246,& ! CH2ONO2***, Compernolle et al., 2009*** 
1.8971,& ! CHONO2***
1.6697,& ! CONO2***
1.5869,& ! CH2OOH, hydrogen peroxide***
1.3594,& ! CHOOH, hydrogen peroxide***
1.132,& ! COOH, hydrogen peroxide***
1.7025,& ! C(=O)OOH***
2.0,& ! CHnOOCHm, peroxide***
2.6217/) ! PAN*** 

Qt=(/0.848,& ! CH3, hexane, Hansen et al., 1991, EAIM
0.540,& !CH2, hexane, Hansen et al., 1991, EAIM
0.228,& ! CH, hexane, Hansen et al., 1991, EAIM
0.0,& !C, hexane, Hansen et al., 1991, EAIM
0.848,& ! CH3[OH],[bonded to OH-group] Marcolli and Peter, 2005
0.540,& ! CH2[OH],[bonded to OH-group] Marcolli and Peter, 2005
0.228,& ! CH[OH],[bonded to OH-group] Marcolli and Peter, 2005
0.0,& ! C[OH],[bonded to OH-group] Marcolli and Peter, 2005
0.848,& ! CH3[alc-tail], Marcolli and Peter, 2005
0.540,& ! CH2[alc-tail], Marcolli and Peter, 2005
0.228,& ! CH[alc-tail], Marcolli and Peter, 2005
0.0,& ! C[alc-tail], Marcolli and Peter, 2005
0.848,& ! CH3[alc], Marcolli and Peter, 2005
0.540,& ! CH2[alc], Marcolli and Peter, 2005
0.228,& ! CH[alc], Marcolli and Peter, 2005
0.0,& ! C[alc], Marcolli and Peter, 2005
1.176,& ! CH2=CH, Hansen et al., 1991, EAIM
0.867,& ! CH=CH, Hansen et al., 1991, EAIM
0.988,& ! CH2=C, Hansen et al., 1991, EAIM
0.676,& ! CH=C, Hansen et al., 1991, EAIM
0.485,& ! C=C, Hansen et al., 1991, EAIM
0.4,& ! ACH, Hansen et al., 1991, EAIM
0.12,& ! AC, Hansen et al., 1991, EAIM
0.968,& ! ACCH3, Hansen et al., 1991, EAIM
0.660,& ! ACCH2, Hansen et al., 1991, EAIM
0.348,& ! ACCH, Hansen et al., 1991, EAIM
1.2,& ! OH, Hansen et al., 1991, EAIM
1.4,& ! H2O, Hansen et al., 1991, EAIM
0.68,& ! ACOH, Hansen et al., 1991, EAIM
1.448,& ! CH3CO, Hansen et al., 1991, EAIM
1.18,& ! CH2CO, Hansen et al., 1991, EAIM
0.948,& ! HCO, Hansen et al., 1991, EAIM
1.728,& ! CH3COO, Hansen et al., 1991, EAIM
1.42,& ! CH2COO, Hansen et al., 1991, EAIM
1.188,& ! HCOO, Hansen et al., 1991, EAIM
1.088,& ! CH3O, Hansen et al., 1991, EAIM
0.78,& ! CH2O, Hansen et al., 1991, EAIM
0.468,& ! CHO, Hansen et al., 1991, EAIM
1.224,& ! COOH, Hansen et al., 1991, EAIM
1.868,& ! CH3NO2, Hansen et al., 1991, EAIM
1.56,& ! CH2NO2, Hansen et al., 1991, EAIM
1.248,& ! CHNO2, Hansen et al., 1991, EAIM
1.104,& ! ACNO2, Hansen et al., 1991, EAIM
1.8682,& ! CH2ONO2***, Compernolle et al., 2009*** 
1.5562,& ! CHONO2***
1.3282,& ! CONO2***
1.437,& ! CH2OOH***
1.125,& ! CHOOH, hydrogen peroxide***
0.897,& ! COOH, hydrogen peroxide***
1.5217,& ! C(=O)OOH***
1.0,& ! CHnOOCHm, peroxide***
2.2887/) ! PAN***  

main_subgroup_index=(/1,& ! CHn
1,& ! CHn 
1,& ! CHn
1,& ! CHn
2,& ! CHn[OH]
2,& ! CHn[OH]
2,& ! CHn[OH]
2,& ! CHn[OH]
3,& ! CHn[alc-tail]
3,& ! CHn[alc-tail]
3,& ! CHn[alc-tail]
3,& ! CHn[alc-tail]
4,& ! CHn[alc]
4,& ! CHn[alc]
4,& ! CHn[alc]
4,& ! CHn[alc]
5,& ! C=C
5,& ! C=C
5,& ! C=C
5,& ! C=C
5,& ! C=C
6,& ! ACH
6,& ! ACH
7,& ! ACCH2
7,& ! ACCH2
7,& ! ACCH2
8,& ! OH
9,& ! H2O
10,& ! ACOH
11,& ! CH2CO
11,& ! CH2CO
12,& ! CHO
13,& ! CCOO
13,& ! CCOO
14,& ! HCOO
15,& ! CH2O
15,& ! CH2O
15,& ! CH2O
16,& ! COOH
17,& ! CNO2
17,& ! CNO2
17,& ! CNO2
18,& ! ACNO2
19,& ! CHnONO2, nitrate
19,& ! CHnONO2, nitrate
19,& ! CHnONO2, nitrate
20,& ! hydrogen peroxide
20,& ! hydrogen peroxide
20,& ! hydrogen peroxide
21,& ! peroxyacid
22,& ! peroxide
23/) ! PAN

! Main groups: CHn, CHn[OH], CHn[alc-tail], CHn[alc], C=C, AromaticC, AromC-alkane, OH, H2O, AromaticC-OH, Carbonyl(keton), Aldehyde Acetate Formate ether, COOH, CNO2, ACNO2, nitrate***, hydroperoxide***, peroxyacid***, peroxide***, PAN***
! anm from EAIM, origionaly from Hansen et al., 1991, Peng et al., 2001,
! Marcolli and Peter, 2005 and Compernolle et al., 2009*** 

anm(1,:)=(/0.0, 0.0, 0.0, 0.0, 245.21, 61.13, 76.5, 986.5, 1318., 1333., 476.4, 677., &
232.1, 507., 251.5, 663.5, 661.5, 543., 500.95, 977.56, 1331., 297.24, 528.5/)
anm(2,:)=(/0.0, 0.0, 0.0, 0.0, 245.21, 61.13, 76.5, 986.5, 2314., 1333., 476.4, 677., &
232.1, 507., 251.5, 663.5, 661.5, 543., 500.95, 977.56, 1331., 297.24, 528.5/)
anm(3,:)=(/0.0, 0.0, 0.0, 0.0, 245.21, 61.13, 76.5, 986.5, 1325., 1333., 476.4, 677., &
232.1, 507., 251.5, 663.5, 661.5, 543., 500.95, 977.56, 1331., 297.24, 528.5/)
anm(4,:)=(/0.0, 0.0, 0.0, 0.0, 245.21, 61.13, 76.5, 986.5, 1890., 1333., 476.4, 677., &
232.1, 507., 251.5, 663.5, 661.5, 543., 500.95, 977.56, 1331., 297.24, 528.5/)
anm(5,:)=(/-35.360, -35.360, -35.360, -35.360, 0.0, 38.81, 74.15, 524.1, 270.6, 526.1, &
182.6, 448.8, 37.85, 333.5, 214.5, 318.9, 357.5, 0.0, 10326., 475.91, 742.38, 606.71, 469.27/)
anm(6,:)=(/-11.120, -11.120, -11.120, -11.120, 3.446, 0.0, 167., 636.1, 903.8, 1329., &
25.77, 347.3, 5.994, 287.1, 32.14, 537.4, 168., 194.9, 500.95, 977.56, 1331., 297.24, 528.5/)
anm(7,:)=(/-69.700, -69.700, -69.700, -69.700, -113.6, -146.8, 0.0, 803.2, 5695., 884.9, &
-52.1, 586.6, 5688., 197.8, 213.1, 872.3, 3629., 4448., 500.95, 977.56, 1331., 297.24, 528.5/)
anm(8,:)=(/156.40, 156.40, 156.40, 156.40, 457., 89.6, 25.82, 0.0, 276.4, -259.7, 84., &
-203.6, 101.1, 267.8, 28.06, 224.39, 256.5, 157.1, 37.631, -330.28, 1789., 221.38, -77.526/)
anm(9,:)=(/300.00, -89.71, 362.1, 162.3, 496.1, 362.3, 377.6, -153., 0.0, 324.5, -195.4, &
-116., 72.87, 233.87, 540.5, -69.29, 220.6, 399.5, 142.65, -341.18, -329.81, -7.2937, 76.211/)
anm(10,:)=(/275.8, 275.8, 275.8, 275.8, 217.5, 25.34, 244.2, -451.6, -601.8, 0.0, -356.1, &
-271.1, -449.4, -32.52, -162.9, 408.9, 0.0, -413.48, 37.631, -330.28, 1789., 221.38, -77.526/)
anm(11,:)=(/26.76, 26.76, 26.76, 26.76, 42.92, 140.1, 365.8, 164.5, 472.5, -133.1, 0.0, &
-37.36, -213.7, -190.4, -103.6, 669.4, 137.5, 548.5, -197.93, -350.58, 252.05, -286.39, -3.8839/)
anm(12,:)=(/505.7, 505.7, 505.7, 505.7, 56.3, 23.39, 106., 529., 480.8, -155.6, 128., 0.0, &
-101.3, 766., 304.1, 497.5, 0.0, 0.0, 402., -387.63, 12274., -18.524, 308.97/)
anm(13,:)=(/114.8, 114.8, 114.8, 114.8, 132.1, 85.84, -170., 245.4, 200.8, -36.72, 372.2, &
-185.1, 0.0, -241.8, -235.7, 660.2, -81.13, 0.0, 1273.8, 928.33, 416., -252.22, 426.52/)
anm(14,:)=(/329.3, 329.3, 329.3, 329.3, 110.4, 18.12, 428., 139.4, 124.63, -234.25, 385.4, &
-236.5, 1167., 0.0, -234., -268.1, 0.0, 0.0, 1273.8, 928.33, 416., -252.22, 426.52/)
anm(15,:)=(/83.36, 83.36, 83.36, 83.36, 26.51, 52.13, 65.69, 237.7, -314.7, -178.5, 191., &
-7.838, 461.3, 457.3, 0.0, 664.6, 95.18, 155.11, 1133.1, -438.74, 2221.9, -130.54, 1160.7/)
anm(16,:)=(/315.3, 315.3, 315.3, 315.3, 1264., 62.32, 89.86, -103.03, -145.88, -11., -297.8, &
-165.5, -256.3, 193.9, -338.5, 0.0, 0.0, 0.0, -100.17, -501.23, -579.8, 79.052, -340.95/)
anm(17,:)=(/-32.69, -32.69, -32.69, -32.69, -1.996, 10.38, -97.05, 261.6, 417.9, 0.0, -142.6, &
0.0, 129.3, 0.0, -94.49, 0.0, 0.0, 533.2, 0.0, 0.0, 0.0, 0.0, 0.0/)
anm(18,:)=(/5541., 5541., 5541., 5541., 16.62, 1824., -127.8, 561.6, 360.7, 815.12, -101.5, &
0.0, 0.0, 0.0, 220.66, 0.0, -85.12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
anm(19,:)=(/-75.718, -75.718, -75.718, -75.718, -294.43, -75.718, -75.718, 818.97, 681.78, &
818.97, 188.72, -179.38, -356.25, -356.25, -289.81, 1173.3, 0.0, 0.0, 0.0, 545.66, 551.95, &
-308.16, -239.65/)
anm(20,:)=(/-23.233, -23.233, -23.233, -23.233, -57.949, -23.233, -23.233, 342.92, 795.55, &
342.92, 380.94, 408.88, -355., -355., 490.36, 1479., 0.0, 0.0, -86.279, 0.0, 202.91, -395.81, -147.47/)
anm(21,:)=(/5853.1, 5853.1, 5853.1, 5853.1, 883.78, 5853.1, 5853.1, -457.93, 670.32, -457.93, &
-98.45, -520.9, 131.15, 131.15, -471.67, 1896.1, 0.0, 0.0, 221.82, -62.0167, 0.0, 210.57, 395.33/)
anm(22,:)=(/-151.61, -151.61, -151.61, -151.61, -237.61, -151.61, -151.61, 820.86, 483.553, &
820.86, 587.21, 509.17, 449.04, 449.04, 142.65, 1043.9, 0.0, 0.0, 676.62, 1088.8, 537.7, 0.0, -2.0795/)
anm(23,:)=(/333.07, 333.07, 333.07, 333.07, 86.307, 333.07, 333.07, 612.05, 319.99, 612.05, &
111.76, -187.02, -157.64, -157.64, -208.91, 1207.7, 0.0, 0.0, 474.47, 392.54, -80.543, 339.08, 0.0/)

DO j=1,52
Xmj(j,:)=v(j,:)/SUM(v,DIM=1) ! mole fraction of subgroups in each particle compound (j)
END DO
  
Ymn=EXP(-anm(main_subgroup_index,main_subgroup_index)/T) ! temperature dependant function of the subgroup interaction parameters (anm)

lnTj=0D0
sumQmj=0D0

DO j=1,NCOND+1
Qmj=Qt*Xmj(:,j)/SUM(Qt*Xmj(:,j)) ! Relative surface area fraction of subgroup m
  DO i=1,52
     DO jj=1,52
     sumQmj(jj)=Qmj(jj)*Ymn(i,jj)/SUM(Qmj*Ymn(:,jj))
     END DO
  lnTj(i,j)=Qt(i)*(1D0-LOG(SUM(Qmj*Ymn(:,i)))-SUM(sumQmj)) 
  END DO
END DO

END SUBROUTINE UNIFAC_parameters

END MODULE thermodynamicsDMS
