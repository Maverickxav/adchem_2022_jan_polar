MODULE acidityDMS
USE constants
USE second_Precision, ONLY : dp    ! KPP Numerical type
IMPLICIT NONE
PRIVATE 

PUBLIC :: ACIDITY_DMS

CONTAINS
SUBROUTINE ACIDITY_DMS(mNV,mNIII,mSVI,mClI,mNa,mCH3SO3I,T,yH,yHSO4,ySO4,yNH4,yNO3,&
yCl,pCO2,mH,fHSO4,fSO4,fNO3,fCl,fCH3SO3,fNH4,mHCO3,mCO3,mOH)

REAL(dp), INTENT(in) :: mNV,mNIII,mSVI,mClI,mNa,mCH3SO3I,T,yH,yHSO4,ySO4,&
yNH4,yNO3,yCl,pCO2
REAL(dp), INTENT(out) :: mH,fHSO4,fSO4,fNO3,fCl,fCH3SO3,fNH4,mHCO3,mCO3,mOH
REAL(dp) :: Kw,K1,K2,K4,K5,K6,H_NH3,KNH3,KHCl,KCH3SO3H,KHNO3,Ka,yHCO3,yCO3,yOH,yCOO,yCH3SO3,&
A,B,C,D,E,F,G,I,J,K,L,M,N,O,P,Q,R,S,T2,U,X,X2,X3,X4,X5,X6,Y,Y2,Y3,Y4,Y5,Y6,mHSO4,mSO4
REAL(dp), DIMENSION(9) :: alfa
!REAL(dp), DIMENSION(8) :: root
 
! Solve the ion ballance equation:
! 2014-01-26: Calculate fNO3,fCl and fNH4 and give as output and inpuit to the charge
! ballance equation for NH3/NH4 equilibration in the condensation solver.

! yH,yHSO4,ySO4 ... is the single solute activity coefficients of H+, HSO4- and SO4 2-


Kw=1D-14 ! H2O <-> OH- + H+ 
K1=1.015D-2*EXP(8.85*(298./T-1D0)+25.14*(1D0+LOG(298./T)-298./T)) ! equilibrium constant HSO4- <-> H+ + SO42-
K2=1.03D11*EXP(34.81*(298./T-1D0)-5.39*(1D0+LOG(298./T)-298./T)) ! equilibrium constant NH3(g)+H+ <->NH4+
K4=3.4D-2 ! atm^-1 Henrys constant for CO2
K5=4.3D-7 ! H2CO3<->HCO3- + H+
K6=4.7D-11 ! HCO3- <-> Co32- + H+

H_NH3=57.6*EXP(13.79*(298./T-1D0)-5.39*(1D0+LOG(298./T)-298./T)) ! mol/(kg atm) NH3(g)<->NH3(aq)
KNH3=K2/H_NH3 ! kg/mol NH3(aq)+H+<->NH4+
KHCl=1.72D6*EXP(23.15*(298./T-1D0)) ! mol/kg HCl(aq)<->H+ + Cl-
!KCH3SO3H=10**(-1.9) ! mol/kg CH3SO3H(aq)<->H+ + CH3SO3-
KCH3SO3H=10**(2.93) ! COSMOtherm
KHNO3=12.*EXP(29.17*(298./T-1D0)+16.83*(1D0+LOG(298./T)-298./T)) ! mol/kg HNO3(aq)<->H+ + NO3-

!Ka=10**(-4.6) ! COOH <-> H+ COO-  (pKa=4.6 for PINONIC and PINIC ACID)

yHCO3=1D0 
yCO3=1D0 
yOH=yH
yCOO=1D0
yCH3SO3=1.0

A= yOH*yH
B= K1*yHSO4*mSVI
C= ySO4*yH
D= K1*yHSO4
E= KNH3*mNIII
F= yNH4*yH
G= KNH3
I= KHCl*mClI
J= yH*yCl
K= KHCl
L= KHNO3*mNV
M= yH*yNO3
N= KHNO3
O= K5*K4*pCO2
P= yH*yHCO3
Q= K6*K5*K4*pCO2
R= yH**2D0*yCO3
S=mSVI
T2=KCH3SO3H*mCH3SO3I
U=yH*yCOO


Y=C*N+D*M
X=J*KCH3SO3H+K*U  
X2=M*K+N*J
Y2=U*G+KCH3SO3H*F
X3=C*K+D*J
Y3=U*G+KCH3SO3H*F
X4=C*N+D*M
Y4=U*G+KCH3SO3H*F
X5=C*N+D*M
Y5=J*G+K*F
X6=C*N+D*M
Y6=J*KCH3SO3H+K*U

  alfa(9)=-A*P*R*C*M*J*U*F                                                      ! a(1)*x^8
  alfa(8)=S*A*P*R * C*M*J*U*F - mNa*A*P*R*C*M*J*U*F &
   -A*P*R* (C*M*J*U*G + C*M*X*F + Y*J*U*F)                                      ! a(2)*x^7
  alfa(7)=Kw*P*R*C*M*J*U*F + O*A*R*C*M*J*U*F &
      +S*A*P*R * (C*M*J*U*G + C*M*X*F + Y*J*U*F) &
      +B*A*P*R*U*F*M*J + L*A*P*R*U*F*C*J &
      +I*A*P*R*U*F*C*M + T2*A*P*R*J*F*C*M &
      -mNa*A*P*R* (C*M*J*U*G + C*M* X*F + Y*J*U*F) - E*A*P*R*J*U*C*M &
      -A*P*R* ( C*M* X*G + C*M* F*K*KCH3SO3H + Y*J*U*G + Y* X*F + D*N*J*U*F)           ! a(3)*x^6
  alfa(6)=Kw*P*R*(C*M*J*U*G + C*M*X*F + Y*J*U*F) + &
      O*A*R*(C*M*J*U*G + C*M*X*F + Y*J*U*F) + &
      2D0*Q*A*P * C*M*J*U*F + &
      S*A*P*R *( C*M* X*G + C*M* F*K*KCH3SO3H + Y*J*U*G + Y* X*F + D*N*J*U*F) + &
      B*A*P*R*(M*J*Y2 + X2*U*F) + L*A*P*R*(C*J*Y3 + X3*U*F) + &
      I*A*P*R*(C*M*Y4 + X4*U*F) + T2*A*P*R*(C*M*Y5 + X5*J*F) - &
      mNa*A*P*R* ( C*M*X*G + C*M*F*K*KCH3SO3H + Y*J*U*G + Y*X*F + D*N*J*U*F) - &
      E*A*P*R*(C*M*Y6 + X6*J*U) - A*P*R* ( C*M*K*KCH3SO3H*G + &
      Y* X*G + Y*F*K*KCH3SO3H + D*N*J*U*G + D*N*X*F)                                  ! a(4)*x^5
  alfa(5)=Kw*P*R*( C*M*X*G + C*M*F*K*KCH3SO3H + Y*J*U*G + Y*X*F + D*N*J*U*F) + &  
      O*A*R*(C*M*X*G + C*M*F*K*KCH3SO3H + Y*J*U*G + Y*X*F + D*N*J*U*F) + &
      2D0*Q*A*P*(C*M*J*U*G + C*M*X*F + Y*J*U*F) + &
      S*A*P*R*(C*M*K*KCH3SO3H*G + Y*X*G + Y*F*K*KCH3SO3H + D*N*J*U*G + D*N*X*F) + &
      B*A*P*R*(M*J*KCH3SO3H*G + X2*Y2 + N*K*U*F) + L*A*P*R*(C*J*KCH3SO3H*G + X3*Y3 + D*K*U*F) + &
      I*A*P*R*(C*M*KCH3SO3H*G + X4*Y4 + D*N*U*F) + T2*A*P*R*(C*M*K*G + X5*Y5 + D*N*J*F) - &
      mNa*A*P*R* (C*M*K*KCH3SO3H*G + Y*X*G + Y*F*K*KCH3SO3H + D*N*J*U*G + D*N*X*F) - &
      E*A*P*R*(C*M*K*KCH3SO3H + X6*Y6 + D*N*J*U) - &
      A*P*R* (Y*K*KCH3SO3H*G + D*N*X*G + D*N*F*K*KCH3SO3H)                                  ! a(5)*x^4
  alfa(4)=Kw*P*R*(C*M*K*KCH3SO3H*G + Y*X*G + Y*F*K*KCH3SO3H + D*N*J*U*G + D*N*X*F) + &
      O*A*R*(C*M*K*KCH3SO3H*G + Y*X*G + Y*F*K*KCH3SO3H + D*N*J*U*G + D*N*X*F) + &
      2D0*Q*A*P*(C*M*X*G + C*M*F*K*KCH3SO3H + Y*J*U*G + Y*X*F + D*N*J*U*F) + &
      S*A*P*R*(Y*K*KCH3SO3H*G + D*N*X*G + D*N*F*K*KCH3SO3H) + &
      B*A*P*R*(X2*KCH3SO3H*G + N*K*Y2) + L*A*P*R*(X3*KCH3SO3H*G + D*K*Y3) + &
      I*A*P*R*(X4*KCH3SO3H*G + D*N*Y4) + T2*A*P*R*(X5*K*G + D*N*Y5) - &
      mNa*A*P*R* (Y*K*KCH3SO3H*G + D*N*X*G + D*N*F*K*KCH3SO3H) - &
      E*A*P*R*(X6*K*KCH3SO3H + D*N*Y6) - A*P*R*D*N*K*KCH3SO3H*G                             ! a(6)*x^3
  alfa(3)=Kw*P*R*( Y*K*KCH3SO3H*G + D*N*X*G + D*N*F*K*KCH3SO3H) + &
  O*A*R*(Y*K*KCH3SO3H*G + D*N*X*G + D*N*F*K*KCH3SO3H) + &
      2D0*Q*A*P*(C*M*K*KCH3SO3H*G + Y*X*G + Y*F*K*KCH3SO3H + D*N*J*U*G + D*N*X*F) + &
      S*A*P*R*D*N*K*KCH3SO3H*G + B*A*P*R*N*K*KCH3SO3H*G + L*A*P*R*D*K*KCH3SO3H*G + I*A*P*R*D*N*KCH3SO3H*G + &
      T2*A*P*R*D*N*K*G - mNa*A*P*R*D*N*K*KCH3SO3H*G - E*A*P*R*D*N*K*KCH3SO3H                ! a(7)*x^2
  alfa(2)=Kw*P*R*D*N*K*KCH3SO3H*G + O*A*R*D*N*K*KCH3SO3H*G + &
      2D0*Q*A*P * ( Y*K*KCH3SO3H*G + D*N*X*G + D*N*F*K*KCH3SO3H)                            ! a(8)*x^1
  alfa(1)=2D0*Q*A*P*D*N*K*KCH3SO3H*G                                                  ! a(9)

CALL root_finder(alfa,mH)

mHCO3=K5*K4*pCO2/mH  ! seinfekd and pandis EQ. 7.2
mCO3=K6*K5*K4*pCO2/(mH**2D0)
mSO4=K1*yHSO4*mSVI/(K1*yHSO4+mH*yH*ySO4)
mHSO4=mSVI-mSO4

IF (mSVI>0) THEN
fSO4=mSO4/mSVI
fHSO4=mHSO4/mSVI
ELSE
fSO4=1D0
fHSO4=0D0
END IF

fNO3=KHNO3/(mH*yH*yNO3+KHNO3) ! [NO3-]/([HNO3]+[NO3-]), Assuming yHNO3=1
fCl=KHCl/(mH*yH*yCl+KHCl) ! [NO3-]/([HNO3]+[NO3-]), Assuming yHCl=1
fNH4=KNH3*mH*yH/(KNH3*mH*yH+yNH4) ! [NH4+]/([NH3]+[NH4+]), Assuming yNH3=1
fCH3SO3=KCH3SO3H/(mH*yH*yCH3SO3+KCH3SO3H) ! [CH3SO3-]/([CH3SO3H]+[CH3SO3-]), Assuming yCH3SO3=1

mOH=Kw/mH
!mCOO=Ka*mCOOHtot/(mH*yH*yCOO+Ka)

END SUBROUTINE ACIDITY_DMS


SUBROUTINE root_finder(a,mH)
! NUMERICAL RECIPES IN FORTRAN 77: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43064-X)
! Copyright (C) 1986-1992 by Cambridge University Press.
! Programs Copyright (C) 1986-1992 by Numerical Recipes Software.

INTEGER :: m,j
REAL(dp), DIMENSION(9), INTENT(in) :: a ! Coefficients of polynomial 
COMPLEX(dp), DIMENSION(8) :: roots
REAL(dp), DIMENSION(8) :: real_roots
REAL(dp), INTENT(out) :: mH 
LOGICAL :: polish

polish = .TRUE.
m=8

CALL zroots(a,m,roots,polish)

DO j=1,m
IF (aimag(roots(j))==0) THEN
real_roots(j)=roots(j)
ELSE
real_roots(j)=0D0
END IF
END DO
mH = MAXVAL(real_roots) 

END SUBROUTINE root_finder

SUBROUTINE zroots(a,m,roots,polish)
INTEGER, INTENT(in) :: m
INTEGER, PARAMETER :: MAXM = 101 ! maximum anticipated value of m+1
REAL(dp), INTENT(in) :: a(m+1) 
COMPLEX(dp), INTENT(out) :: roots(m)
LOGICAL, INTENT(in) :: polish
REAL(dp), PARAMETER :: EPS=1.e-6 ! A small number.

!USES laguer
! Given the degree m and the complex coeffcients a(1:m+1) of the polynomial
!SUM m+1 i=1 a(i)xi−1, this routine successively calls laguer and finds all m complex roots. The logical variable
! polish should be input as .true. if polishing (also by Laguerre's method) is desired,
! .false. if the roots will be subsequently polished by other means.
INTEGER :: i,j,jj,its
COMPLEX(dp) :: x,b,c
REAL(dp) :: ad(MAXM)

DO j=1,m+1 ! Copy of coeffcients for successive deflation.
ad(j)=a(j)
END DO

DO j=m,1,-1 ! Loop over each root to be found.
x=cmplx(0.,0.) ! Start at zero to favour convergence to smallest remaining root.
CALL laguer(ad,j,x,its) ! Find the root.

IF(abs(aimag(x)).le.2.*EPS**2*abs(real(x))) x=cmplx(real(x),0.)
		
		
roots(j)=x
b=ad(j+1) ! Forward deflation.

DO  jj=j,1,-1
c=ad(jj)
ad(jj)=b
b=x*b+c
END DO
END DO

IF (polish) THEN
DO j=1,m ! Polish the roots using the undeflated coffcients.
CALL laguer(a,m,roots(j),its)
END DO
END IF

DO j=2,m ! Sort roots by their real parts by straight insertion.
x=roots(j)
DO i=j-1,1,-1
IF (real(roots(i)).le.real(x)) goto 10
roots(i+1)=roots(i)
END DO
i=0
10 roots(i+1)=x
END DO
RETURN
END SUBROUTINE zroots

SUBROUTINE laguer(a,m,x,its)

INTEGER, INTENT(in) ::  m
INTEGER, INTENT(out) :: its
INTEGER, PARAMETER :: MR=8, MT=10, MAXIT=MT*MR
REAL(dp), PARAMETER :: EPSS=2D-7
REAL(dp), INTENT(in) :: a(m+1)
COMPLEX(dp), INTENT(inout) :: x

! Given the degree m and the complex coeffcients a(1:m+1) of the polynomial SUM m+1 i=1 a(i)xi−1,
! and given a complex value x, this routine improves x by Laguerre's method until it converges,
! within the achievable rounoff limit, to a root of the given polynomial. The number
! of iterations taken is returned as its.
! Parameters: EPSS is the estimated fractional roundoff error. We try to break (rare) limit
!cycles with MR different fractional values, once every MT steps, for MAXIT total allowed
!iterations.

INTEGER :: iter,j
REAL(dp) :: abx,abp,abm,err,frac(MR)
COMPLEX(dp) :: dx,x1,b,d,f,g,h,sq,gp,gm,g2

SAVE frac
DATA frac /.5,.25,.75,.13,.38,.62,.88,1./ ! Fractions used to break a limit cycle.
	
DO iter=1,MAXIT  ! Loop over iterations up to allowed maximum.
its=iter
b=a(m+1)
err=abs(b)
d=cmplx(0.,0.)
f=cmplx(0.,0.)
abx=abs(x)

DO j=m,1,-1 ! Efficient computation of the polynomial and its first two derivatives.
f=x*f+d 
d=x*d+b
b=x*b+a(j)
err=abs(b)+abx*err
END DO

err=EPSS*err ! Estimate of roundoff error in evaluating polynomial.
if(abs(b).le.err) THEN ! We are on the root.
RETURN

ELSE ! The generic case: use Laguerre's formula.
g=d/b
g2=g*g
h=g2-2.*f/b
sq=sqrt((m-1)*(m*h-g2))
gp=g+sq
gm=g-sq
abp=abs(gp)
abm=abs(gm)

IF(abp.lt.abm) gp=gm
IF (max(abp,abm).gt.0.) THEN
dx=m/gp
ELSE
dx=exp(cmplx(log(1.+abx),float(iter)))
END IF
END IF

x1=x-dx
IF(x.eq.x1) RETURN ! Converged.

IF (mod(iter,MT).ne.0) THEN
x=x1
ELSE ! Every so often we take a fractional step, to break any limit cycle (itself a rare occurrence).
x=x-dx*frac(iter/MT)
END IF

END DO
write(*,*) m
write(*,*) 'too many iterations in laguer' ! Very unusual, can occur only for complex roots.
RETURN ! Try a different starting guess for the root.
END SUBROUTINE laguer

END MODULE acidityDMS
