!=======================================================================
!  REL_EM_ACT.EXT
!  This include file contains "production and loss within canopy"
!  factors.
!
!
!  MEGAN v2.0
!
!  Created by Tan 11/30/06
!=======================================================================

MODULE M2_REL_EM_ACT

      INTEGER, parameter ::      N_CAT=5
      REAL    ::       Anew(5)
      REAL    ::       Agro(5)
      REAL    ::      Amat(5)
      REAL    ::       Aold(5)

      DATA    Anew(  1),  Agro(  1),  Amat(  1),  Aold(  1) &
          /  1.0      ,  1.0      ,  1.0      ,  1.0       /

      DATA    Anew(  2),  Agro(  2),  Amat(  2),  Aold(  2) &
          /  2.0      ,  1.8      ,  0.95     ,  1.0       /

      DATA    Anew(  3),  Agro(  3),  Amat(  3),  Aold(  3) &
          /  0.4      ,  0.6      ,  1.075    ,  1.0       /

      DATA    Anew(  4),  Agro(  4),  Amat(  4),  Aold(  4) &
          /  3.0      ,  2.6      ,  0.85     ,  1.0       /

      DATA    Anew(  5),  Agro(  5),  Amat(  5),  Aold(  5) &
          /  0.05     ,  0.6      ,  1.125    ,  1.0       /
     
END MODULE M2_REL_EM_ACT
