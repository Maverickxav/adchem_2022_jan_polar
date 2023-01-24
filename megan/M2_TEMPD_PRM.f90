!=======================================================================
!  TEMPD_PRM.EXT
!  This include file contains "temperature dependent" parameter for 
!  light-independent emissions 
!
!
!  This is an input file for MEGAN v2.0
!  This is version v210 of this file
!  Created by A. Guenther 8/11/07
!=======================================================================

MODULE M2_TEMPD_PRM

      INTEGER, parameter     ::   N_TDF_SPC=20
      CHARACTER*16 ::  TDF_SPC(20)
      REAL        ::   TDF_PRM(20)
      INTEGER     ::   TDF_MAP(20)

      DATA     TDF_SPC(  1)      , TDF_PRM(  1), TDF_MAP(  1) & 
            / 'ISOP            ', 0.09        , 1           /
      DATA     TDF_SPC( 2)      , TDF_PRM( 2), TDF_MAP( 2) & 
            / 'MBO             ', 0.09        , 2          /      

      DATA     TDF_SPC(  3)      , TDF_PRM(  3), TDF_MAP(  3) & 
            / 'MYRC            ', 0.15         , 3           /

      DATA     TDF_SPC(  4)      , TDF_PRM(  4), TDF_MAP(  4) &
            / 'SABI            ', 0.15         , 4           /

      DATA     TDF_SPC(  5)      , TDF_PRM(  5), TDF_MAP(  5) & 
            / 'LIMO            ', 0.15         , 5           /

      DATA     TDF_SPC(  6)      , TDF_PRM(  6), TDF_MAP(  6) & 
            / '3CAR            ', 0.15         , 6           /

      DATA     TDF_SPC(  7)      , TDF_PRM(  7), TDF_MAP(  7) & 
            / 'OCIM            ', 0.15         , 7           /

      DATA     TDF_SPC(  8)      , TDF_PRM(  8), TDF_MAP(  8) &
            / 'BPIN            ', 0.15         , 8           /

      DATA     TDF_SPC(  9)      , TDF_PRM(  9), TDF_MAP(  9) & 
            / 'APIN            ', 0.15         , 9           /

      DATA     TDF_SPC( 10)      , TDF_PRM( 10), TDF_MAP( 10) & 
            / 'FARN            ', 0.17        , 10          /

      DATA     TDF_SPC( 11)      , TDF_PRM( 11), TDF_MAP( 11) & 
            / 'BCAR            ', 0.17        , 11          /

      DATA     TDF_SPC( 12)      , TDF_PRM( 12), TDF_MAP( 12) & 
            / 'MEOH            ', 0.08        , 12          /

      DATA     TDF_SPC( 13)      , TDF_PRM( 13), TDF_MAP( 13) & 
            / 'ACTO            ', 0.11        , 15          /
      
      DATA     TDF_SPC( 14)      , TDF_PRM( 14), TDF_MAP( 14) &
            / 'ACTA            ', 0.13        , 14          /

      DATA     TDF_SPC( 15)      , TDF_PRM( 15), TDF_MAP( 15) &
            / 'FORM            ', 0.09        , 15          /

      DATA     TDF_SPC( 16)      , TDF_PRM( 16), TDF_MAP( 16) &
            / 'CH4             ', 0.05        , 16          /

      DATA     TDF_SPC( 17)      , TDF_PRM( 17), TDF_MAP( 17) &
            / 'NO              ', 0.11        , 17          /

      DATA     TDF_SPC( 18)      , TDF_PRM( 18), TDF_MAP( 18) & 
            / 'OMTP            ', 0.15         , 18           /
            
      DATA     TDF_SPC( 19)      , TDF_PRM( 19), TDF_MAP( 19) & 
            / 'OSQT            ', 0.17        , 19          /

      DATA     TDF_SPC( 20)      , TDF_PRM( 20), TDF_MAP( 20) &
            / 'CO              ', 0.09        , 20          /

END MODULE M2_TEMPD_PRM
