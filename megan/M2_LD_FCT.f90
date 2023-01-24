!=======================================================================
!  LD_FCT.EXT
!  This include file contains "light dependent" factors.
!
!
!
!  This is an input file for MEGAN v2.0
!  This is version v210 of this file
!  Created by A. Guenther 8/11/07
!=======================================================================

Module M2_LD_FCT

      INTEGER, parameter ::   N_LDF_SPC=20
      CHARACTER*16       ::   LDF_SPC(20)
      REAL               ::   LDF_FCT(20)
      INTEGER            ::   LDF_MAP(20)

      DATA     LDF_SPC(  1)      , LDF_FCT(  1), LDF_MAP(  1) &
            / 'ISOP            ', 0.9999      , 1            /
            
      DATA     LDF_SPC( 2)      , LDF_FCT( 2), LDF_MAP( 2) &
            / 'MBO             ', 0.9999      , 2           /      

      DATA     LDF_SPC(  3)      , LDF_FCT(  3), LDF_MAP(  3) &
            / 'MYRC            ', 0.6        , 3            /

      DATA     LDF_SPC(  4)      , LDF_FCT(  4), LDF_MAP(  4) &
            / 'SABI            ', 0.6         , 4            /

      DATA     LDF_SPC(  5)      , LDF_FCT(  5), LDF_MAP(  5) &
            / 'LIMO            ', 0.2        , 5            /

      DATA     LDF_SPC(  6)      , LDF_FCT(  6), LDF_MAP(  6) &
            / '3CAR            ', 0.2        , 6            /

      DATA     LDF_SPC(  7)      , LDF_FCT(  7), LDF_MAP(  7) &
            / 'OCIM            ', 0.8         , 7            /

      DATA     LDF_SPC(  8)      , LDF_FCT(  8), LDF_MAP(  8) &
            / 'BPIN            ', 0.2         , 8            /

      DATA     LDF_SPC(  9)      , LDF_FCT(  9), LDF_MAP(  9) &
            / 'APIN            ', 0.6         , 9            /

      DATA     LDF_SPC( 10)      , LDF_FCT( 10), LDF_MAP( 10) &
            / 'FARN            ', 0.5         , 10           /

      DATA     LDF_SPC( 11)      , LDF_FCT( 11), LDF_MAP( 11) &
            / 'BCAR            ', 0.5         , 11           /
    
      DATA     LDF_SPC( 12)      , LDF_FCT( 12), LDF_MAP( 12) &
            / 'MEOH            ', 0.8        , 12           /

      DATA     LDF_SPC( 13)      , LDF_FCT( 13), LDF_MAP( 13) &
            / 'ACTO            ', 0.2        , 13           /
            
      DATA     LDF_SPC( 14)      , LDF_FCT( 14), LDF_MAP( 14) &
            / 'ACTA            ', 0.2         , 14           /

      DATA     LDF_SPC( 15)      , LDF_FCT( 15), LDF_MAP( 15) &
            / 'FORM            ', 0.5         , 15           /      

      DATA     LDF_SPC( 16)      , LDF_FCT( 16), LDF_MAP( 16) &
            / 'CH4             ', 0.75        , 16           /

      DATA     LDF_SPC( 17)      , LDF_FCT( 17), LDF_MAP( 17) &
            / 'NO              ', 0.0         , 17           /
      
      DATA     LDF_SPC( 18)      , LDF_FCT( 18), LDF_MAP( 18) &
            / 'OMTP            ', 0.6         , 18            /
      
      DATA     LDF_SPC( 19)      , LDF_FCT( 19), LDF_MAP( 19) &
            / 'OSQT            ', 0.5         , 19           /
    
      DATA     LDF_SPC( 20)      , LDF_FCT( 20), LDF_MAP( 20) &
             / 'CO              ', 0.9999         , 20           /

end module M2_LD_FCT
