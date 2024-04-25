#ifdef TIDE
         subroutine equ_tide

#include <def-undef.h>
         use precision_mod
         use param_mod
         use pconst_mod
         use dyn_mod
         implicit none

         INTEGER, PARAMETER :: TN=6
!                          S2         M2       N2           K1          P1           O1      
         REAL(R8), PARAMETER,DIMENSION(TN) :: &
              amp = (/0.112841_r8, 0.242334_r8, 0.046398_r8, 0.141565_r8,&
                      0.046848_r8, 0.100514_r8/)
         REAL(R8), PARAMETER,DIMENSION(TN) :: &
           period = (/43200.0_r8, 44712.0_r8, 45570.0_r8, &
                      86164.0_r8, 86637.0_r8, 92950.0_r8/)
         REAL(R8), PARAMETER,DIMENSION(TN) :: &
           Klove = (/0.302_r8, 0.302_r8, 0.302_r8, 0.256_r8, 0.000_r8, 0.298_r8/)
         REAL(R8), PARAMETER,DIMENSION(TN) :: &
           Hlove = (/0.602_r8, 0.602_r8, 0.602_r8, 0.520_r8, 0.000_r8, 0.603_r8/)  

         REAL(R8), PARAMETER :: atm_amp=0.0113_r8,atm_period=43200.0_r8,atm_alpha=112.0_r8

         REAL(R8),DIMENSION(TN) :: freq,atm_freq,beta
         REAL(R8) :: PHAI_IJ,atm_alpha1
!         REAL(R8),DIMENSION(IMT,JMT_global) :: elevat_io

         DO I=1,TN
         beta(I) = 1.0_r8+KLOVE(I)-HLOVE(I)
         END DO

         DO I=1,TN
         freq(I) = 2.*PI/PERIOD(I)
         END DO

!         atm_freq = 2.*PI/atm_PERIOD
!         atm_alpha1  = atm_alpha*2.0*PI/360.0

         DO J=1,JMT
         DO I=1,IMT
         elevat(i,j)=0.0_r8
!         atm_elevat(i,j)=0.0_r8

           PHAI_IJ = lon_tidal(I)*2.*PI/360.0
         DO K=1,TN
          IF(K <= 3)THEN
            elevat(i,j) = vit(i,j,1)*beta(k)*amp(k)*COS(lat_tidal(j)*2.*PI/360.0)**2  &
                    *COS(freq(k)*time_tidal*IDTB+2.0*PHAI_IJ) + elevat(i,j)
          ELSE
            elevat(i,j) = vit(i,j,1)*beta(k)*amp(k)*SIN(2.0*lat_tidal(j)*2.*PI/360.0) &
                    *COS(freq(k)*time_tidal*IDTB+PHAI_IJ) + elevat(i,j)
          END IF
         END DO

!            atm_elevat(i,j) = atm_amp*COS(atm_freq*time_tidal*IDTB+2.0*PHAI_IJ-atm_alpha1)
!            elevat(i,j) = elevat(i,j)+atm_elevat(i,j)
            elevat(i,j) = elevat(i,j)

         END DO
         END DO

!        call local_to_global_4d(elevat,elevat_io,1,1)
!        if (mytid.eq.0) then
!        print*,elevat_io
!        stop 9988
!        endif

         end subroutine equ_tide
#else
SUBROUTINE equ_tide
RETURN
END SUBROUTINE equ_tide
#endif
