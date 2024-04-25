#ifdef TIDE
         subroutine equ_tide_mom

#include <def-undef.h>
         use precision_mod
         use param_mod
         use pconst_mod
         use dyn_mod
         implicit none

!                               w(1/day) 1 beta h a (m) 
!1 K1 Luni-solar declinational 0.7292117 0.736 0.141565 
!2 O1 Principal lunar declinational 0.6759774 0.695 0.100661 
!3 P1 Principal solar declinational 0.7252295 0.706 0.046848 
!4 Q1 Larger lunar elliptic 0.6495854 0.695 0.019273 
!5 M2 Principal lunar 1.405189 0.693 0.242334 
!6 S2 Principal solar 1.454441 0.693 0.112743 
!7 N2 Largerl lunar elliptic 1.378797 0.693 0.046397 
!8 K2 Luni-solar declinational 1.458423 0.693 0.030684 
         INTEGER, PARAMETER :: TN=8
         REAL(R8), PARAMETER,DIMENSION(TN) :: &
              amp = (/0.141565_r8, 0.100661_r8, 0.046848_r8, 0.019273_r8,&
                      0.242334_r8, 0.112743_r8,0.046397_r8,0.030684_r8/)
         REAL(R8), PARAMETER,DIMENSION(TN) :: &
           freq = (/0.7292117_r8,0.6759774_r8,0.7252295_r8,0.6495854_r8,&
                    1.405189_r8,1.454441_r8,1.378797_r8,1.458423_r8/)

         REAL(R8), PARAMETER,DIMENSION(TN) :: &
           beta = (/0.736_r8,0.695_r8,0.706_r8,0.695_r8,&
                    0.693_r8,0.693_r8,0.693_r8,0.693_r8/)

         REAL(R8) :: PHAI_IJ
!         REAL(R8),DIMENSION(IMT,JMT_global) :: elevat_io

         elevat=0.0_r8

         DO J=1,JMT
         DO I=1,IMT

            PHAI_IJ = lon_tidal(I)*2.*PI/360.0
         DO K=1,TN
          IF(K > 4)THEN
            elevat(i,j) = vit(i,j,1)*beta(k)*amp(k)*COS(lat_tidal(j)*2.*PI/360.0)**2  &
                    *COS(freq(k)*1e-4*time_tidal*IDTB+2.0*PHAI_IJ) + elevat(i,j)
          ELSE
            elevat(i,j) = vit(i,j,1)*beta(k)*amp(k)*SIN(2.0*lat_tidal(j)*2.*PI/360.0) &
                    *COS(freq(k)*1e-4*time_tidal*IDTB+PHAI_IJ) + elevat(i,j)
          END IF
         END DO


         END DO
         END DO

!        call local_to_global_4d(elevat,elevat_io,1,1)
!        if (mytid.eq.0) then
!        print*,elevat_io
!        stop 9988
!        endif

         end subroutine equ_tide_mom

#else
SUBROUTINE equ_tide_mom
RETURN    
END SUBROUTINE equ_tide_mom
#endif
