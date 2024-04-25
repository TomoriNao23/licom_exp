#define LOGMSG()
!write(mytid+600,'(a,3i4)')"SMUV",__LINE__,k,j
!  CVS: $Id: smuvh.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     ========================
      SUBROUTINE SMUV (X,Z,KK,fil_lat)
!     ========================
!     1-D zonal smoother

#include <def-undef.h>
use precision_mod
use param_mod
#ifdef SPMD
use msg_mod
#endif
use pconst_mod, only: sinu,pi,torad
      IMPLICIT NONE
      INTEGER :: JFS1,JFS2,JFN1,JFN2,KK, NN,NCY
      REAL(r8)    :: fil_lat
      REAL(r8)    :: X (IMT,JMT,KK),XS (IMT),Z (IMT,JMT,KM)


!lhl      fil_lat=55.D0
!
!$OMP PARALLEL DO PRIVATE (K,J,nn,xs)
      DO K = 1,KK
         do j =jst , jmt
            if (sinu(j).le.cos(fil_lat*torad)) then
               nn=int(cos(fil_lat*torad)/sinu(j))*2

               DO NCY = 1,NN
                  DO I = 1,IMT
                     XS (I)= X (I,J,K)* Z (I,J,K)
                  END DO
               DO I = 2,IMM
                  X (I,J,K)= (0.5D0*XS(I)+0.25D0*(XS(I-1)+XS(I+1)))*Z(I,J,K)
               END DO
                  X (1,J,K)= X (IMM,J,K)
                  X (IMT,J,K)= X (2,J,K)
               END DO
!
            end if
         END DO
     END DO

      RETURN
      END SUBROUTINE SMUV



!     ========================
      SUBROUTINE SMTS (X,Z,KK,fil_lat)
!     ========================
!     1-D zonal smoother

#include <def-undef.h>
use precision_mod
use param_mod
#ifdef SPMD
use msg_mod
#endif
use pconst_mod,only: sint,PI,torad
      IMPLICIT NONE

      INTEGER :: JFS1,JFS2,JFN1,JFN2,KK, NN,NCY
      REAL(r8)    :: fil_lat
      REAL(r8)    :: X (IMT,JMT,KK),XS (IMT),Z (IMT,JMT,KM)


!lhl      fil_lat=56.D0
!
!$OMP PARALLEL DO PRIVATE (K,J,nn,xs)
      DO K = 1,KK
         do j =jst , jmt
            if (sint(j).le.cos(fil_lat*torad))then
               nn=int(cos(fil_lat*torad)/sint(j))*2

               DO NCY = 1,NN
                  DO I = 1,IMT
                     XS (I)= X (I,J,K)* Z (I,J,K)
                  END DO
               DO I = 2,IMM
                  X(I,J,K)=(XS(I)*(1.0D0-0.25D0*Z(I-1,J,K)-0.25D0*Z(I+1,J,K)) &
                            +0.25D0*(XS(I-1)+XS(I+1)))*Z(I,J,K)
               END DO
                  X (1,J,K)= X (IMM,J,K)
                  X (IMT,J,K)= X (2,J,K)
               END DO
!
            end if
         END DO
     END DO


      RETURN
      END SUBROUTINE SMTS



!     ========================
      SUBROUTINE SMZ0 (X,Z,fil_lat)
!     ========================
!     1-D zonal smoother

#include <def-undef.h>
use precision_mod
      use param_mod
#ifdef SPMD
use msg_mod
#endif
use pconst_mod,only:sint,pi,torad
      IMPLICIT NONE

      INTEGER :: JFS1,JFS2,JFN1,JFN2,NN, NCY
      REAL(r8)    :: fil_lat
      REAL(r8)    :: X (IMT,JMT),XS (IMT),Z (IMT,JMT,KM)


!lhl      fil_lat=54.D0
!
!$OMP PARALLEL DO PRIVATE (K,J,nn,xs)
         do j =jst , jmt
            if (sint(j).le.cos((fil_lat-1.0D0)*torad))then
               nn=int(cos((fil_lat-1.0D0)*torad)/sint(j))*2

               DO NCY = 1,NN
                  DO I = 1,IMT
                     XS (I)= X (I,J)* Z (I,J,1)
                  END DO
               DO I = 2,IMM
                  X(I,J)=(XS(I)*(1.0D0-0.25D0*Z(I-1,J,1)-0.25D0*Z(I+1,J,1)) &
                            +0.25D0*(XS(I-1)+XS(I+1)))*Z(I,J,1)

               END DO
                  X (1,J)= X (IMM,J)
                  X (IMT,J)= X (2,J)
               END DO
!
            end if
         END DO


      RETURN
      END SUBROUTINE SMZ0


