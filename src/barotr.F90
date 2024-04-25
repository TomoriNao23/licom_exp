#define LOGMSG()
!call ZLOG(__LINE__ )
  SUBROUTINE ZLOG(line)
		use msg_mod
		use param_mod
		implicit none
		integer line
    call mpi_barrier(mpi_comm_ocn,ierr)
	end
!  CVS: $Id: barotr.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      SUBROUTINE BAROTR
!     =================
#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use dyn_mod
use work_mod
#ifdef SPMD
use msg_mod
#endif
      IMPLICIT NONE

      INTEGER :: IEB,NC,IEB_LOOP
      real(r8)    :: gstar ,am_viv,fil_lat1,fil_lat2
!      real(r8)    :: CDBAR(IMT,JMT)
!
!---------------------------------------------------------------------
!      Define the threthold latitute for zonal smoother
       fil_lat1=69.0D0
       fil_lat2=69.0D0
!---------------------------------------------------------------------
!     INITIALIZE WORK ARRAYS
!---------------------------------------------------------------------
      wka=0
      work=0
!---------------------------------------------------------------------
!     EULER BACKWARD SCHEME IS USED FOR THE FIRST STEP OF EVERY MONTH
!     IEB=0: LEAP-FROG SCHEME; IEB=1: EULER BACKWARD SCHEME
!---------------------------------------------------------------------
      IEB = 0 ; IEB_LOOP=0

      IF (ISB == 0)  THEN
         IEB = 1 ; IEB_LOOP=1
      END IF
       baro_loop : DO NC = 1,NBB+IEB_LOOP

#if ( defined TIDE )
!lhl20100801
             call equ_tide_mom
!             call equ_tide
#ifdef SPMD
       call exchange_2d(elevat)
#endif
!lhl20100801
#endif

      if (IEB==1.or.ISB>1) then

!---------------------------------------------------------------------
!     COMPUTE THE "ARTIFICIAL" HORIZONTAL VISCOSITY
!---------------------------------------------------------------------

#if ( defined SMAG1)
!$OMP PARALLEL DO PRIVATE (J,I)
         DO J = JST,JET
            DO I = 1,IMT
               WKA (I,J,11)= UBP (I,J)
               WKA (I,J,12)= VBP (I,J)
            END DO
         END DO
         CALL SMAG2 (1)
#if (defined SMAG_FZ)
         DO J = JSM,JEM
            DO I = 2,IMM
               WKA (I,J,5)= VIV (I,J,1)* (0.5* OUX (J)* (WKA (I +1,J,7) &
                            - WKA (I -1,J,7)) &
               - R2E (J)* WKA (I,J +1,8) + R2F (J)* WKA (I,J -1,8))
               WKA (I,J,6)= VIV (I,J,1)* (0.5* OUX (J)* (WKA (I +1,J,9) &
                            - WKA (I -1,J,9)) &
               - R3E (J)* WKA (I,J +1,10) + R3F (J)* WKA (I,J -1,10)    &
                            + R4E (J)* WKA (I,J,7))
            END DO
         END DO

!-new
#else
!$OMP PARALLEL DO PRIVATE (J,I)
         DO J = JSM,JEM
            DO I = 2,IMM
               WKA (I,J,5)= VIV (I,J,1)* (0.5* OUX (J)* (WKA (I +1,J,7) &
                            - WKA (I -1,J,7)) &
               - R2E (J)* WKA (I,J +1,8) + R2F (J)* WKA (I,J -1,8))
               WKA (I,J,6)= VIV (I,J,1)* (0.5* OUX (J)* (WKA (I +1,J,9) &
                            - WKA (I -1,J,9)) &
               - R3E (J)* WKA (I,J +1,10) + R3F (J)* WKA (I,J -1,10)    &
                            + R4E (J)* WKA (I,J,7))
            END DO
         END DO
#endif
#else
#if (defined BIHAR)
!$OMP PARALLEL DO PRIVATE (J,I)
         DO J = JSM,JEM
            DO I = 2,IMM
               WKA (I,J,7)= AM3(I,J,1)*VIV(I,J,1)*(R1D(J)*(UBP(I,J+1)-UBP(I,J))&
                           -R1C (J)*(UBP(I,J)-UBP(I,J-1)) &
                           +SOUX (J)* (UBP(I+1,J)-2.0*UBP (I,J)+UBP(I-1,J)) &
                           +CV1(J)*UBP(I,J)+CV2(J)*(VBP(I+1,J)-VBP(I-1,J)))
               WKA (I,J,8)= AM3 (I,J,1)* VIV(I,J,1)*(R1D(J)*(VBP(I,J+1)-VBP (I,J)) &
                           -R1C(J)*(VBP(I,J)-VBP(I,J-1))&
                           +SOUX (J)* (VBP(I +1,J) -2.0* VBP (I,J)+VBP (I -1,J)) &
                           + CV1 (J)*VBP(I,J)-CV2(J)*(UBP(I+1,J)-UBP(I-1,J)))
            END DO
            wka(1,j,7)=wka(imm,j,7)
            wka(imt,j,7)=wka(2,j,7)
            wka(1,j,8)=wka(imm,j,8)
            wka(imt,j,8)=wka(2,j,8)
         END DO
#ifdef SPMD
         call exchange_2d(wka(1,1,7))
         call exchange_2d(wka(1,1,8))
#endif

!$OMP PARALLEL DO PRIVATE (J,I)
         DO J = JSM,JEM
            DO I = 2,IMM
               WKA (I,J,5)= VIV (I,J,1)* ( R1D (J)* (WKA (I,J +1,7)     &
                           -WKA(I,J,7))-R1C (J)*(WKA(I,J,7)-WKA(I,J-1,7))+  &
                           SOUX (J)* (WKA (I+1,J,7) -2.0*               &
               WKA (I,J,7) + WKA (I -1,J,7)) + CV1 (J)* WKA (I,J,7)     &
                            + CV2 (J)*(WKA (I+1,J,8)-WKA (I-1,J,8)))
               WKA (I,J,6)= VIV (I,J,1)* ( R1D (J)* (WKA (I,J +1,8)     &
                            - WKA (I,J,8)) &
               - R1C (J)* (WKA (I,J,8) - WKA (I,J -1,8)) + SOUX (J)* (  &
                           WKA (I +1,J,8) -2.0*&
               WKA (I,J,8) + WKA (I -1,J,8)) + CV1 (J)* WKA (I,J,8)     &
                            - CV2 (J)* (WKA (I +1,J,7)-WKA (I-1,J,7)))
            END DO
            wka(1,j,5)=wka(imm,j,5)
            wka(imt,j,5)=wka(2,j,5)
            wka(1,j,6)=wka(imm,j,6)
            wka(imt,j,6)=wka(2,j,6)
         END DO
!!!!!!
#else
!$OMP PARALLEL DO PRIVATE (J,I,am_viv)
         DO J = JSM,JEM
            DO I = 2,IMM
               am_viv=AM3 (I,J,1)* VIV (I,J,1)
               WKA (I,J,5)=am_viv*(R1D(J)*(UBP(I,J+1)-UBP(I,J))-R1C(J)*  &
                           (UBP(I,J)-UBP(I,J-1))+SOUX(J)*(UBP(I+1,J)-    &
                           2.0*UBP(I,J)+UBP(I-1,J))+CV1(J)*UBP(I,J)+     &
                           CV2(J)*(VBP(I+1,J)-VBP(I-1,J)))
               WKA (I,J,6)=am_viv*(R1D(J)*(VBP(I,J+1)-VBP(I,J))-R1C(J)*  &
                           (VBP(I,J)-VBP(I,J-1))+SOUX(J)*(VBP(I+1,J)-    &
                           2.0* VBP (I,J)+VBP(I-1,J))+CV1(J)*VBP (I,J)-  &
                           CV2(J)*(UBP(I+1,J)-UBP(I-1,J)))
            END DO
         END DO
#endif
#endif

            IF (mod(isb,36)  == 1 ) THEN
!$OMP PARALLEL DO PRIVATE (J,I)
            DO J = JSM,JEM
               DO I = 2,IMM
                  DLUB (I,J)= DLUB (I,J) + WKA (I,J,5)
                  DLVB (I,J)= DLVB (I,J) + WKA (I,J,6)
               END DO
            END DO
            END IF
         END IF

!---------------------------------------------------------------------
!     + (g'-1)g*dH/dr
!---------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE (J,I,gstar)
         DO J = JSM,JEM
            DO I = 2,IMM
               gstar=(WGP (I,J) -1.0)*G *0.5
!lhl20100801
#if ( defined TIDE )
               WKA (I,J,1) = WKA (I,J,5) &
               + gstar*OUX(J)*(H0 (I,J) - H0 (I -1,J) + H0 (I,J +1) - H0 (I -1,J +1))*1.052 &
!               + gstar*OUX(J)*(H0 (I,J) - H0 (I -1,J) + H0 (I,J +1) - H0 (I -1,J +1))*0.948 &
!               + gstar*OUX(J)*(H0 (I,J) - H0 (I -1,J) + H0 (I,J +1) - H0 (I -1,J +1))*0.052 &
               - gstar*OUX(J)*(elevat (I,J) - elevat (I -1,J) + elevat (I,J +1) - elevat (I -1,J +1))
               WKA (I,J,2) = WKA (I,J,6) &
               + gstar*OUY(J)*(H0 (I,J +  1) - H0 (I,J) + H0 (I -1,J +1) - H0 (I -1,J))*1.052 &
!               + gstar*OUY(J)*(H0 (I,J +  1) - H0 (I,J) + H0 (I -1,J +1) - H0 (I -1,J))*0.948 &
!               + gstar*OUY(J)*(H0 (I,J +  1) - H0 (I,J) + H0 (I -1,J +1) - H0 (I -1,J))*0.052 &
               - gstar*OUY(J)*(elevat (I,J + 1) - elevat (I,J) + elevat (I -1,J +1) - elevat (I -1,J))
#else
               WKA (I,J,1) = WKA (I,J,5) &
               + gstar*OUX(J)*(H0 (I,J) - H0 (I -1,J) + H0 (I,J +1) - H0 (I -1,J +1)) 
               WKA (I,J,2) = WKA (I,J,6) &
               + gstar*OUY(J)*(H0 (I,J +  1) - H0 (I,J) + H0 (I -1,J +1) - H0 (I -1,J)) 
#endif
!lhl20100801
            END DO
         END DO


!---------------------------------------------------------------------
!     COMPUTING H0 AT U/V POINTS
!---------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE (J,I)
         DO J = JSM,JEM
            DO I = 2,IMM
               WORK (I,J) = 0.25* (H0 (I,J) + H0 (I -1,J) + H0 (I,J +1) &
                            + H0 (I -1,J +1))
            END DO

!     SET CYCLIC CONDITIONS ON EASTERN AND WESTERN BOUNDARY
!     DO J=2,JMM
            WORK (1,J) = WORK (IMM,J)
            WORK (IMT,J) = WORK (2,J)
         END DO

#ifdef SPMD
       call exchange_2d(work)
#endif

!Yu
!---------------------------------------------------------------------
!     COMPUTING DU & DV
!---------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE (J,I)
         DO J = JSM,JEM
            DO I = 2,IMM
               WKA (I,J,1)= VIV (I,J,1)* ( WKA (I,J,1) + DLUB (I,J)     &
                              - FF (J)* VBP (I,J) + &
               PAX (I,J) + PXB (I,J) - WORK (I,J)* WHX (I,J) )
               WKA (I,J,2)= VIV (I,J,1)* ( WKA (I,J,2) + DLVB (I,J)     &
                              + FF (J)* UBP (I,J) + &
               PAY (I,J) + PYB (I,J) - WORK (I,J)* WHY (I,J) )
            END DO
         END DO

!---------------------------------------------------------------------
!     CORIOLIS ADJUSTMENT
!---------------------------------------------------------------------

         IF (ISB == 0) THEN
!$OMP PARALLEL DO PRIVATE (J,I)
            DO J = JSM,JEM
               DO I = 2,IMM
                  WKA (I,J,3)= EBEA (J)* WKA (I,J,1) - EBEB (J)* WKA (I,J,2)
                  WKA (I,J,4)= EBEA (J)* WKA (I,J,2) + EBEB (J)* WKA (I,J,1)
               END DO
            END DO
         ELSE
!$OMP PARALLEL DO PRIVATE (J,I)
            DO J = JSM,JEM
               DO I = 2,IMM
                  WKA (I,J,3)= EBLA (J)* WKA (I,J,1) - EBLB (J)* WKA (I,J,2)
                  WKA (I,J,4)= EBLA (J)* WKA (I,J,2) + EBLB (J)* WKA (I,J,1)
               END DO
            END DO
         END IF


!     SET CYCLIC CONDITIONS ON EASTERN AND WESTERN BOUNDARY
!$OMP PARALLEL DO PRIVATE (J)
         DO J = JSM,JEM
            WKA (1,J,3) = WKA (IMM,J,3)
            WKA (IMT,J,3) = WKA (2,J,3)
            WKA (1,J,4) = WKA (IMM,J,4)
            WKA (IMT,J,4) = WKA (2,J,4)
         END DO


!---------------------------------------------------------------------
!     COMPUTING DH0
!---------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE (J,I)
         DO J = JST,JET
            DO I = 1,IMT
               WKA (I,J,1)= UB (I,J)* (DZPH (I,J) + WORK (I,J))
               WKA (I,J,2)= VB (I,J)* (DZPH (I,J) + WORK (I,J))
            END DO
         END DO


!$OMP PARALLEL DO PRIVATE (J,I)
         DO J = JSM,JEM
            DO I = 2,IMM
               WORK (I,J)=VIT(I,J,1)*(-1)*(0.5*OTX(J)*(WKA(I+1,J,1)+WKA(I+1,J-1,1) &
                          -WKA(I,J,1)-WKA(I,J-1,1))+2.0*R2A(J)*(WKA(I,J,2) + &
                          WKA(I+1,J,2))-2.0*R2B(J)*(WKA (I,J-1,2)+WKA (I+1,J-1,2)))
            END DO
!     ENDDO

!     SET CYCLIC CONDITIONS ON EASTERN AND WESTERN BOUNDARY

            WORK (1,J) = WORK (IMM,J)
            WORK (IMT,J) = WORK (2,J)
         END DO


!---------------------------------------------------------------------
!     PREDICTING VB , UB & H0
!---------------------------------------------------------------------
!#if ( defined TIDE )
!     DO J = 1,JMT
!        DO I = 1,IMT
!       if (abs(WORK(I,J)).gt.1.e-5) then
!       CDBAR(I,J)=0.95*1e-4
!       else if  ((abs(WORK(I,J)).gt.1.e-7).and.(abs(WORK(I,J)).le.1e-5)) then
!       CDBAR(I,J)=0.95*10*abs(WORK(I,J))
!       else
!       CDBAR(I,J)=0.95*1e-6
!       endif
!       WKA(I,J,3)=WKA(I,J,3)-CDBAR(I,J)*(UBP(I,J)**2+VBP(I,J)**2)**0.5*UBP(I,J)
!       WKA(I,J,4)=WKA(I,J,4)-CDBAR(I,J)*(UBP(I,J)**2+VBP(I,J)**2)**0.5*VBP(I,J)
!        END DO
!     ENDDO
!#endif
!YU  Oct. 24,2005
         CALL SMUV (WKA(1,1,3) ,VIV,1,fil_lat1)
         CALL SMUV (WKA(1,1,4) ,VIV,1,fil_lat1)
         CALL SMZ0 (WORK,VIT,fil_lat1)
!YU  Oct. 24,2005
!
         IF (ISB < 1) THEN

!
!$OMP PARALLEL DO PRIVATE (J,I)
         DO J = JSM,JEM
            DO I = 1,IMT
               UB (I,J)= UBP (I,J) + WKA (I,J,3)* DTB
               VB (I,J)= VBP (I,J) + WKA (I,J,4)* DTB
               H0 (I,J)= H0P (I,J) + WORK (I,J) * DTB
            END DO
         END DO

!---------------------------------------------------------------------
!     FILTER FORCING AT HIGT LATITUDES
!---------------------------------------------------------------------

            CALL SMUV (UB ,VIV,1,fil_lat1)
            CALL SMUV (VB ,VIV,1,fil_lat1)
            CALL SMZ0 (H0,VIT,fil_lat1)

#ifdef SPMD
         call exchange_pack(ub,vb,h0)
#endif

         IF (IEB == 0) THEN
            ISB = ISB +1
!$OMP PARALLEL DO PRIVATE (J,I)
!Yu         DO J = JSM,JEM
            DO J = JST,JET ! Dec. 4, 2002, Yongqiang YU
               DO I = 1,IMT
                  H0F (I,J) = H0F (I,J) + H0 (I,J)
                  H0BF (I,J) = H0BF (I,J) + H0 (I,J)
               END DO
            END DO
            cycle baro_loop
         END IF

         IEB = 0

         cycle baro_loop

      ELSE


!$OMP PARALLEL DO PRIVATE (J,I)
         DO J = JSM,JEM
            DO I = 1,IMT
               WKA (I,J,1) = UBP (I,J) + WKA (I,J,3)* DTB2
               WKA (I,J,2) = VBP (I,J) + WKA (I,J,4)* DTB2
               WORK(I,J)   = H0P (I,J) + WORK (I,J)* DTB2
            END DO
         END DO

!---------------------------------------------------------------------
!     FILTER FORCING AT HIGT LATITUDES
!---------------------------------------------------------------------


#ifdef SPMD
         call exchange_pack(wka(1,1,1),wka(1,1,2),work)
#endif
!$OMP PARALLEL DO PRIVATE (J,I)
         DO J = JST,JET
            DO I = 1,IMT
               UBP (I,J) = AFB2* UB (I,J) + AFB1* (UBP (I,J) + WKA (I,J,1))
               UB (I,J) = WKA (I,J,1)*VIV(I,J,1)
               VBP (I,J) = AFB2* VB (I,J) + AFB1* (VBP (I,J) + WKA (I,J,2))
               VB (I,J) = WKA (I,J,2)*VIV(I,J,1)
               H0P (I,J) = AFB2* H0 (I,J) + AFB1* (H0P (I,J) + WORK(I,J))
               H0 (I,J) = WORK (I,J)
            END DO
         END DO

!YU  Oct. 24,2005
!lhl0711         IF (MOD(ISB,1200)==0) THEN
         IF (MOD(ISB,720)==1) THEN
            CALL SMUV (UB ,VIV,1,fil_lat2)
            CALL SMUV (VB ,VIV,1,fil_lat2)
            CALL SMZ0 (H0 ,VIT,fil_lat2)
            CALL SMUV (UBP,VIV,1,fil_lat2)
            CALL SMUV (VBP,VIV,1,fil_lat2)
            CALL SMZ0 (H0P,VIT,fil_lat2)
         END IF
!        IF ( isb == 1 .and. mod(month,2) ==0 ) THEN
!           CALL SMOOTH (UB ,VIV(1,1,1),0.5,0)
!           CALL SMOOTH (VB ,VIV(1,1,1),0.5,0)
!           CALL SMOOTH (UBP,VIV(1,1,1),0.5,0)
!           CALL SMOOTH (VBP,VIV(1,1,1),0.5,0)
!           CALL SMOOTH (H0 ,VIT(1,1,1),0.5,1)
!           CALL SMOOTH (H0P,VIT(1,1,1),0.5,1)
!
!           call exchange_2d(ub)
!           call exchange_2d(vb)
!           call exchange_2d(ubp)
!           call exchange_2d(vbp)
!           call exchange_2d(h0)
!           call exchange_2d(h0p)
!        END IF

!YU  Oct. 24,2005

         ISB = ISB +1
      END IF



!$OMP PARALLEL DO PRIVATE (J,I)
         DO J = JST,JET ! Dec. 4, 2002, Yongqiang YU
            DO I = 1,IMT
               H0F (I,J) = H0F (I,J) + H0 (I,J)
               H0BF (I,J) = H0BF (I,J) + H0 (I,J)
            END DO
         END DO

#ifdef TIDE
!lhl20100801
      IF (ISB >= 1) time_tidal=time_tidal+1.0
!      IF (ISB >= 1) time_tidal=time_tidal+IDTB
!      if (mytid.eq.0) print*,time_tidal
!lhl20100801
#endif

      END DO baro_loop

      deallocate(dlub,dlvb)
      RETURN
      END SUBROUTINE BAROTR


