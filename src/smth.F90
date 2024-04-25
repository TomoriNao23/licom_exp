!  CVS: $Id: smth.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =========================
      SUBROUTINE SMOOTH (X,VVV,S,itv)
!     =========================
!   9-point smooth
use precision_mod
use param_mod
      IMPLICIT NONE
      integer :: itv
      REAL(r8)    :: X (imt,jmt),XS (imt,jmt),VVV(imt,jmt)
      real(r8) :: IT5,IT6,IT7,IT8,IT1,IT2,IT3,IT4
      REAL(r8)    :: XS5,XS6,XS7,XS8,XS1,XS2,XS3,XS4
      REAL(r8)    :: S1,BT0,BT1,BT2,S
      S1 = 1.0D0- S
 
      BT0 = S1* S1
      BT1 = 0.5D0* S * S1
      BT2 = 0.25D0* S * S
      DO J = 1,jmt
         DO I = 1,imt
            XS (I,J)= X (I,J)*vvv(i,j)
         END DO
 
      END DO
 
!
    if ( itv == 1) then 
      JJJ : DO J = 2,JMM
         III : DO I = 2,IMM
            IF (vvv (I,J)  < 0.5 ) CYCLE III
            IT5 = vvv (I +1,J +1)
            XS5 = XS (I +1,J +1)
            IT6 = vvv (I -1,J +1)
            XS6 = XS (I -1,J +1)
            IT7 = vvv (I -1,J -1)
            XS7 = XS (I -1,J -1)
            IT8 = vvv (I +1,J -1)
            XS8 = XS (I +1,J -1)
            IT1 = vvv (I +1,J)
            XS1 = XS (I +1,J)
            IT2 = vvv (I -1,J)
            XS2 = XS (I -1,J)
            IT3 = vvv (I,J +1)
            XS3 = XS (I,J +1)
            IT4 = vvv (I,J -1)
            XS4 = XS (I,J -1)
            IF (IT1 < 0.5 ) XS1 = XS (I,J)
            IF (IT2 < 0.5) XS2 = XS (I,J)
            IF (IT3 < 0.5) XS3 = XS (I,J)
            IF (IT4 < 0.5) XS4 = XS (I,J)
            IF (IT5 < 0.5) XS5 = XS (I,J)
            IF (IT6 < 0.5) XS6 = XS (I,J)
            IF (IT7 < 0.5) XS7 = XS (I,J)
            IF (IT8 < 0.5) XS8 = XS (I,J)
            X (I,J)= BT0* XS (I,J) + BT1* (XS1+ XS2+ XS3+ XS4) + BT2* ( &
                 XS5+ XS6+ XS7+ XS8)
         END DO III
      END DO JJJ
      DO J = 2,JMM
         X (1,J) = X (IMM,J)
         X (IMT,J) = X (2,J)
      END DO
   else
      KKK : DO J = 2,JMM
         LLL : DO I = 2,IMM
            IF (vvv (I,J)  < 0.5 ) CYCLE LLL
            XS5 = XS (I +1,J +1)
            XS6 = XS (I -1,J +1)
            XS7 = XS (I -1,J -1)
            XS8 = XS (I +1,J -1)
            XS1 = XS (I +1,J)
            XS2 = XS (I -1,J)
            XS3 = XS (I,J +1)
            XS4 = XS (I,J -1)
            X (I,J)= BT0* XS (I,J) + BT1* (XS1+ XS2+ XS3+ XS4) + BT2* ( &
                 XS5+ XS6+ XS7+ XS8)
         END DO LLL
      END DO KKK
      DO J = 2,JMM
         X (1,J) = X (IMM,J)
         X (IMT,J) = X (2,J)
      END DO

   end if
 
      RETURN
      END SUBROUTINE SMOOTH
 

 
