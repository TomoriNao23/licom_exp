!  CVS: $Id: intfor.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      SUBROUTINE INTFOR
!     =================
!     INTERPOLATES MONTH MEAN FILEDS
 
#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use forc_mod
#if ( defined SPMD )
use msg_mod
#endif
      IMPLICIT NONE
 
      INTEGER :: IPT1,IPT2
      REAL(r8):: FACTOR
      INTEGER :: IYEAR,IREC
      CHARACTER (LEN=180) :: FNAME 
      CHARACTER (LEN=4)   :: FHEAD 
      CHARACTER (LEN=3)   :: FTAIL(12)
      REAL(r8),dimension(:,:),allocatable :: sx,sy
!      REAL(r4),dimension(imt,jmt) :: sx_in,sy_in 
!      REAL(r8),dimension(imt,jmt) :: sx,sy 
!      REAL(r8) :: EK
      REAL(r8) :: epsln
      allocate(sx(imt,jmt))
      allocate(sy(imt,jmt))
!
      epsln = 1.0D-25
!
      FTAIL=RESHAPE((/'jan','feb','mar','apr','may','jun', &
                   'jul','aug','sep','oct','nov','dec'/),(/12/)) 
!
! 
      IF ( IDAY <= 15) THEN
         IPT1 = MON0-1
         IF (IPT1 == 0 ) IPT1 = 12
         IPT2 = MON0
         FACTOR = FLOAT (IDAY -15)/ FLOAT (NMONTH (IPT1)) + 1
      ELSE
         IPT1 = MON0
         IPT2 = MOD (MON0,12) +1
         FACTOR = FLOAT (IDAY -15)/ FLOAT (NMONTH (IPT1))
      END IF
 
 
!lhl0711
            PSA3=1000.
!lhl0711

!$OMP PARALLEL DO PRIVATE (J,I)
      DO J = JST,JET
         DO I = 1,IMT
            SU (I,J)= ( (SU3 (I,J,IPT2) - SU3 (I,J,IPT1))               &
                    * FACTOR + SU3 (I,J,IPT1))
            SV (I,J)= ( (SV3 (I,J,IPT2) - SV3 (I,J,IPT1))               &
                    * FACTOR + SV3 (I,J,IPT1))
            PSA (I,J)= (PSA3 (I,J,IPT2) - PSA3 (I,J,IPT1))              &
                    * FACTOR + PSA3 (I,J, IPT1)
            TSA (I,J)= (TSA3 (I,J,IPT2) - TSA3 (I,J,IPT1))              &
                    * FACTOR + TSA3 (I,J, IPT1)
            SSS (I,J)= (SSS3 (I,J,IPT2) - SSS3 (I,J,IPT1))              &
                    * FACTOR + SSS3 (I,J, IPT1)
            SWV (I,J)= ( (SWV3 (I,J,IPT2) - SWV3 (I,J,IPT1))            &
                    * FACTOR + SWV3 (I,J, IPT1))
            UVA (I,J)= ( (UVA3 (I,J,IPT2) - UVA3 (I,J,IPT1))            &
                    * FACTOR + UVA3 (I,J, IPT1))
            QAR (I,J)= ( (QAR3 (I,J,IPT2) - QAR3 (I,J,IPT1))            &
                    * FACTOR + QAR3 (I,J, IPT1))
            CLD (I,J)= ( (CLD3 (I,J,IPT2) - CLD3 (I,J,IPT1))            &
                    * FACTOR + CLD3 (I,J, IPT1))
            SST (I,J)= ( (SST3 (I,J,IPT2) - SST3 (I,J,IPT1))            &
                    * FACTOR + SST3 (I,J, IPT1))
!lhl
           NSWV (I,J)= ( (NSWV3(I,J,IPT2) - NSWV3(I,J,IPT1))            &
                    * FACTOR + NSWV3 (I,J, IPT1))
           DQDT (I,J)= ( (DQDT3(I,J,IPT2) - DQDT3(I,J,IPT1))            &
                    * FACTOR + DQDT3(I,J, IPT1))
!lhl
         END DO
      END DO
!lhl1204
!
! Calculate the friction velocity at T-grid
!
      SX=0.0D0
      SY=0.0D0
         DO J = JSM,JET
            DO I = 2,IMM
               sx(i,j)= VIT(I,J,1)*(su(i,j  )+su(i+1,j  )&
                                   +su(i,j-1)+su(i+1,j-1)) &
                                 /(VIV(i,j  ,1)+VIV(i+1,j  ,1) &
                                  +VIV(i,j-1,1)+VIV(i+1,j-1,1) + epsln)
               sy(i,j)= VIT(I,J,1)*(sv(i,j  )+sv(i+1,j  )&
                                   +sv(i,j-1)+sv(i+1,j-1)) &
                                 /(VIV(i,j  ,1)+VIV(i+1,j  ,1) &
                                  +VIV(i,j-1,1)+VIV(i+1,j-1,1) + epsln)
            END DO
               sx(IMT,j)= sx(2,j)
               sy(IMT,j)= sy(2,j)
               sx(1,j)= sx(imm,j)
               sy(1,j)= sy(imm,j)
        END DO
!

      DO J = JST,JET
         DO I = 1,IMT
         USTAR(I,J)=sqrt(sqrt(sx(I,J)*sx(I,J)+sy(I,J)*sy(I,J))*OD0)*vit(i,j,1)
         END DO
      END DO
#ifdef DEBUG
      call chk_var2d(sx,"sx",1)
      call chk_var2d(sy,"sy",1)
      call chk_var2d(ustar,"us",1)
#endif
!
!!lhl1204

        DO J = JST,JET
         DO I = 1,IMT
           seaice(I,J)= ( (seaice3 (I,J,IPT2) - seaice3 (I,J,IPT1))     &
                    * FACTOR + seaice3 (I,J,IPT1))
           END DO
        END DO

        DO J = JST,JET
         DO I = 1,IMT
           runoff(I,J)= runoff3(i,j,1)
           END DO
        END DO

#if (defined SOLARCHLORO) 
        DO J = JST,JET
         DO I = 1,IMT
           chloro(I,J)= ( (chloro3 (I,J,IPT2) - chloro3 (I,J,IPT1))     &
                    * FACTOR + chloro3 (I,J,IPT1))
           END DO
        END DO
#endif

#if (defined WAVEFORCE)
         DO J = JST,JET
         DO I = 1,IMT
         DO K=1,KM
           tbv(I,J,K)= ( (tbv3 (I,J,K,IPT2) - tbv3 (I,J,K,IPT1))     &
                    * FACTOR + tbv3 (I,J,K,IPT1))
         END DO
         END DO
         END DO
#endif
 
 
 
 
#if ( defined FRC_DAILY)
#ifdef SPMD
      if (mytid==0)then
      IYEAR = 1078+ IYFM
!      PRINT *,IYEAR,MON0,IDAY
      write (FHEAD,'(i4)') IYEAR
      FNAME ="/export/home/lhl/LICOM2.0/forcing/ew"//FHEAD//FTAIL(MON0)//".dat"
!      PRINT *,FNAME
      OPEN (110,FILE = FNAME,FORM ='UNFORMATTED',ACCESS ='DIRECT',      &
                      RECL = JMT_global * IMT*4)
      IREC = IDAY
      READ (110,REC = IREC) ( (SU_in_io (I,J),I = 1,IMT),J = 1,JMT_global)
      CLOSE (110)

      FNAME ="/export/home/lhl/LICOM2.0/forcing/ns"//FHEAD//FTAIL(MON0)//".dat"
      OPEN (110,FILE = FNAME,FORM ='UNFORMATTED',ACCESS ='DIRECT',      &
                      RECL = JMT_global * IMT*4 )
      IREC = IDAY
      READ (110,REC = IREC) ( (SV_in_io (I,J),I = 1,IMT),J = 1,JMT_global)
      CLOSE (110)
       SU_io=SU_in_io
       SV_io=SV_in_io
!
!      DO J = 1,JMT_global
!         DO I = 1,IMT
!            SU_io (I,J)=  SU_io (I,J)*VIV_global(I,J,1)
!            SV_io (I,J)=  SV_io (I,J)*VIV_global(I,J,1)
!         END DO
!      END DO
!
      end if
!
!      if (mytid==0)then
!      OPEN (111,FILE = 'aa.dat',FORM ='UNFORMATTED',ACCESS ='DIRECT',      &
!                    RECL = JMT_global * IMT*4)
!      write (111,REC = 1) ( (SU_io (I,J),I = 1,IMT),J = 1,JMT_global)
!      write (111,REC = 2) ( (SV_io (I,J),I = 1,IMT),J = 1,JMT_global)
!      CLOSE (111)
!	stop 99999
!      end if
!
!Yu
!      call mpi_bcast(su_io,imt*jmt_global,mpi_real,0,mpi_comm_ocn,ierr)
!      call mpi_bcast(sv_io,imt*jmt_global,mpi_real,0,mpi_comm_ocn,ierr)
!      call mpi_barrier(mpi_comm_ocn,ierr)
!      call global_to_local_2d(su_io,sx)
!      call global_to_local_2d(sv_io,sy)
       call global_to_local_4d(su_io,sx,1,1)
       call global_to_local_4d(sv_io,sy,1,1)
!Yu
#else
      IYEAR = 1078+ IYFM
      PRINT *,IYEAR,MON0,IDAY
      write (FHEAD,'(i4)') IYEAR
      FNAME ="/export/home/lhl/LICOM2.0/forcing/ew"//FHEAD//FTAIL(MON0)//".dat"
      OPEN (110,FILE = FNAME,FORM ='UNFORMATTED',ACCESS ='DIRECT',      &
                      RECL = JMT * IMT*4)
      IREC = IDAY
      READ (110,REC = IREC) ( (SX_in (I,J),I = 1,IMT),J = 1,JMT)
      CLOSE (110)
 
      FNAME ="/export/home/lhl/LICOM2.0/forcing/ns"//FHEAD//FTAIL(MON0)//".dat"
      OPEN (110,FILE = FNAME,FORM ='UNFORMATTED',ACCESS ='DIRECT',      &
                      RECL = JMT * IMT*4)
      IREC = IDAY
      READ (110,REC = IREC) ( (SY_in (I,J),I = 1,IMT),J = 1,JMT)
      CLOSE (110)
       SX=SX_in
       SY=SY_in
#endif
 
      DO J = 1,JMT
         DO I = 1,IMT
            SU (I,J)=   SX (I,J)
            SV (I,J)= - SY (I,J)
         END DO
      END DO
 
#endif
 
      deallocate(sx)
      deallocate(sy)
      RETURN
      END SUBROUTINE INTFOR
 
 
