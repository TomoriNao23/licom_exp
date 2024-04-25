!  CVS: $Id: rdriver.F90,v 1.7 2003/08/25 07:47:52 lhl Exp $
!     ======================
      SUBROUTINE RDRIVER
!     ======================
 
#include <def-undef.h>
use param_mod
use pconst_mod
use forc_mod
use work_mod
#ifdef SPMD
use msg_mod
#endif

#ifdef COUP
use shr_sys_mod
#endif 

      IMPLICIT NONE
#include <netcdf.inc>
!
!     Define Variables.
      integer*4   :: ncid, iret
      integer*4,  dimension(4) :: start(4)
      integer*4,  dimension(4) :: count(4)
      integer*4 ::K1
 
!      REAL    :: WCOE (JMT),ABC
 
      if (mytid==0)then
      write(6,*)"Beginning------RDRIVER! "
#ifdef COUP
      call shr_sys_flush(6)
#endif
      endif 

!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = 1,JMT
            DO I = 1,IMT
               WKA (I,J,K)= 0.0
            END DO
         END DO
      END DO
 
 
!-----------------------------------------------------------------------
!     SU3  : Sea surface zonal wind stress         (N/M**2)
!     SV3  : Sea surface meridional wind stress    (N/M**2)
!     PSA3 : Sea surface air pressure              (Pa)
!     SWV3 : Total net downward solar radiation    (W/M**2)
!    NSWV3 : None Solar flux                       (Wm-2)
!    DQDT3 : Dq/Dt                                 (WK-1m-2)
!     SST3 : Sea surface temperature               (Celsius)
!     SSS3 : Sea surface salinity                  (psu)
!   chloro3:chlorophll concentration               (mg m-3)
!-----------------------------------------------------------------------
 
!     READ FORCING FIELD
#ifdef SPMD
      if(mytid==0)then
!
#if (!defined CDFIN)
      OPEN (90,FILE ='MODEL.FRC',STATUS ='OLD',FORM ='UNFORMATTED')
      READ (90) SWV3_IO,NSWV3_IO,DQDT3_IO,SU3_IO,SV3_IO,SST3_IO,SSS3_IO
      CLOSE (90)
#if (defined SOLARCHLORO)
      OPEN (91,FILE ='MODEL_CHLFRC',STATUS ='OLD',FORM ='UNFORMATTED')
      READ (91) chloro3_io
      CLOSE (91)
#endif

#if (defined WAVE_FORCE)
      OPEN (91,FILE ='MODEL_WAVEFRC',STATUS ='OLD',FORM ='UNFORMATTED')
      READ (91) tbv3_io
      CLOSE (91)
#endif

#else
!
!----------------------------------------------------
! Open netCDF file.
!----------------------------------------------------
      iret=nf_open('MODEL.FRC',nf_nowrite,ncid)
      call check_err (iret)
!----------------------------------------------------
!   Retrieve data
!----------------------------------------------------
      start(1)=1 ; count(1)=imt
      start(2)=1 ; count(2)=jmt_global
      start(3)=1 ; count(3)=1
      start(4)=1 ; count(4)=12

      iret=nf_get_vara_double(ncid,   5,start,count,swv3_io)
      call check_err (iret)
      iret=nf_get_vara_double(ncid,   6,start,count,nswv3_io)
      call check_err (iret)
      iret=nf_get_vara_double(ncid,   7,start,count,dqdt3_io)
      call check_err (iret)
      iret=nf_get_vara_double(ncid,   8,start,count,su3_io)
      call check_err (iret)
      iret=nf_get_vara_double(ncid,   9,start,count,sv3_io)
      call check_err (iret)
      iret=nf_get_vara_double(ncid,  10,start,count,sst3_io)
      call check_err (iret)
      iret=nf_get_vara_double(ncid,  11,start,count,sss3_io)
      call check_err (iret)
!
      iret = nf_close (ncid)
      call check_err (iret)
!===============================================
!input the chlorophyll concentration

#if (defined SOLARCHLORO)
      iret=nf_open('MODEL_CHLFRC',nf_nowrite,ncid)
      call check_err (iret)
      iret=nf_get_vara_double(ncid,      5,start,count,chloro3_io)
!     swv3_io must be change
      call check_err (iret)
      iret = nf_close (ncid)
      call check_err (iret)
#endif
!==============================================
!      write(*,*) seaice3_io(180,:,12)

!==============================================
!input annual mean runoff
!
!===============================================
!input the wave force
#if (defined WAVEFORCE)

      start(1)=1 ; count(1)=imt
      start(2)=1 ; count(2)=jmt_global
      start(3)=1 ; count(3)=km
      start(4)=1 ; count(4)=12

      iret=nf_open('MODEL_WAVEFRC',nf_nowrite,ncid)
      call check_err (iret)
      iret=nf_get_vara_double(ncid,      5,start,count,tbv3_io)
      call check_err (iret)
      iret = nf_close (ncid)
      call check_err (iret)
      
        DO K = 1,12
         DO J = 1,jmt_global
          DO I = 1,IMT
           do K1=1,KM
            tbv3_io(i,j,k1,k)=tbv3_io(i,j,k1,k)*vit_global(i,j,k1)
           end do
           end do
           end do
           end do
      
#endif
!==============================================

#endif
      end if
!
      call global_to_local_4d(su3_io,su3,12,1)
      call global_to_local_4d(sv3_io,sv3,12,1)
      call global_to_local_4d(sss3_io,sss3,12,1)
      call global_to_local_4d(swv3_io,swv3,12,1)
      call global_to_local_4d(nswv3_io,nswv3,12,1)
      call global_to_local_4d(dqdt3_io,dqdt3,12,1)
      call global_to_local_4d(sst3_io,sst3,12,1)
      call global_to_local_4d(chloro3_io,chloro3,12,1)
      call global_to_local_4d(seaice3_io,seaice3,12,1)
      call global_to_local_4d(runoff3_io,runoff3,1,1)
      call global_to_local_4d(tbv3_io,tbv3,12,km)
!
#else
!Yu
#if (!defined CDFIN)
      OPEN (90,FILE ='MODEL.FRC',STATUS ='OLD',FORM ='UNFORMATTED')
      READ (90) SWV3_io,NSWV3_io,DQDT3_io,SU3_io,SV3_io,SST3_io,SSS3_io
      CLOSE (90)
      OPEN (91,FILE ='MODEL_CHLFRC',STATUS ='OLD',FORM ='UNFORMATTED')
      READ (91) chloro3_io
      CLOSE (91)
      OPEN (91,FILE ='MODEL_WAVEFRC',STATUS ='OLD',FORM ='UNFORMATTED')
      READ (91) tbv3_io
      CLOSE (91)
!
#else
!----------------------------------------------------
! Open netCDF file.
!----------------------------------------------------
      iret=nf_open('MODEL.FRC',nf_nowrite,ncid)
      call check_err (iret)

!----------------------------------------------------
!   Retrieve data
!----------------------------------------------------
      start(1)=1 ; count(1)=imt
      start(2)=1 ; count(2)=jmt_global
      start(3)=1 ; count(3)=1
      start(4)=1 ; count(4)=12

      iret=nf_get_vara_double(ncid,   5,start,count,swv3_io)
      call check_err (iret)
      iret=nf_get_vara_double(ncid,   6,start,count,nswv3_io)
      call check_err (iret)
      iret=nf_get_vara_double(ncid,   7,start,count,dqdt3_io)
      call check_err (iret)
      iret=nf_get_vara_double(ncid,   8,start,count,su3_io)
      call check_err (iret)
      iret=nf_get_vara_double(ncid,   9,start,count,sv3_io)
      call check_err (iret)
      iret=nf_get_vara_double(ncid,  10,start,count,sst3_io)
      call check_err (iret)
      iret=nf_get_vara_double(ncid,  11,start,count,sss3_io)
      call check_err (iret)
!
      iret = nf_close (ncid)
      call check_err (iret)
#if (defined SOLARCHLORO)
      iret=nf_open('MODEL_CHLFRC',nf_nowrite,ncid)
      call check_err (iret)
      iret=nf_get_vara_double(ncid,      5,start,count,chloro3_io)
      call check_err (iret)
      iret = nf_close (ncid)
      call check_err (iret)
#endif

#if (defined WAVEFORCE)

      start(1)=1 ; count(1)=imt
      start(2)=1 ; count(2)=jmt_global
      start(3)=1 ; count(3)=km
      start(4)=1 ; count(4)=12
      
      iret=nf_open('MODEL_WAVEFRC',nf_nowrite,ncid)
      call check_err (iret)
      iret=nf_get_vara_double(ncid,      5,start,count,tbv3_io)
      call check_err (iret)
      iret = nf_close (ncid)
      call check_err (iret)
#endif
!
#endif

      SWV3=SWV3_io
      NSWV3=NSWV3_io
      DQDT3=DQDT3_io
      SU3=SU3_io
      SV3=SV3_io
      SST3=SST3_io
      SSS3=SSS3_io
      chloro3=chloro3_io
      tbv3=tbv3_io
#endif
!Yu
      if (mytid==0) then
#ifdef COUP
      call shr_sys_flush(6)
#endif
      end if
 
!-----------------------------------------------------------------------
!     land grids of the forcing fields assigned to 0 
!-----------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,12
         DO J = 1,JMT
            DO I = 1,IMT
                SWV3(I,J,K)= SWV3(I,J,K)*VIT(I,J,1)
               NSWV3(I,J,K)=NSWV3(I,J,K)*VIT(I,J,1)
               DQDT3(I,J,K)=DQDT3(I,J,K)*VIT(I,J,1)
                 SU3(I,J,K)=  SU3(I,J,K)*VIV(I,J,1)
                 SV3(I,J,K)=  SV3(I,J,K)*VIV(I,J,1)
                SST3(I,J,K)= SST3(I,J,K)*VIT(I,J,1)
                SSS3(I,J,K)= SSS3(I,J,K)*VIT(I,J,1)
                SEAICE3(I,J,K)= SEAICE3(I,J,K)*VIT(I,J,1)
#if (defined SOLARCHLORO)
	       chloro3(I,J,K)= chloro3(I,J,K)*VIT(I,J,1) 
#endif
            
#if (defined WAVEFORCE)
            DO K1 = 1,KM
	       tbv3(I,J,K1,K)= tbv3(I,J,K1,K)*VIT(I,J,K1)
	     END DO
#endif
              END DO
           END DO
        END DO

!-----------------------------------------------------------------------
!     salinity = (psu-35)*0.001
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,12
         DO J = 1,JMT
            DO I = 1,IMT
               SSS3 (I,J,K) = (SSS3 (I,J,K) -35.0D0)*0.001D0
            END DO
         END DO
      END DO
 
!-----------------------------------------------------------------------
!     reverse VS (southward is positive)
!     notice: the former program VS is reversed during preparing forcing field
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,12
         DO J = 1,JMT
            DO I = 1,IMT
               SV3 (I,J,K) = -SV3 (I,J,K)
            END DO
         END DO
      END DO
 
 
!-----------------------------------------------------------------------
!     CALCULATING THE ANNUAL MEAN FORCING FIELD 
!-----------------------------------------------------------------------
 
#if (defined FRC_ANN)
 
      DO M = 1,12
!$OMP PARALLEL DO PRIVATE (J,I)
         DO J = 1,JMT
            DO I = 1,IMT
               WKA (I,J,1)= WKA (I,J,1) + SU3 (I,J,M)
               WKA (I,J,2)= WKA (I,J,2) + SV3 (I,J,M)
               WKA (I,J,3)= WKA (I,J,3) + SSS3 (I,J,M)
               WKA (I,J,4)= WKA (I,J,4) + SWV3 (I,J,M)
               WKA (I,J,5)= WKA (I,J,5) + SST3 (I,J,M)
               WKA (I,J,6)= WKA (I,J,6) + NSWV3 (I,J,M)
               WKA (I,J,7)= WKA (I,J,7) + DQDT3 (I,J,M)
               WKA (I,J,8)= WKA (I,J,8) + seaice3 (I,J,M)
#if (defined SOLARCHLORO)
               WKA (I,J,9)= WKA (I,J,9) + chloro3 (I,J,M)
#endif
            END DO
         END DO
         
!     DO J = 1,JMT
!       DO I = 1,IMT
!       DO k=1,km    
!#if (defined WAVEFORCE)
!               WKA (I,J,8)= WKA (I,J,8) + tbv3 (I,J,M)
!#endif
!        END DO
!        END DO  
!        END DO
      END DO  
 
!$OMP PARALLEL DO PRIVATE (M,J,I)
      DO M = 1,12
         DO J = 1,JMT
            DO I = 1,IMT
               SU3 (I,J,M) = WKA (I,J,1)/12.0D0
               SV3 (I,J,M) = WKA (I,J,2)/12.0D0
              SSS3 (I,J,M) = WKA (I,J,3)/12.0D0
              SWV3 (I,J,M) = WKA (I,J,4)/12.0D0
              SST3 (I,J,M) = WKA (I,J,5)/12.0D0
             NSWV3 (I,J,M) = WKA (I,J,6)/12.0D0
             DQDT3 (I,J,M) = WKA (I,J,7)/12.0D0
             seaice3 (I,J,M) = WKA (I,J,8)/12.0D0
#if (defined SOLARCHLORO)
             chloro3 (I,J,M) = WKA (I,J,9)/12.0D0
#endif

            END DO
         END DO
      END DO

#endif
!
      if (mytid==0)then
      write(6,*)"END-----------RDRIVER!"
#ifdef COUP
      call shr_sys_flush(6)
#endif
      endif 
 
      RETURN
      END SUBROUTINE RDRIVER
 
