!  CVS: $Id: grids.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     ================
      SUBROUTINE GRIDS
!     ================
!     TOPOGRAPHY & GRIDS
!-----------------------------------------------------------------------
!
! Purpose: Set up some constants dependent on model grids.
!
! Author: Yongqiang Yu and Hailong Liu, Dec. 31, 2002
!
!
!-----------------------------------------------------------------------


#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use pmix_mod
use work_mod
use cdf_mod, only : start1,count1,start2,count2,start3,count3,start4,count4
#ifdef SPMD
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
#endif
      IMPLICIT NONE
#include <netcdf.inc>

!
!     Define Variables.
      integer*4   :: ncid, iret
!
      REAL(r8)    :: ZKP_IN (KMP1)
!      REAL(r8)    :: ZKP (KMP1)
      REAL(r8)    :: OMEGA,DX,ABCD,YU,CURU,EPS
      REAL(r8)    :: rpart,efold1,efold2,swarg1,swarg2
      REAL(r8)    :: AJQ,rscl1,rscl2
      REAL(r8)    :: ivk(imt,jmt_global)
      INTEGER :: IREC
      REAL(r4)    :: am3_io(imt,jmt_global,km)
      REAL(r8)    :: am3_tmp(imt,jmt_global,km)
      REAL(r8)    :: am3_u(imt,jmt_global,km)
      REAL(r8)    :: lat_r8(jmt_global)

#ifdef SPMD

      allocate(work_global(imt,jmt_global))
!YU   allocate(oux_global(jmt_global),ouy_global(jmt_global))
!YU   allocate(otx_global(jmt_global),ff_global(jmt_global))
!YU   allocate(sotx_global(jmt_global),soux_global(jmt_global))
!YU   allocate(cv1_global(jmt_global),cv2_global(jmt_global))
!YU   allocate(snlat_global(jmt_global))
!YU   allocate(sint_global(jmt_global),sinu_global(jmt_global))
!YU   allocate(dyr_global(jmt_global),dyt_global(jmt_global))
!YU   allocate(dxdyu_global(jmt_global),dxdyt_global(jmt_global))
!YU   allocate(r1a_global(jmt_global),r1b_global(jmt_global))
!YU   allocate(r2a_global(jmt_global),r2b_global(jmt_global))
!YU   allocate(r1c_global(jmt_global),r1d_global(jmt_global))
!YU   allocate(r2c_global(jmt_global),r2d_global(jmt_global))
!YU   allocate(ebea_global(jmt_global),ebeb_global(jmt_global))
!YU   allocate(ebla_global(jmt_global),eblb_global(jmt_global))
!YU   allocate(epea_global(jmt_global),epeb_global(jmt_global))
!YU   allocate(epla_global(jmt_global),eplb_global(jmt_global))

!YU   allocate(hbx_global(imt,jmt_global),hby_global(imt,jmt_global))
!YU   allocate(ohbu_global(imt,jmt_global),ohbt_global(imt,jmt_global))
!YU   allocate(dzph_global(imt,jmt_global))
!Yu
!YU   allocate(vit_global(imt,jmt_global,km),viv_global(imt,jmt_global,km))

#if  ( defined SMAG)
      allocate(cxt_global(jmt_global),cxu_global(jmt_global))
      allocate(cyt_global(jmt_global),cyu_global(jmt_global))
      allocate(r1e_global(jmt_global),r1f_global(jmt_global))
      allocate(r2e_global(jmt_global),r2f_global(jmt_global))
      allocate(r3e_global(jmt_global),r3f_global(jmt_global))
      allocate(r4e_global(jmt_global),r4f_global(jmt_global))
#endif
      allocate(cost_global(jmt_global),cosu_global(jmt_global))
#endif
      if (mytid==0)then
      write(6,*)"Beginning------GRIDS !"
      endif

!XC
#if (defined TSPAS)
      if (mytid==0)then
      write(6,*)"Use TSPAS advection scheme!"
      endif
#else
      if (mytid==0)then
      write(6,*)"Use CTCS advection scheme!"
      endif
#endif
!XC

!--------------------------------------------------------------
!     SET LOCAL CONSTANTS
!--------------------------------------------------------------
!     AM_TRO=2.0E+3
!     AM_EXT=2.0E+5
!     DLAM  =0.5
      PI = 4.0* ATAN (1.0)
      TORAD = PI /180.0D0
      RADIUS = 6371000.0D0
      OMEGA = 0.7292D-4

!--------------------------------------------------------------
!     ZONAL RESOULTION
!--------------------------------------------------------------
      DX = RADIUS * DLAM * TORAD
!lhl
#if (defined NETCDF) || (defined ALL)
      DO i = 1,imt
         lon (i)= 0.0+ (i -1)*DLAM
      END DO

#endif
!lhl
#ifdef TIDE
!lhl20100801
      lon_tidal=lon
!lhl20100801
#endif

!--------------------------------------------------------------
!     READ SCRIPPS T63 TOPOGRAPHY AND PRODUCE V-TYPE FIELDS
!--------------------------------------------------------------
#ifdef SPMD
      if (mytid==0) then
#if (!defined CDFIN)
      OPEN (71,FILE ='INDEX.DATA',FORM ='UNFORMATTED',STATUS ='OLD')
       READ (71)VIT_IN_GLOBAL,DYR_IN_GLOBAL,ZKP_IN
      CLOSE (71)
#else
!----------------------------------------------------
! Open netCDF file.
!----------------------------------------------------
      iret=nf_open('INDEX.DATA',nf_nowrite,ncid)
      call check_err (iret)
!----------------------------------------------------
!   Retrieve data
!----------------------------------------------------
      start3(1)=1 ; count3(1)=imt
      start3(2)=1 ; count3(2)=jmt_global
      start3(3)=1 ; count3(3)=km
      iret=nf_get_vara_double(ncid,   5,start3,count3,vit_in_global)
      call check_err (iret)
!
      start1(1)=1 ; count1(1)=jmt_global
      iret=nf_get_vara_double(ncid,   6,start1,count1,dyr_in_global)
      call check_err (iret)
      
!
      start1(1)=1 ; count1(1)=kmp1
      iret=nf_get_vara_double(ncid,   7,start1,count1,zkp_in)
      call check_err (iret)
!
      iret = nf_close (ncid)
      call check_err (iret)
!YU
!-------------------------------------------------------
!    Read background vertical diffusion coefficient
!-------------------------------------------------------
#ifdef COUP
#ifdef SPMD
      ahv_back =0.0D0
      open (27,file="ahv_back.txt",form="formatted")
      read(27,'(5x,10D10.2)') ahv_back_global
      close(27)

#endif
#endif
!YU
!-------------------------------------------------------
!    Read Basin Index Field
!-------------------------------------------------------
   !   iret=nf_open('BASIN.nc',nf_nowrite,ncid)
   !   call check_err (iret)
!
    !  start2(1)=1 ; count2(1)=imt
    !  start2(2)=1 ; count2(2)=jmt_global
    !  iret=nf_get_vara_int(ncid,   3,start2,count2,basin)
    !  call check_err (iret)
    !  iret = nf_close (ncid)
    !  call check_err (iret)
!
#endif
       VIT_GLOBAL=VIT_IN_GLOBAL
       DYR_GLOBAL=DYR_IN_GLOBAL
       ZKP=ZKP_IN
!
!Yu
      do k=1,km
      do j=1,jst_global
      do i=1,imt
         vit_global(i,j,k)=0.0
      end do
      end do
      end do
!!lhl100913
!      do k=1,km
!      do j=123,126
!      do i=143,144
!         vit_global(i,j,k)=0.0
!      end do
!      end do
!      end do
!      do k=1,km
!      do j=130,135
!      do i=138,142
!         vit_global(i,j,k)=0.0
!      end do
!      end do
!      end do
!!lhl100913
!      do k=1,25
!      do j=53,57
!      do i=283,288
!         vit_global(i,j,k)=1.0
!      end do
!      end do
!      end do

!      do k=1,km
!       vit_global(70,186,k)=0.0
!       vit_global(70,187,k)=0.0
!       vit_global(70,188,k)=0.0
!       vit_global(71,187,k)=0.0
!       vit_global(71,188,k)=0.0
!       vit_global(72,187,k)=0.0
!      end do
!
      do k=1,km
!
      do j=1,jmt_global-1
      do i=2,imt
      ivk(i,j)=vit_global(i-1,j,k)*vit_global(i-1,j+1,k)*vit_global(i,j,k)*vit_global(i,j+1,k)
      enddo
      ivk(1,j)=ivk(imt-1,j)
      enddo
      do i=1,imt
      ivk(i,jmt_global)=0
      enddo
!
      do j=2,jmt_global-1
      do i=2,imt-1
#if (defined D_PRECISION)
      if(ivk (i,j).gt.dmax1(ivk (i+1,j),ivk (i-1,j),ivk (i,j+1),ivk (i,j-1)))then
#else
      if(ivk (i,j).gt.max(ivk (i+1,j),ivk (i-1,j),ivk (i,j+1),ivk (i,j-1)))then
#endif
      ivk (i,j)=0
      endif
      enddo
      enddo
!
      do j=2,jmt_global
      do i=1,imt-1
#if (defined D_PRECISION)
      vit_global(i,j,k)=dmin1(1.D0,ivk(i,j-1)+ivk(i,j)+ivk(i+1,j-1)+ivk(i+1,j))
#else
      vit_global(i,j,k)=min(1.,ivk(i,j-1)+ivk(i,j)+ivk(i+1,j-1)+ivk(i+1,j))
#endif
      enddo
      vit_global(imt,j,k)=vit_global(2,j,k)
      enddo
!
      enddo

!lhl060629
      end if
      call mpi_barrier(mpi_comm_ocn,ierr)
      call mpi_bcast(vit_global,imt*jmt_global*km,MPI_PR,0,mpi_comm_ocn,ierr)
      call mpi_bcast(basin,imt*jmt_global,mpi_integer,0,mpi_comm_ocn,ierr)
      call mpi_bcast(dyr_global,jmt_global,MPI_PR,0,mpi_comm_ocn,ierr)
      call mpi_bcast(zkp,km+1,MPI_PR,0,mpi_comm_ocn,ierr)
      call mpi_bcast(ahv_back,jmt_global,MPI_PR,0,mpi_comm_ocn,ierr)
      call mpi_barrier(mpi_comm_ocn,ierr)
!--------------------------------------------------------------
!     Decompose model domain for the each MPI task.
!--------------------------------------------------------------
      call boundary
#else
#if (!defined CDFIN)
      OPEN (71,FILE ='INDEX.DATA',FORM ='UNFORMATTED',STATUS ='OLD')
      READ (71)VIT_IN,DYR_IN,ZKP_IN
      CLOSE (71)
#else
!----------------------------------------------------
! Open netCDF file.
!----------------------------------------------------
      iret=nf_open('INDEX.DATA',nf_nowrite,ncid)
      call check_err (iret)
!----------------------------------------------------
!   Retrieve data
!----------------------------------------------------
      start3(1)=1 ; count3(1)=imt
      start3(2)=1 ; count3(2)=jmt_global
      start3(3)=1 ; count3(3)=km

      iret=nf_get_vara_double(ncid,   5,start3,count3,vit_in)
      call check_err (iret)
!
      start1(1)=1 ; count1(1)=jmt_global
      iret=nf_get_vara_double(ncid,   6,start1,count1,dyr_in)
      call check_err (iret)
!
      start1(1)=1 ; count1(1)=kmp1
      iret=nf_get_vara_double(ncid,   7,start1,count1,zkp_in)
      call check_err (iret)

      iret = nf_close (ncid)
      call check_err (iret)
!
!-------------------------------------------------------
!    Read Basin Index Field
!-------------------------------------------------------
      iret=nf_open('BASIN.nc',nf_nowrite,ncid)
      call check_err (iret)
!
      start2(1)=1 ; count2(1)=imt
      start2(2)=1 ; count2(2)=jmt_global
      iret=nf_get_vara_int(ncid,   3,start2,count2,basin)
      call check_err (iret)
      iret = nf_close (ncid)
      call check_err (iret)
!
#endif
       vit=vit_in
       dyr=dyr_in
       zkp=zkp_in

!Yu
      do k=1,km
      do j=1,jst
      do i=1,imt
         vit(i,j,k)=0.0
      end do
      end do
      end do
!Yu

!lhl100913
!      do k=1,5
!      do j=123,126
!      do i=143,144
!         vit(i,j,k)=1.0
!      end do
!      end do
!      end do
!!lhl100913

!lhl060629
      do k=1,km
      do j=1,jmm
      do i=2,imt
      ivk(i,j)=vit(i-1,j,k)*vit(i-1,j+1,k)*vit(i,j,k)*vit(i,j+1,k)
      enddo
      ivk(1,j)=ivk(imm,j)
      enddo
      do i=1,imt
      ivk(i,jmt)=0.D0
      enddo
!
      do j=2,jmm
      do i=2,imm
#if (defined D_PRECISION)
      if(ivk (i,j).gt.dmax1(ivk (i+1,j),ivk (i-1,j),ivk (i,j+1),ivk (i,j-1)))then
#else
      if(ivk (i,j).gt.max(ivk (i+1,j),ivk (i-1,j),ivk (i,j+1),ivk (i,j-1)))then
#endif
!      print*,i,j,k
      ivk (i,j)=0.D0
      endif
      enddo
      enddo
!
      do j=2,jmt
      do i=1,imm
#if (defined D_PRECISION)
      vit(i,j,k)=dmin1(1.D0,ivk(i,j-1)+ivk(i,j)+ivk(i+1,j-1)+ivk(i+1,j))
#else
      vit(i,j,k)=min(1.0,ivk(i,j-1)+ivk(i,j)+ivk(i+1,j-1)+ivk(i+1,j))
#endif
      enddo
      vit(imt,j,k)=vit(2,j,k)
      enddo
      enddo
!
!lhl060629

      DO K = 1,KM
         DO J = 1,JMM
            DO I = 2,IMT
               VIV (I,J,K)= VIT (I -1,J,K)* VIT (I -1,J +1,K)* VIT (I,J,&
                   K)* VIT (I,J +1,K)
            END DO
            VIV (1,J,K)= VIV (IMM,J,K)
         END DO
         DO I = 1,IMT
            VIV (I,JMT,K)= 0.0
         END DO
      END DO
!
!      do k=1,km
!      do j=2,jmt_global-1
!      do i=2,imt-1
!      if ((vit(i,j,k).gt.0.5).and.(vit(i-1,j,k).lt.0.5).and.(vit(i+1,j,k).lt.0.5)) print*,i,j
!      enddo
!      enddo
!      enddo
!      do k=1,km
!      do j=2,jmt_global-1
!      do i=2,imt-1
!      if ((vit(i,j,k).gt.0.5).and.(viv(i,j,k)+vit(i+1,j,k)+viv(i,j-1,k)+vit(i+1,j-1,k).lt.0.5)) print*,i,j
!      enddo
!      enddo
!      enddo
#endif

!--------------------------------------------------------------
!     VERTICAL LAYERED PARAMETERS
!--------------------------------------------------------------
!     ZKP    DEPTHS OF BOX BOTTOMS ON MODEL GRID "T" BOXES (M)
!            (0 -25 -50 -75 ... -5600)
!     ZKT    DEPTHS OF BOX BOTTOMS ON MODEL GRID "W" BOXES (M)
!     DZP    D(ZKP)
!     ODZP   1/DZP
!     ODZT   1/DZT

#if (defined NETCDF) || (defined ALL)
      lev1=zkp
#endif
      DO K = 1,KM
         DZP (K) = ZKP (K) - ZKP (K +1)
         ZKT (K) = (ZKP (K) + ZKP (K +1))*0.5D0
!lhl
#if (defined NETCDF) || (defined ALL)
         lev (k)= (ZKP (K) + ZKP (K +1))*0.5D0
#endif
!lhl
         ODZP (K)= 1.0D0/ DZP (K)
      END DO

      ODZT (1)= 2.0D0* ODZP (1)
      DO K = 2,KM
         ODZT (K)= 1.0D0/ (ZKT (K -1) - ZKT (K))
      END DO


#ifdef SPMD
!--------------------------------------------------------------
!     SOME CONSTANTS RELATED TO THICKNESS
!--------------------------------------------------------------
      DO J = 1,JMT
         DO I = 1,IMT
            ITNU (I,J)= 0
            DO K = 1,KM
               IF (VIT_global (I,J_global(j),K) > 0.0)ITNU (I,J)= ITNU (I,J) +1
            END DO
         END DO
      END DO
!
!lhl1204
       DO J = 1,JMT
          DO I = 1,IMT
             NA(I,J)=MAX(1,ITNU (I,J)-1)
          END DO
       END DO
!lhl1204
!
      DO J = 1,JMT_global
         DO I = 1,IMT
            ABCD = 0.0
            DO K = 1,KM
               ABCD = ABCD+ VIT_global (I,J,K)* DZP (K)
            END DO
            WORK_global (I,J)= ABCD
            IF (ABCD > 0.0)THEN
               OHBT_global (I,J)= 1.0D0/ ABCD
            ELSE
               OHBT_global (I,J)= 0.0D0
            END IF
         END DO
      END DO

      DO J = 1,JMT_global
         DO I = 1,IMT
            ABCD = 0.0
            DO K = 1,KM
               ABCD = ABCD+ VIV_global (I,J,K)* DZP (K)
            END DO
            DZPH_global (I,J)= ABCD
            IF (ABCD > 0.0)THEN
               OHBU_global (I,J)= 1.0D0/ ABCD
            ELSE
               OHBU_global (I,J)= 0.0D0
            END IF
         END DO
      END DO

!--------------------------------------------------------------
!     J-DEPENDENT PARAMETER
!--------------------------------------------------------------
!     NORTH POLAR WILL BE TREATED AS AN ISLAND
!     CO-LATITUDE: 0 TO 169.2706 (FROM 90N TO 79.2706S)
      WKJ_global(1)= 0.0D0
      DO J = 2,JMT_global
         WKJ_global(J)= WKJ_global(J -1) + DYR_global(J -1)* TORAD
      END DO
!
!lhl060506
      DO J = 1,JMT_global
         FF1_global (J) = 2.0D0* OMEGA * COS (WKJ_global (J))
      END DO
!lhl060506
!
      DO J = 1,JMM_global
         AJQ = 90.0D0-(WKJ_global (J) + WKJ_global(J +1))*0.5D0/ TORAD
#if (defined D_PRECISION)
         SNLAT_global(J)= DSIGN (1.0D0, AJQ)
#else
         SNLAT_global(J)= SIGN (1.0, AJQ)
#endif
      END DO

      SNLAT_global(JMT_global)= SNLAT_global(JMM_global)
!lhl
#if (defined NETCDF) || (defined ALL)
      lat(1)=90.
      DO J = 2,jmt_global
         lat(j)=  lat(j-1)-dyr_global(j)
      END DO
      lat_r8=lat
#endif

      DO J = 1,JMT_global
         COST_global (J) = COS(WKJ_global(J))
      END DO

      DO J = 2,JMT_global
         SINT_global(J) = SIN(WKJ_global (J))
      END DO
      SINT_global (1) = 1.0D-25

      DO J = 1,JMT_global
         DYR_global(J) = RADIUS * DYR_global (J)* TORAD
      END DO

      DO J = 2,JMT_global
         DYT_global (J) = (DYR_global(J) + DYR_global(J -1))*0.5D0
      END DO

      DYT_global(1) = DYT_global (2)

      DO J = 1,JMM_global
         YU = (WKJ_global (J) + WKJ_global (J +1))*0.5D0
         SINU_global (J) = SIN(YU)
         COSU_global(J) = COS (YU)
         FF_global (J) = 2.0D0* OMEGA * COS (YU)
         CURU = COS(YU)/ SIN (YU)/ RADIUS
         CV1_global (J) = 1.0D0/ (RADIUS * RADIUS) - CURU * CURU
         CV2_global (J) = CURU / (SIN (YU)* DX)
      END DO

      SINU_global (JMT_global)= SINU_global (JMM_global)
      COSU_global (JMT_global)= COSU_global (JMM_global)
      FF_global (JMT_global) = FF_global (JMM_global)
      CV1_global (JMT_global) = CV1_global (JMM_global)
      CV2_global (JMT_global) = CV2_global (JMM_global)

      DO J = 1,JMT_global
#if ( defined SMAG)
! parameter for Smagrinsky scheme.
!lhl         CXT_global(J) = (KARMAN * min (SINT_global(J)* DX,DYR_global(J))/ PI)**2
!lhl         CXU_global(J) = (KARMAN * min (SINU_global(J)* DX,DYR_global(J))/ PI)**2
!lhl         CYT_global(J) = (KARMAN * min (SINT_global(J)* DX,DYR_global(J))/ PI)**2
!lhl         CYU_global(J) = (KARMAN * min (SINU_global(J)* DX,DYR_global(J))/ PI)**2
         CXT_global(J) = (KARMAN * SINT_global(J)* DX)**2*0.5D0
         CXU_global(J) = (KARMAN * SINU_global(J)* DX)**2*0.5D0
         CYT_global(J) = (KARMAN * DYR_global(J))**2*0.5D0
         CYU_global(J) = (KARMAN * DYR_global(J))**2*0.5D0
#endif
         OUY_global(J) = 1.0D0/ DYR_global(J)
         OTX_global(J) = 1.0D0/ (SINT_global(J)* DX)
         OUX_global(J) = 1.0D0/ (SINU_global(J)* DX)
         SOTX_global(J) = OTX_global(J)* OTX_global(J)
         SOUX_global(J) = OUX_global(J)* OUX_global(J)
      END DO


      DO J = 1,JMM_global
         R1A_global(J)= SINT_global(J) / (SINU_global(J)* DYR_global(J))*0.25D0
         R1B_global(J)= SINT_global(J +1)/ (SINU_global(J)* DYR_global(J))*0.25D0
      END DO


      R1A_global (JMT_global)= R1A_global (JMM_global)
      R1B_global (JMT_global)= R1B_global(JMM_global)

      DO J = 2,JMT_global
         R2A_global (J)= SINU_global (J) / (SINT_global(J)* DYT_global(J))*0.25D0
         R2B_global (J)= SINU_global (J-1)/(SINT_global(J)* DYT_global(J))*0.25D0
      END DO

      R2A_global (1)= R2A_global (2)
      R2B_global (1)= R2B_global (2)

!XC
#if (defined TSPAS)
      DO J = 2,JMT_global
         dtdy_global(j)=dts/dyr_global(j)
         dtdx_global(j)=dts*otx_global(j)
         RAA_global(j) =4*R2A_global(j)*dyt_global(j)
         RBB_global(j) =4*R2B_global(j)*dyt_global(j)
      ENDDO

        dtdy_global(1)=dtdy_global(2)
        dtdx_global(1)=dtdx_global(2)
        RAA_global(1) =RAA_global(2)
        RBB_global(1) =RBB_global(2)
#endif
!XC

      DO J = 1,JMM_global
         R1C_global (J)= SINT_global (J )/ (DYT_global (J )*SINU_global (J)*DYR_global(J))
         R1D_global (J)= SINT_global (J+1)/(DYT_global(J+1)*SINU_global(J)* DYR_global(J))
      END DO

      R1C_global(JMT_global)= R1C_global(JMM_global)
      R1D_global(JMT_global)= R1D_global(JMM_global)

      DO J = 2,JMT_global
         R2C_global(J)= SINU_global(J -1)/ (DYR_global(J -1)* SINT_global(J)* DYT_global(J))
         R2D_global(J)= SINU_global(J )/ (DYR_global(J )* SINT_global(J)* DYT_global(J))
      END DO

      R2C_global(1)= R2C_global(2)
      R2D_global(1)= R2D_global(2)

#if ( defined SMAG)
#if ( defined SMAG_FZ )
      DO J = 2,JMT_global
         R1E_global(J)= 0.5* SINT_global(J)/ (DYT_global(J)* SINU_global(J ))
         R1F_global(J)= 0.5* SINT_global(J)/ (DYT_global(J)* SINU_global(J -1))
      END DO

      R1E_global(1)= R1E_global(2)
      R1F_global(1)= R1F_global(2)
#else
      DO J = 2,JMM_global
         R1E_global(J)= SINU_global(J)/ (DYT_global(J) + DYT_global(J +1))/ SINU_global(J +1)
         R1F_global(J)= SINU_global(J)/ (DYT_global(J) + DYT_global(J +1))/ SINU_global(J -1)
      END DO

      R1E_global(1)= R1E_global(2)
      R1F_global(1)= R1F_global(2)
      R1E_global(JMT)= R1E_global(JMM)
      R1F_global(JMT)= R1F_global(JMM)
#endif

#if ( defined SMAG_FZ )
!-new
      DO J = 2,JMT_global
         R2E_global(J)= SINU_global(J)* SINU_global(J)/ DYR_global(J)/ &
                        (SINT_global(J)* SINT_global(J))
         R2F_global(J)= SINU_global(J -1)* SINU_global(J -1)/ DYR_global(J)/&
                        (SINT_global(J)* SINT_global(J))
      END DO

      R2E_global(1)= R2E_global(2)
      R2F_global(1)= R2F_global(2)
!-new
#else
      DO J = 2,JMM_global
         R2E_global(J)= SINU_global(J+1)*SINU_global(J +1)/(DYR_global(J)+DYR_global(J+1))/(  &
                   SINU_global(J)* SINU_global(J))
         R2F_global(J)= SINU_global(J-1)*SINU_global(J-1)/(DYR_global(J)+DYR_global(J+1))/ (   &
                   SINU_global(J)* SINU_global(J))
      END DO
      R2E_global(1)= R2E_global(2)
      R2F_global(1)= R2F_global(2)
      R2E_global(JMT_global)= R2E_global(JMM_global)
      R2F_global(JMT_global)= R2F_global(JMM_global)
#endif

#if ( defined SMAG_FZ )
      DO J = 2,JMT_global
         R3E_global(J)= SINU_global(J)/ DYR_global(J)/ SINT_global(J)
         R3F_global(J)= SINU_global(J -1)/ DYR_global(J)/ SINT_global(J)
      END DO
      R3E_global(1)= R3E_global(2)
      R3F_global(1)= R3F_global(2)
#else
      DO J = 2,JMM_global
         R3E_global(J)= SINT_global(J +1)/ (DYR_global(J) + DYR_global(J +1))/  &
                        (SINU_global(J)* SINU_global(J))
         R3F_global(J)= SINT_global(J -1)/ (DYR_global(J) + DYR_global(J +1))/ &
                        (SINU_global(J)* SINU_global(J))
      END DO
      R3E_global(1)= R3E_global(2)
      R3F_global(1)= R3F_global(2)
      R3E_global(JMT_global)= R3E_global(JMM_global)
      R3F_global(JMT_global)= R3F_global(JMM_global)
#endif

#if ( defined SMAG_FZ )
      DO J = 1,JMT_global
         R4E_global(J)= COSU_global(J)/ (RADIUS * SINU_global(J))
      END DO

#else
      DO J = 1,JMT_global
         R4E_global(J)= COSU_global(J)/ (RADIUS * SINU_global(J))
      END DO

#endif
#endif

!--------------------------------------------------------------

      DO J = 1,JMT_global
         EPS = 0.5D0* FF_global(J)* DTC
         EPEA_global(J) = 1.0D0/ (1.0D0+ EPS * EPS)
         EPEB_global(J) = EPS / (1.0D0+ EPS * EPS)
         EPS = FF_global(J)* DTC
         EPLA_global(J) = 1.0D0/ (1.0D0+ EPS * EPS)
         EPLB_global(J) = EPS / (1.0D0+ EPS * EPS)
         EPS = 0.5D0* FF_global(J)* DTB
         EBEA_global(J) = 1.0D0/ (1.0D0+ EPS * EPS)
         EBEB_global(J) = EPS / (1.0D0+ EPS * EPS)
         EPS = FF_global(J)* DTB
         EBLA_global(J) = 1.0D0/ (1.0D0+ EPS * EPS)
         EBLB_global(J) = EPS / (1.0D0+ EPS * EPS)
      END DO


      DO J = 1,JMT_global
         DO I = 1,IMT
            HBY_global(I,J)= 0.0D0
            HBX_global(I,J)= 0.0D0
         END DO
      END DO

      DO J = 2,JMM_global
         DO I = 2,IMT
            HBY_global(I,J)= VIV_global(I,J,1)* OUY_global(J)*0.5D0* &
            (WORK_global(I,J +1) - WORK_global(I,J) + WORK_global(I -1,J +1) - WORK_global(I -1,J))
            HBX_global(I,J)= VIV_global(I,J,1)* OUX_global(J)*0.5D0* &
            (WORK_global(I,J) - WORK_global(I -1,J) + WORK_global(I,J +1) - WORK_global(I -1,J +1))
         END DO
         HBY_global(1,J)= HBY_global(IMM,J)
         HBX_global(1,J)= HBX_global(IMM,J)
      END DO

      DO J = 1,JMT_global
         DXDYU_global(J)= DX * DYR_global(J)* SINU_global(J)
         DXDYT_global(J)= DX * DYT_global(J)* SINT_global(J)
      END DO
      ASEA = 0.D0
      DO J = 1,JMT_global
         DO I = 2,IMM
            ASEA = ASEA + DYT_global(J)* SINT_global(J)* VIT_global(I,J,1)
         END DO
      END DO

      VSEA = 0.D0
      DO K = 1,KM
         DO J = 1,JMT_global
            DO I = 2,IMM
               VSEA = VSEA + DZP(K)* DXDYT_global(J)* VIT_global(I,J,K)
            END DO
         END DO
      END DO

!--------------------------------------------------------------
!     Shortwave penetration is a double exponential as follows
!--------------------------------------------------------------

#if (defined SOLAR)

      rpart = 0.58D0
      efold1 = 0.35D0
      efold2 = 23.0D0
      rscl1 = 1.0D0/ efold1
      rscl2 = 1.0D0/ efold2

      DO k = 1,kmm1
#if (defined D_PRECISION)
         swarg1 = dmax1 (ZKP (k +1)* rscl1, -70.0D0)
         swarg2 = dmax1 (ZKP (k +1)* rscl2, -70.0D0)
#else
         swarg1 = max (ZKP (k +1)* rscl1, -70.0)
         swarg2 = max (ZKP (k +1)* rscl2, -70.0)
#endif
         pen (k) = rpart * exp (swarg1) + (1.0D0- rpart)* exp (swarg2)
      END DO


      DO k = 1,kmm1
         pen (k) = pen (k)* OD0CP
      END DO

!      write(*,'(5D25.15)')pen
#endif

!--------------------------------------------------------------
!     compute boundary of area where vertical mixing coefficients
!     are depended on Richardson number.
!--------------------------------------------------------------

      ABCD = 60.0D0* TORAD
      DO J = 2,JMM_global
         IF (ASIN (SINU_global(J -1)) < ABCD.AND.ASIN (SINU_global(J)) >= ABCD) RUST = J
         IF (ASIN (SINU_global(J +1)) < ABCD.AND.ASIN (SINU_global(J)) >= ABCD) RUEND = J
         IF (ASIN (SINT_global(J -1)) < ABCD.AND.ASIN (SINT_global(J)) >= ABCD) RTST = J
         IF (ASIN (SINT_global(J +1)) < ABCD.AND.ASIN (SINT_global(J)) >= ABCD) RTEND = J
      END DO
!Yu

#if  ( defined SMAG)
      call global_to_local_1d(cxt_global,cxt)
      call global_to_local_1d(cxu_global,cxu)
      call global_to_local_1d(cyt_global,cyt)
      call global_to_local_1d(cyu_global,cyu)
      call global_to_local_1d(r1e_global,r1e)
      call global_to_local_1d(r1f_global,r1f)
      call global_to_local_1d(r2e_global,r2e)
      call global_to_local_1d(r2f_global,r2f)
      call global_to_local_1d(r3e_global,r3e)
      call global_to_local_1d(r3f_global,r3f)
      call global_to_local_1d(r4e_global,r4e)
      call global_to_local_1d(r4f_global,r4f)

#endif
      call global_to_local_1d(cosu_global,cosu)
      call global_to_local_1d(cost_global,cost)

      call global_to_local_1d(oux_global,oux)
      call global_to_local_1d(ouy_global,ouy)
      call global_to_local_1d(otx_global,otx)
      call global_to_local_1d(sotx_global,sotx)
      call global_to_local_1d(soux_global,soux)
      call global_to_local_1d(ff_global,ff)
!lhl
      call global_to_local_1d(ff1_global,ff1)
!lhl
#ifdef TIDE
!lhl20100801
      call global_to_local_1d(lat_r8,lat_tidal)
!lhl20100801
#endif
      call global_to_local_1d(cv1_global,cv1)
      call global_to_local_1d(cv2_global,cv2)
      call global_to_local_1d(snlat_global,snlat)
      call global_to_local_1d(sinu_global,sinu)
      call global_to_local_1d(sint_global,sint)
      call global_to_local_1d(wkj_global,wkj)
      call global_to_local_1d(dyr_global,dyr)
      call global_to_local_1d(dyt_global,dyt)
      call global_to_local_1d(dxdyu_global,dxdyu)
      call global_to_local_1d(dxdyt_global,dxdyt)
      call global_to_local_1d(r1a_global,r1a)
      call global_to_local_1d(r1b_global,r1b)
      call global_to_local_1d(r2a_global,r2a)
      call global_to_local_1d(r2b_global,r2b)
!XC
#if (defined TSPAS)
      call global_to_local_1d(dtdy_global,dtdy)
      call global_to_local_1d(dtdx_global,dtdx)
      call global_to_local_1d(raa_global ,raa)
      call global_to_local_1d(rbb_global ,rbb)
#endif
!XC
      call global_to_local_1d(r1c_global,r1c)
      call global_to_local_1d(r1d_global,r1d)
      call global_to_local_1d(r2c_global,r2c)
      call global_to_local_1d(r2d_global,r2d)
      call global_to_local_1d(ebea_global,ebea)
      call global_to_local_1d(ebeb_global,ebeb)
      call global_to_local_1d(ebla_global,ebla)
      call global_to_local_1d(eblb_global,eblb)
      call global_to_local_1d(epea_global,epea)
      call global_to_local_1d(epeb_global,epeb)
      call global_to_local_1d(epla_global,epla)
      call global_to_local_1d(eplb_global,eplb)
!Yu
      call global_to_local_4d(hbx_global,hbx,1,1)
      call global_to_local_4d(hby_global,hby,1,1)
      call global_to_local_4d(ohbt_global,ohbt,1,1)
      call global_to_local_4d(ohbu_global,ohbu,1,1)
      call global_to_local_4d(dzph_global,dzph,1,1)
!
!Yu
      call global_to_local_4d(vit_global,vit,km,1)
      call global_to_local_4d(viv_global,viv,km,1)
      call global_to_local_4d(ahv_back_global,ahv_back,1,1)

#else
!--------------------------------------------------------------
!     SOME CONSTANTS RELATED TO THICKNESS
!--------------------------------------------------------------
      DO J = 1,JMT
         DO I = 1,IMT
            ITNU (I,J)= 0
            DO K = 1,KM
               IF (VIT (I,J,K) > 0.0)ITNU (I,J)= ITNU (I,J) +1
            END DO
         END DO
      END DO
!lhl1204
       DO J = 1,JMT
          DO I = 1,IMT
             NA(I,J)=MAX(1,ITNU (I,J)-1)
          END DO
       END DO
!lhl1204

      DO J = 1,JMT
         DO I = 1,IMT
            ABCD = 0.0D0
            DO K = 1,KM
               ABCD = ABCD+ VIT (I,J,K)* DZP (K)
            END DO
            WORK (I,J)= ABCD
            IF (ABCD > 0.0)THEN
               OHBT (I,J)= 1.0D0/ ABCD
            ELSE
               OHBT (I,J)= 0.0D0
            END IF
         END DO
      END DO

      DO J = 1,JMT
         DO I = 1,IMT
            ABCD = 0.0
            DO K = 1,KM
               ABCD = ABCD+ VIV (I,J,K)* DZP (K)
            END DO
            DZPH (I,J)= ABCD
            IF (ABCD > 0.0)THEN
               OHBU (I,J)= 1.0D0/ ABCD
            ELSE
               OHBU (I,J)= 0.0D0
            END IF
         END DO
      END DO

!--------------------------------------------------------------
!     J-DEPENDENT PARAMETER
!--------------------------------------------------------------
!     NORTH POLAR WILL BE TREATED AS AN ISLAND
!     CO-LATITUDE: 0 TO 169.2706 (FROM 90N TO 79.2706S)
      WKJ (1)= 0.0
      DO J = 2,JMT
         WKJ (J)= WKJ (J -1) + DYR (J -1)* TORAD
      END DO
!lhl060506
      DO J = 1,JMT
         FF1 (J) = 2.0D0* OMEGA * COS (WKJ(J))
      END DO
!lhl060506

      DO J = 1,JMM
         AJQ = 90.0D0- (WKJ (J) + WKJ (J +1))*0.5D0/ TORAD
#if (defined D_PRECISION)
         SNLAT(J)= DSIGN (1.0D0, AJQ)
#else
         SNLAT(J)= SIGN (1.0, AJQ)
#endif
      END DO

      SNLAT (JMT)= SNLAT (JMM)
!lhl
#if (defined NETCDF) || (defined ALL)
      lat(1)=90.
      DO J = 2,jmt
         lat(j)=  lat(j-1)-dyr(j)
      END DO
#endif
!lhl
      DO J = 1,JMT
         COST (J) = COS (WKJ (J))
      END DO

      DO J = 2,JMT
         SINT (J) = SIN (WKJ (J))
      END DO
      SINT (1) = 1.0D-25

      DO J = 1,JMT
         DYR (J) = RADIUS * DYR (J)* TORAD
      END DO

      DO J = 2,JMT
         DYT (J) = (DYR (J) + DYR (J -1))*0.5D0
      END DO

      DYT (1) = DYT (2)

      DO J = 1,JMM
         YU = (WKJ (J) + WKJ (J +1))*0.5D0
         SINU (J) = SIN (YU)
         COSU(J)=COS(YU)
#if ( defined SMAG)
         COSU (J) = COS (YU)
#endif
         FF (J) = 2.0D0* OMEGA * COS (YU)
         CURU = COS (YU)/ SIN (YU)/ RADIUS
         CV1 (J) = 1.0D0/ (RADIUS * RADIUS) - CURU * CURU
         CV2 (J) = CURU / (SIN (YU)* DX)
      END DO

      SINU (JMT)= SINU (JMM)
      COSU (JMT)=COSU(JMM)
#if ( defined SMAG)
      COSU (JMT)= COSU (JMM)
#endif
      FF (JMT) = FF (JMM)
      CV1 (JMT) = CV1 (JMM)
      CV2 (JMT) = CV2 (JMM)

      DO J = 1,JMT
#if ( defined SMAG)
! parameter for Smagrinsky scheme.
!lhl         CXT (J) = (KARMAN * min (SINT (J)* DX,DYR (J))/ PI)**2
!lhl         CXU (J) = (KARMAN * min (SINU (J)* DX,DYR (J))/ PI)**2
!lhl         CYT (J) = (KARMAN * min (SINT (J)* DX,DYR (J))/ PI)**2
!lhl         CYU (J) = (KARMAN * min (SINU (J)* DX,DYR (J))/ PI)**2
         CXT(J) = (KARMAN * SINT(J)* DX)**2*0.5D0
         CXU(J) = (KARMAN * SINU(J)* DX)**2*0.5D0
         CYT(J) = (KARMAN * DYR(J)))**2*0.5D0
         CYU(J) = (KARMAN * DYR(J)))**2*0.5D0
#endif
         OUY (J) = 1.0/ DYR (J)
         OTX (J) = 1.0/ (SINT (J)* DX)
         OUX (J) = 1.0/ (SINU (J)* DX)
         SOTX (J) = OTX (J)* OTX (J)
         SOUX (J) = OUX (J)* OUX (J)
      END DO


      DO J = 1,JMM
         R1A (J)= SINT (J) / (SINU (J)* DYR (J))*0.25D0
         R1B (J)= SINT (J +1)/ (SINU (J)* DYR (J))*0.25D0
      END DO

      R1A (JMT)= R1A (JMM)
      R1B (JMT)= R1B (JMM)

      DO J = 2,JMT
         R2A (J)= SINU (J) / (SINT (J)* DYT (J))*0.25D0
         R2B (J)= SINU (J -1)/ (SINT (J)* DYT (J))*0.25D0
      END DO

      R2A (1)= R2A (2)
      R2B (1)= R2B (2)

!XC
#if (defined TSPAS)
      DO J = 2,JMT
         dtdy(j)=dts/dyr(j)
         dtdx(j)=dts*otx(j)
         RAA(j) =4*R2A(j)*dyt(j)
         RBB(j) =4*R2B(j)*dyt(j)
      ENDDO

        dtdy(1)=dtdy(2)
        dtdx(1)=dtdx(2)
        RAA(1) =RAA(2)
        RBB(1) =RBB(2)
#endif
!XC

      DO J = 1,JMM
         R1C (J)= SINT (J )/ (DYT (J )* SINU (J)* DYR (J))
         R1D (J)= SINT (J +1)/ (DYT (J +1)* SINU (J)* DYR (J))
      END DO

      R1C (JMT)= R1C (JMM)
      R1D (JMT)= R1D (JMM)

      DO J = 2,JMT
         R2C (J)= SINU (J -1)/ (DYR (J -1)* SINT (J)* DYT (J))
         R2D (J)= SINU (J )/ (DYR (J )* SINT (J)* DYT (J))
      END DO

      R2C (1)= R2C (2)
      R2D (1)= R2D (2)

#if ( defined SMAG)
#if ( defined SMAG_FZ )
      DO J = 2,JMT
         R1E (J)= 0.5D0* SINT (J)/ (DYT (J)* SINU (J ))
         R1F (J)= 0.5D0* SINT (J)/ (DYT (J)* SINU (J -1))
      END DO

      R1E (1)= R1E (2)
      R1F (1)= R1F (2)
#else
      DO J = 2,JMM
         R1E (J)= SINU (J)/ (DYT (J) + DYT (J +1))/ SINU (J +1)
         R1F (J)= SINU (J)/ (DYT (J) + DYT (J +1))/ SINU (J -1)
      END DO

      R1E (1)= R1E (2)
      R1F (1)= R1F (2)
      R1E (JMT)= R1E (JMM)
      R1F (JMT)= R1F (JMM)
#endif

#if ( defined SMAG_FZ )
      DO J = 2,JMT
         R2E (J)= SINU (J)* SINU (J)/ DYR (J)/ (SINT (J)* SINT (J))
         R2F (J)= SINU (J -1)* SINU (J -1)/ DYR (J)/ (SINT (J)* SINT (J))
      END DO

      R2E (1)= R2E (2)
      R2F (1)= R2F (2)
!-new
#else
      DO J = 2,JMM
         R2E (J)= SINU (J +1)* SINU (J +1)/ (DYR (J) + DYR (J +1))/ (   &
                   SINU (J)* SINU (J))
         R2F (J)= SINU (J -1)* SINU (J -1)/ (DYR (J) + DYR (J +1))/ (   &
                   SINU (J)* SINU (J))
      END DO
      R2E (1)= R2E (2)
      R2F (1)= R2F (2)
      R2E (JMT)= R2E (JMM)
      R2F (JMT)= R2F (JMM)
#endif

#if ( defined SMAG_FZ )
      DO J = 2,JMT
         R3E (J)= SINU (J)/ DYR (J)/ SINT (J)
         R3F (J)= SINU (J -1)/ DYR (J)/ SINT (J)
      END DO
      R3E (1)= R3E (2)
      R3F (1)= R3F (2)
#else
      DO J = 2,JMM
         R3E (J)= SINT (J +1)/ (DYR (J) + DYR (J +1))/ (SINU (J)* SINU (J))
         R3F (J)= SINT (J -1)/ (DYR (J) + DYR (J +1))/ (SINU (J)* SINU (J))
      END DO
      R3E (1)= R3E (2)
      R3F (1)= R3F (2)
      R3E (JMT)= R3E (JMM)
      R3F (JMT)= R3F (JMM)
#endif

#if ( defined SMAG_FZ )
      DO J = 1,JMT
         R4E (J)= COSU (J)/ (RADIUS * SINU (J))
      END DO

#else
      DO J = 1,JMT
         R4E (J)= COSU (J)/ (RADIUS * SINU (J))
      END DO

#endif
#endif

      DO J = 1,JMT
         EPS = 0.5D0* FF (J)* DTC
         EPEA (J) = 1.0D0/ (1.0D0+ EPS * EPS)
         EPEB (J) = EPS / (1.0D0+ EPS * EPS)
         EPS = FF (J)* DTC
         EPLA (J) = 1.0D0/ (1.0D0+ EPS * EPS)
         EPLB (J) = EPS / (1.0D0+ EPS * EPS)
         EPS = 0.5D0* FF (J)* DTB
         EBEA (J) = 1.0D0/ (1.0D0+ EPS * EPS)
         EBEB (J) = EPS / (1.0D0+ EPS * EPS)
         EPS = FF (J)* DTB
         EBLA (J) = 1.0D0/ (1.0D0+ EPS * EPS)
         EBLB (J) = EPS / (1.0D0+ EPS * EPS)
      END DO


      DO J = 1,JMT
         DO I = 1,IMT
            HBY (I,J)= 0.0D0
            HBX (I,J)= 0.0D0
         END DO
      END DO

      DO J = 2,JMM
         DO I = 2,IMT
            HBY (I,J)= VIV (I,J,1)* OUY (J)*0.5D0* &
            (WORK (I,J +1) - WORK (I,J) + WORK (I -1,J +1) - WORK (I -1,J))
            HBX (I,J)= VIV (I,J,1)* OUX (J)*0.5D0* &
            (WORK (I,J) - WORK (I -1,J) + WORK (I,J +1) - WORK (I -1,J +1))
         END DO
         HBY (1,J)= HBY (IMM,J)
         HBX (1,J)= HBX (IMM,J)
      END DO

      DO J = 1,JMT
         DXDYU (J)= DX * DYR (J)* SINU (J)
         DXDYT (J)= DX * DYT (J)* SINT (J)
      END DO

      ASEA = 0.D0
      DO J = 1,JMT
         DO I = 2,IMM
            ASEA = ASEA + DYT (J)* SINT (J)* VIT (I,J,1)
         END DO
      END DO

      VSEA = 0.D0
      DO K = 1,KM
         DO J = 1,JMT
            DO I = 2,IMM
               VSEA = VSEA + DZP (K)* DXDYT (J)* VIT (I,J,K)
            END DO
         END DO
      END DO

!--------------------------------------------------------------
!     Shortwave penetration is a double exponential as follows
!--------------------------------------------------------------

#if (defined SOLAR)
      rpart = 0.58D0
      efold1 = 0.35D0
      efold2 = 23.0D0
      rscl1 = 1.0D0/ efold1
      rscl2 = 1.0D0/ efold2

      DO k = 1,kmm1
#if (defined D_PRECISION)
         swarg1 = dmax1 (ZKP (k +1)* rscl1, -70.0D0)
         swarg2 = dmax1 (ZKP (k +1)* rscl2, -70.0D0)
#else
         swarg1 = max (ZKP (k +1)* rscl1, -70.0)
         swarg2 = max (ZKP (k +1)* rscl2, -70.0)
#endif
         pen (k) = rpart * exp (swarg1) + (1.0D0- rpart)* exp (swarg2)
      END DO

!      if (mytid==0) print*,pen

      DO k = 1,kmm1
!lhl        pen(k) = pen(k)*solar0*OD0CP
         pen (k) = pen (k)* OD0CP
      END DO

#endif

!--------------------------------------------------------------
!     compute boundary of area where vertical mixing coefficients
!     are depended on Richardson number.
!--------------------------------------------------------------

      ABCD = 60.0D0* TORAD
      DO J = 2,JMM
         IF (ASIN (SINU (J -1)) < ABCD.AND.ASIN (SINU (J)) >= ABCD) RUST = J
         IF (ASIN (SINU (J +1)) < ABCD.AND.ASIN (SINU (J)) >= ABCD) RUEND = J
         IF (ASIN (SINT (J -1)) < ABCD.AND.ASIN (SINT (J)) >= ABCD) RTST = J
         IF (ASIN (SINT (J +1)) < ABCD.AND.ASIN (SINT (J)) >= ABCD) RTEND = J
      END DO


#endif
!lhl060506
!--------------------------------------------------------------
!     Rossby radiu of deformation for Large, Danabasoglu and Doney scheme
!     for tapering isopycnal mixing at surface.
!    C=2.0 m/s, RRD between 15km and 100km
!--------------------------------------------------------------
! at T point
      DO J=1,JMT
      IF (FF(J)/=0.0) THEN
      ABCD=2.D0/ABS(FF(J))
      ELSE
      ABCD=100000.D0
      ENDIF
#ifdef D_PRECISION
      RRD1(J)=DMAX1(15000.D0,DMIN1(100000.D0,ABCD))
#else
      RRD1(J)=MAX(15000.,MIN(100000.0,ABCD))
#endif
      ENDDO
! at U point
      DO J=1,JMT
      IF (FF(J)/=0.0) THEN
      ABCD=2.D0/ABS(FF1(J))
      ELSE
      ABCD=100000.D0
      ENDIF
#ifdef D_PRECISION
      RRD2(J)=DMAX1(15000.D0,DMIN1(100000.D0,ABCD))
#else
      RRD2(J)=MAX(15000.0,MIN(100000.0,ABCD))
#endif
      ENDDO

!       print*,mytid,j_global
!       stop 777
!lhl060506

!lhl--------------------------------------------------------------
!     Set initial laternal diffussion coeffcient
!--------------------------------------------------------------
#if ( defined SMAG)
      DO K = 1,KM
         DO J = 1,JMT
            DO I = 1,IMT
!XC               AH3 (I,J,K)= 0.0
               AH3 (I,J,K)= 1.0D+3
            END DO
         END DO
      END DO

#else
#if (defined BIHAR)
      DO K = 1,KM
         DO J = 1,JMT
            DO I = 1,IMT
               AH3 (I,J,K)= -10.0D+11
            END DO
         END DO
      END DO

#else
      ah3 =0.0d0
!     DO K = 1,10
!        DO J = 1,JMT
!           DO I = 1,IMT
!YU
!          AH3 (I,J,K)= 1.0D+3
!YU
!           END DO
!        END DO
!     END DO

#endif
#endif

!lhl
!--------------------------------------------------------------
!     Set laternal viscosity coeffcient depended on J
!--------------------------------------------------------------

#if ( defined SMAG)

      DO J = 1,JMT
         AM (J)= 5.0D+3
      END DO


#else
#if (defined BIHAR)
      DO J = 1,JMT
!lhl
         AM (J)= -10.0D+11
      END DO

      DO K = 1,KM
         DO J = 1,JMT
            DO I = 1,IMT
#if (defined ACOS)
               AM3 (I,J,K)= AM (J)* SINU (J)* SINU (J)* SINU (J)
#else
               AM3 (I,J,K)= AM (J)
#endif
            END DO
         END DO
      END DO

!!!!!
      DO J = 1,JMT
         AM (J)= AM_EXT
      END DO
      ABCD = 40.0D0* TORAD
      DO J = 1,JMT
         IF (ASIN (SINU (J)) >= ABCD) AM (J)= AM_TRO
      END DO

!!!!!
#else
      DO J = 1,JMT
         AM (J)= AM_EXT
      END DO


      ABCD = 40.0D0* TORAD
      DO J = 1,JMT
         IF (ASIN (SINU (J)) >= ABCD) AM (J)= AM_TRO
      END DO

      DO K = 1,KM
         DO J = 1,JMT
            DO I = 1,IMT
#if (defined ACOS)
               AM3 (I,J,K)= AM (J)* SINU (J)
#else
               AM3 (I,J,K)= AM (J)
#endif
            END DO
         END DO
      END DO

#endif
#endif


!      if (mytid==0)then
!
!        open (111,file="viscosity2.dat",&
!        access="direct",form="unformatted",recl=imt*jmt_global*4)
!        do k=1,km
!        read(111,rec=k) ((am3_io(i,j,k),i=1,imt),j=1,jmt_global)
!        enddo
!       close(111)
!
!! south to north
!      do k = 1,km
!      do j = 1,jmt_global
!      do i = 1,imt
!       am3_tmp(i,j,k)=am3_io(i,jmt_global+1-j,k)*vit_global(i,j,k)
!      end do
!      end do
!      end do
!
!!     to u grids
!      do k = 1,km
!         do j = 1,jmt_global-1
!            do i = 2,imm
!               am3_u (i,j,k)= 0.25d0* (am3_tmp (i,j,k) + am3_tmp (i -1,j,k)&
!                    + am3_tmp (i,j +1,k) + am3_tmp (i -1,j +1,k))*viv_global(i,j,k)
!            end do
!              am3_u(1,j,k)=am3_u(imm,j,k)
!              am3_u(imt,j,k)=am3_u(2,j,k)
!         end do
!            do i = 1,imt
!               am3_u (i,jmt_global,k)= am3_u (i,jmt_global-1,k)
!            end do
!      end do
!
!!      do k = 1,km
!!      do j = 1,jmt_global
!!      do i = 1,imt
!!      am3_u(i,j,k)=(am3_u(i,j,k)+8000.0)*viv_global(i,j,k)
!!      end do
!!      end do
!!      end do
!
!      end if
!!
!       call global_to_local_4d(am3_u,am3,km,1)
!!
#ifdef SPMD
      deallocate(work_global)
!YU   deallocate(oux_global,ouy_global)
!YU   deallocate(otx_global,ff_global)
!YU   deallocate(sotx_global,soux_global)
!YU   deallocate(cv1_global,cv2_global)
!YU   deallocate(snlat_global)
!YU   deallocate(sint_global,sinu_global)
!YU   deallocate(dyr_global,dyt_global)
!YU   deallocate(dxdyu_global,dxdyt_global)
!YU   deallocate(r1a_global,r1b_global)
!YU   deallocate(r2a_global,r2b_global)
!YU   deallocate(r1c_global,r1d_global)
!YU   deallocate(r2c_global,r2d_global)
!YU   deallocate(ebea_global,ebeb_global)
!YU   deallocate(ebla_global,eblb_global)
!YU   deallocate(epea_global,epeb_global)
!YU   deallocate(epla_global,eplb_global)

!YU   deallocate(hbx_global,hby_global)
!YU   deallocate(ohbu_global,ohbt_global)
!YU   deallocate(dzph_global)
!Yu
!YU   deallocate(vit_global,viv_global)

#if  ( defined SMAG)
      deallocate(cxt_global,cxu_global)
      deallocate(cyt_global,cyu_global)
      deallocate(r1e_global,r1f_global)
      deallocate(r2e_global,r2f_global)
      deallocate(r3e_global,r3f_global)
      deallocate(r4e_global,r4f_global)
      deallocate(cost_global,cosu_global)
#endif
     deallocate(cosu_global)
#endif
      if (mytid==0)then
      write(6,*)"END------------GRIDS !"
      endif
      RETURN
      END SUBROUTINE GRIDS


