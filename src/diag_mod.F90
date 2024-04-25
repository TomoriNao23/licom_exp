module diag_mod
#include <def-undef.h>
!-------------------------------------------------------------------------------
!
! Author: Yongqiang YU  ( 1  Dec, 2003)
!
!-------------------------------------------------------------------------------
use param_mod 
use pconst_mod 
use output_mod 
use tracer_mod
use work_mod 
!
      implicit none
      logical :: diag_msf,diag_bsf,diag_budget,diag_mth
!
      contains
!

!     ===================
      SUBROUTINE msf(va)
!     ===================
!
      implicit none
!
      real(r4) :: va(imt,jmt_global,km)
      integer :: nin,nta(jmt_global)
!
      if (diag_msf.eqv..false.) return
!
!     allocate (work_1(imt,jmt_global,km),work_2(imt,jmt_global,km))
!
!
      DO  NIN=1,2
!
      IF(NIN.eq.1)then
!
!$OMP PARALLEL DO PRIVATE (K,J,I)
      do k=1,km
      do j=1,jmt_global
      do i=1,imt
#ifdef SPMD
      work_1(i,j,k)=vit_global(i,j,k)
      work_2(i,j,k)=viv_global(i,j,k)
#else
      work_1(i,j,k)=vit(i,j,k)
      work_2(i,j,k)=viv(i,j,k)
#endif
      enddo
      enddo
      enddo
!
      ELSE
!$OMP PARALLEL DO PRIVATE (K,J,I)
      do k=1,km
      do j=1,jmt_global
      do i=1,imt
         if (basin(i,j)==1.or.basin(i,j)==4) then
#ifdef SPMD
            work_1(i,j,k)=1.0*vit_global(i,j,k)
#else
            work_1(i,j,k)=1.0*vit(i,j,k)
#endif
         else
            work_1(i,j,k)=0.0
         end if
      enddo
      enddo
      enddo
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K=1,KM
!
      DO J=1,JMT_GLOBAL-1
      DO I=2,IMT
          work_2(I,J,K)=work_1(I-1,J,K)*work_1(I-1,J+1,K)*work_1(I,J,K)*work_1(I,J+1,K)
      ENDDO
          work_2(1,J,K)=work_2(IMM,J,K)
      ENDDO
!
      DO I=1,IMT
          work_2(I,JMT_global,K)=0.0
      ENDDO
!
      ENDDO
      ENDIF
!
!$OMP PARALLEL DO PRIVATE (K,J)
      DO K=1,km+1
      DO J=1,jmt_global
           psi(J,K,NIN)=0.0
      ENDDO
      ENDDO
!
      DO K=2,km+1
      DO J=1,jmt_global
      DO I=2,imm
#ifdef SPMD
           psi(J,K,NIN)=psi(j,k,nin)+work_2(i,j,k-1)*va(I,J,K-1)*dzp(k-1)/oux_global(j)*1.0E-6*(-1.0)
#else
           psi(J,K,NIN)=psi(j,k,nin)+work_2(i,j,k-1)*va(I,J,K-1)*dzp(k-1)/oux(j)*1.0E-6*(-1.0)
#endif
      ENDDO
           psi(J,K,NIN)=psi(J,K-1,NIN)+psi(J,K,NIN)
      ENDDO
      ENDDO
!
      DO J=1,jmt_global
          NTA(J)=0
      DO K=1,km
      DO I=2,imm
          IF(work_1(I,J,K).gt.0.5) THEN
             NTA(J)=k
          END IF
      END DO
      END DO
!
      END DO
!
!$OMP PARALLEL DO PRIVATE (J)
      DO J=1,jmt_global
           IF(NTA(J).LT.1) psi(J,1,NIN)=SPVAL
      ENDDO
!
!$OMP PARALLEL DO PRIVATE (J,K)
      DO J=1,jmt_global
      DO K=2,km+1
           IF(NTA(J).LT.(K-1)) psi(J,K,NIN)=SPVAL
      ENDDO
      ENDDO
!
      END DO  
!
!     deallocate (work_1,work_2)
!
      RETURN
      END SUBROUTINE MSF

!
!====================================
      SUBROUTINE BAROSF(ua)
!====================================
      implicit none
      real(r4) :: ua(imt,jmt_global,km)
!
      if (diag_bsf.eqv..false.) return
!
!     allocate (work_1(imt,jmt_global,1),work_2(imt,jmt_global,1))

!$OMP PARALLEL DO PRIVATE (I,J)
      do j=1,jmt_global
      do i=1,imt
      work_1(i,j,1)=0.0
      work_2(i,j,1)=0.0
      enddo
      enddo
!
!$OMP PARALLEL DO PRIVATE (I,J,K)
      do k=1,km
      do j=1,jmt_global
      do i=1,imt
#ifdef SPMD
      work_2(i,j,1)=work_2(i,j,1)+viv_global(i,j,k)*ua(i,j,k)*dzp(k)*dyr_global(j)
#else
      work_2(i,j,1)=work_2(i,j,1)+viv(i,j,k)*ua(i,j,k)*dzp(k)*dyr(j)
#endif
      enddo
      enddo
      enddo
!
      do j=2,jmt_global
      do i=1,imt
      work_1(i,j,1)=work_1(i,j-1,1)+(work_2(i,j,1)+work_2(i,j-1,1))*0.5*1.0e-6
      enddo
      enddo
!
!$OMP PARALLEL DO PRIVATE (I)
      do i=1,imt
      work_1(i,1,1)=work_1(i,2,1)
      enddo
!
!$OMP PARALLEL DO PRIVATE (I,J)
      do j=1,jmt_global
      do i=1,imt
#ifdef SPMD
      if (viv_global(i,j,1)<0.5) then
#else
      if (viv(i,j,1)<0.5) then
#endif
      bsf(i,j)=spval
      else
      bsf(i,j)=work_1(i,j,1)
      endif
      enddo
      enddo
!
!     deallocate (work_1,work_2)
!
      return
      END SUBROUTINE BAROSF
!

!     =================
      subroutine diag_tracer(NNN) 
!     =================
!     written by liu hai long 2004, jan 
!
      IMPLICIT none
!
#include <netcdf.inc>
 
      integer, parameter :: xx1=121,xx2=281,yy1=61,yy2=170,zz1=1,zz2=6
!      integer, parameter :: xx1=1,xx2=imt,yy1=1,yy2=jmt_global,zz1=1,zz2=km
      CHARACTER ( LEN =   4 ) :: ftail
      CHARACTER ( LEN =  15 ) :: fname
      CHARACTER ( LEN =   8 ) :: dd
      CHARACTER ( LEN =   10 ) :: tt
      CHARACTER ( LEN =   5 ) :: zz
      INTEGER(r4)             :: vv(8)
      INTEGER :: nwmf,NNN

!     error status return
      integer:: iret

!     file id
      integer:: ncid

!     dimension id
      integer:: lon_dim,lat_dim,lev_dim,time_dim

!     dimension lenth
      integer, parameter :: lon_len=xx2-xx1+1,lat_len=yy2-yy1+1,lev_len=zz2-zz1+1,time_len=1


!     variable id
#ifdef ISO
      integer::  lat_id,lon_id,lev1_id,lev_id,time_id,net_id,trend_id,pen_id,&
               ax_id,ay_id,az_id,dx_id,dy_id,dz_id,&
               ax_iso_id,ay_iso_id,az_iso_id,dx_iso_id,dy_iso_id,dz_iso_id
#else
      integer::  lat_id,lon_id,lev1_id,lev_id,time_id,net_id,trend_id,pen_id,&
               ax_id,ay_id,az_id,dx_id,dy_id,dz_id
#endif

!     variable rank
#ifdef ISO
      integer, parameter ::lat_rank=1,lon_rank=1,lev_rank=1,time_rank=1,&
                           net_rank=3,trend_rank=4,pen_rank=4,&
                           ax_rank=4,ay_rank=4,az_rank=4,&
                           dx_rank=4,dy_rank=4,dz_rank=4,&
                       ax_iso_rank=4,ay_iso_rank=4,az_iso_rank=4,&
                       dx_iso_rank=4,dy_iso_rank=4,dz_iso_rank=4
#else
      integer, parameter ::lat_rank=1,lon_rank=1,lev_rank=1,time_rank=1,&
                           net_rank=3,trend_rank=4,pen_rank=4,&
                           ax_rank=4,ay_rank=4,az_rank=4,&
                           dx_rank=4,dy_rank=4,dz_rank=4
#endif
!
!     variable shapes
#ifdef ISO
      integer :: lat_dims(lat_rank),lon_dims(lon_rank),lev_dims(lev_rank),time_dims(time_rank), &
                  net_dims(net_rank),trend_dims(trend_rank),pen_dims(pen_rank),&
                  ax_dims(ax_rank),ay_dims(ay_rank),az_dims(az_rank),&
                  dx_dims(dx_rank),dy_dims(dy_rank),dz_dims(dz_rank),&
                  ax_iso_dims(ax_iso_rank),ay_iso_dims(ay_iso_rank),az_iso_dims(az_iso_rank),&
                  dx_iso_dims(dx_iso_rank),dy_iso_dims(dy_iso_rank),dz_iso_dims(dz_iso_rank)
#else
      integer :: lat_dims(lat_rank),lon_dims(lon_rank),lev_dims(lev_rank),time_dims(time_rank), &
                  net_dims(net_rank),trend_dims(trend_rank),pen_dims(pen_rank),&
                  ax_dims(ax_rank),ay_dims(ay_rank),az_dims(az_rank),&
                  dx_dims(dx_rank),dy_dims(dy_rank),dz_dims(dz_rank)
#endif

!     variables
      integer :: t0_cdf
      real :: t1x_cdf(lon_len) 
      real :: t1y_cdf(lat_len) 
      real :: t1z_cdf(lev_len) 
      real :: t2_cdf(lon_len,lat_len,time_len) 
      real :: t3_cdf(lon_len,lat_len,lev_len,time_len) 

!    start and count
      integer:: start1(1),count1(1)
      integer:: start2(2),count2(2)
      integer:: start3(3),count3(3)
      integer:: start4(4),count4(4)
!
      if (diag_budget.eqv..false.) return
 
!---------------------------------------------------------------------
!     output monthly results
!---------------------------------------------------------------------
!    file name
!
      IF (NNN==1) THEN
      fname(1:5)='MHEAT'
      ELSE
      fname(1:5)='MSALT'
      ENDIF
!
      nwmf=iyfm
      write (ftail,'(i4.4)') nwmf
      fname(6:9)=ftail
      fname(10:10)='-'
      write(fname(11:12),'(i2.2)')mon0
      fname(13:15)='.nc'
! 
!--------------------------------------------------------------
!     cdf output
!--------------------------------------------------------------
!     file defination
! enter define mode
         iret = nf_create (fname, NF_CLOBBER, ncid)
         CALL check_err (iret)
! define dimensions
         iret = nf_def_dim (ncid, 'lon', lon_len, lon_dim)
         CALL check_err (iret)
 
         iret = nf_def_dim (ncid, 'lat', lat_len, lat_dim)
         CALL check_err (iret)

         iret = nf_def_dim (ncid, 'lev', lev_len, lev_dim)
         CALL check_err (iret)
 
!         iret = nf_def_dim (ncid, 'time', NF_UNLIMITED, time_dim)
         iret = nf_def_dim (ncid, 'time', time_len, time_dim)
         CALL check_err (iret)
! define variables
         lon_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'lon', NF_REAL, lon_rank, lon_dims, lon_id)
         CALL check_err (iret)
!
         lat_dims (1) = lat_dim
         iret = nf_def_var (ncid, 'lat', NF_REAL, lat_rank, lat_dims, lat_id)
         CALL check_err (iret)
!
         lev_dims (1) = lev_dim
         iret = nf_def_var (ncid, 'lev', NF_REAL, lev_rank, lev_dims, lev_id)
         CALL check_err (iret)
!
         time_dims (1) = time_dim
         iret = nf_def_var (ncid, 'time', NF_INT, time_rank, time_dims, time_id)
         CALL check_err (iret)
 
         net_dims (3) = time_dim
         net_dims (2) = lat_dim
         net_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'net', NF_REAL, net_rank, net_dims, net_id)
         CALL check_err (iret)
 
         trend_dims (4) = time_dim
         trend_dims (3) = lev_dim
         trend_dims (2) = lat_dim
         trend_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'trend', NF_REAL, trend_rank, trend_dims, trend_id)
         CALL check_err (iret)
!
      IF (NNN==1) THEN
         pen_dims (4) = time_dim
         pen_dims (3) = lev_dim
         pen_dims (2) = lat_dim
         pen_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'pen', NF_REAL, pen_rank, pen_dims, pen_id)
         CALL check_err (iret)
      ENDIF
!
         ax_dims (4) = time_dim
         ax_dims (3) = lev_dim
         ax_dims (2) = lat_dim
         ax_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'ax', NF_REAL, ax_rank, ax_dims, ax_id)
         CALL check_err (iret)
!
         ay_dims (4) = time_dim
         ay_dims (3) = lev_dim
         ay_dims (2) = lat_dim
         ay_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'ay', NF_REAL, ay_rank, ay_dims, ay_id)
         CALL check_err (iret)
!
         az_dims (4) = time_dim
         az_dims (3) = lev_dim
         az_dims (2) = lat_dim
         az_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'az', NF_REAL, az_rank, az_dims, az_id)
         CALL check_err (iret)
!
         dx_dims (4) = time_dim
         dx_dims (3) = lev_dim
         dx_dims (2) = lat_dim
         dx_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'dx', NF_REAL, dx_rank, dx_dims, dx_id)
         CALL check_err (iret)
!
         dy_dims (4) = time_dim
         dy_dims (3) = lev_dim
         dy_dims (2) = lat_dim
         dy_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'dy', NF_REAL, dy_rank, dy_dims, dy_id)
         CALL check_err (iret)
!
         dz_dims (4) = time_dim
         dz_dims (3) = lev_dim
         dz_dims (2) = lat_dim
         dz_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'dz', NF_REAL, dz_rank, dz_dims, dz_id)
         CALL check_err (iret)
!
#ifdef ISO
         ax_iso_dims (4) = time_dim
         ax_iso_dims (3) = lev_dim
         ax_iso_dims (2) = lat_dim
         ax_iso_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'ax_iso', NF_REAL, ax_iso_rank, ax_iso_dims, ax_iso_id)
         CALL check_err (iret)
!
         ay_iso_dims (4) = time_dim
         ay_iso_dims (3) = lev_dim
         ay_iso_dims (2) = lat_dim
         ay_iso_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'ay_iso', NF_REAL, ay_iso_rank, ay_iso_dims, ay_iso_id)
         CALL check_err (iret)
!
         az_iso_dims (4) = time_dim
         az_iso_dims (3) = lev_dim
         az_iso_dims (2) = lat_dim
         az_iso_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'az_iso', NF_REAL, az_iso_rank, az_iso_dims, az_iso_id)
         CALL check_err (iret)
!
         dx_iso_dims (4) = time_dim
         dx_iso_dims (3) = lev_dim
         dx_iso_dims (2) = lat_dim
         dx_iso_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'dx_iso', NF_REAL, dx_iso_rank, dx_iso_dims, dx_iso_id)
         CALL check_err (iret)
!
         dy_iso_dims (4) = time_dim
         dy_iso_dims (3) = lev_dim
         dy_iso_dims (2) = lat_dim
         dy_iso_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'dy_iso', NF_REAL, dy_iso_rank, dy_iso_dims, dy_iso_id)
         CALL check_err (iret)
!
         dz_iso_dims (4) = time_dim
         dz_iso_dims (3) = lev_dim
         dz_iso_dims (2) = lat_dim
         dz_iso_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'dz_iso', NF_REAL, dz_iso_rank, dz_iso_dims, dz_iso_id)
         CALL check_err (iret)
#endif
!
! assign attributes
         iret = nf_put_att_text (ncid, lat_id, 'long_name', 21, 'latitude (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lat_id, 'units', 13, 'degrees_north')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lon_id, 'long_name', 22, 'longitude (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lon_id, 'units', 12, 'degrees_east')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lev_id, 'long_name', 18, 'depth (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lev_id, 'units', 5, 'meter')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, time_id, 'long_name', 4, 'time')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, time_id, 'units', 23, 'months since 1001-01-01')
         CALL check_err (iret)
 
         iret = nf_put_att_text (ncid, net_id, 'long_name', 16, 'net surface flux')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, net_id, 'units', 16, 'K/sec or psu/sec')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, net_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)
 
         iret = nf_put_att_text (ncid, trend_id, 'long_name', 10, 'time trend')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, trend_id, 'units', 16, 'K/sec or psu/sec')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, trend_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
!
      IF (NNN==1) THEN
         iret = nf_put_att_text (ncid, pen_id, 'long_name', 15, 'solar penetrate')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, pen_id, 'units', 13, 'degree/second')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, pen_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
      ENDIF
!
         iret = nf_put_att_text (ncid, ax_id, 'long_name', 11, 'x_advection')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, ax_id, 'units', 16, 'K/sec or psu/sec')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, ax_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
!
         iret = nf_put_att_text (ncid, ay_id, 'long_name', 11, 'y_advection')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, ay_id, 'units', 16, 'K/sec or psu/sec')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, ay_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
!
         iret = nf_put_att_text (ncid, az_id, 'long_name', 11, 'z_advection')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, az_id, 'units', 16, 'K/sec or psu/sec')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, az_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
!
         iret = nf_put_att_text (ncid, dx_id, 'long_name', 11, 'x_diffusion')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, dx_id, 'units', 16, 'K/sec or psu/sec')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, dx_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
!
         iret = nf_put_att_text (ncid, dy_id, 'long_name', 11, 'y_diffusion')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, dy_id, 'units', 16, 'K/sec or psu/sec')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, dy_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
!
         iret = nf_put_att_text (ncid, dz_id, 'long_name', 11, 'z_diffusion')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, dz_id, 'units', 16, 'K/sec or psu/sec')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, dz_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
!
#ifdef ISO
         iret = nf_put_att_text (ncid, ax_iso_id, 'long_name', 35, 'x_advection due to isopycnal mixing')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, ax_iso_id, 'units', 16, 'K/sec or psu/sec')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, ax_iso_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
!
         iret = nf_put_att_text (ncid, ay_iso_id, 'long_name', 35, 'y_advection due to isopycnal mixing')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, ay_iso_id, 'units', 16, 'K/sec or psu/sec')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, ay_iso_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
!
         iret = nf_put_att_text (ncid, az_iso_id, 'long_name', 35, 'z_advection due to isopycnal mixing')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, az_iso_id, 'units', 16, 'K/sec or psu/sec')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, az_iso_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
!
         iret = nf_put_att_text (ncid, dx_iso_id, 'long_name', 35, 'x_diffusion due to isopycnal mixing')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, dx_iso_id, 'units', 16, 'K/sec or psu/sec')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, dx_iso_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
!
         iret = nf_put_att_text (ncid, dy_iso_id, 'long_name', 35, 'y_diffusion due to isopycnal mixing')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, dy_iso_id, 'units', 16, 'K/sec or psu/sec')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, dy_iso_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
!
         iret = nf_put_att_text (ncid, dz_iso_id, 'long_name', 35, 'z_diffusion due to isopycnal mixing')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, dz_iso_id, 'units', 16, 'K/sec or psu/sec')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, dz_iso_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
#endif
!
!   define global attribute
         CALL date_and_time (dd,tt,zz,vv) 
 
         iret = NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'title', 26, 'monthly mean tracer budget')
         CALL check_err (iret)
         iret = NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'history', 20, tt //'  '//dd)
         CALL check_err (iret)
         iret = NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'source', 35, 'LASG/IAP Climate system Ocean Model')
         CALL check_err (iret)
! leave define mode
         iret = nf_enddef (ncid)
         CALL check_err (iret)
 
!----------------------------------------------------------
!     prepare data for storing
!----------------------------------------------------------
!
         DO i = xx1,xx2
         t1x_cdf(i-xx1+1)=lon(i)
         ENDDO
         iret = nf_put_var_real (ncid, lon_id, t1x_cdf)
         CALL check_err (iret)
!
         DO i = yy1,yy2
         t1y_cdf(i-yy1+1)=lat(i)
         ENDDO
         iret = nf_put_var_real (ncid, lat_id, t1y_cdf)
         CALL check_err (iret)
!
         DO k = zz1,zz2
         t1z_cdf(k-zz1+1)=lev(k)
         ENDDO
         iret = nf_put_var_real (ncid, lev_id, t1z_cdf)
         CALL check_err (iret)
!
         t0_cdf = month -1
         start1 (1)= 1
         count1 (1)= time_len
         iret = nf_put_vara_int (ncid, time_id,start1,count1,t0_cdf)
         CALL check_err (iret)
! 
! store variables
         start3 (1)= 1
         start3 (2)= 1
         start3 (3)= 1
         count3 (1)= lon_len
         count3 (2)= lat_len
         count3 (3)= time_len
 
!$OMP PARALLEL DO PRIVATE (J,I)
         DO j = yy1,yy2
            DO i = xx1,xx2
#ifdef SPMD
               IF (vit_global (i,j,1) > 0.5) THEN
                  t2_cdf (i-xx1+1,j-yy1+1,1)= netmon_io (i,j,NNN)/ (nmonth (mon0))
#else
               IF (vit (i,j,1) > 0.5) THEN
                  t2_cdf (i-xx1+1,j-yy1+1,1)= netmon (i,j,NNN)/ (nmonth (mon0))
#endif
               ELSE
                  t2_cdf (i-xx1+1,j-yy1+1,1)= spval
               END IF
            END DO
         END DO
 
         iret = nf_put_vara_real (ncid,net_id,start3, count3, t2_cdf)
         CALL check_err (iret)
! 
         start4 (1)= 1
         start4 (2)= 1
         start4 (3)= 1
         start4 (4)= 1
         count4 (1)= lon_len
         count4 (2)= lat_len
         count4 (3)= lev_len
         count4 (4)= time_len
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO k = zz1,zz2
            DO j = yy1,yy2
               DO i = xx1,xx2
#ifdef SPMD
                  IF (vit_global (i,j,k) > 0.5) THEN
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= trendmon_io (i,j,k,NNN)/ (nmonth (mon0))
#else
                  IF (vit (i,j,k) > 0.5) THEN
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= trendmon (i,j,k,NNN)/ (nmonth (mon0))
#endif
                  ELSE
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= spval
                  END IF
               END DO
            END DO
         END DO
 
         iret = nf_put_vara_real (ncid,trend_id, start4, count4, t3_cdf)
         CALL check_err (iret)
! 
      IF (NNN==1) THEN
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO k = zz1,zz2
            DO j = yy1,yy2
               DO i = xx1,xx2
#ifdef SPMD
                  IF (vit_global (i,j,k) > 0.5) THEN
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= penmon_io (i,j,k)/ (nmonth (mon0))
#else
                  IF (vit (i,j,k) > 0.5) THEN
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= penmon (i,j,k)/ (nmonth (mon0))
#endif
                  ELSE
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= spval
                  END IF
               END DO
            END DO
         END DO
 
         iret = nf_put_vara_real (ncid,pen_id, start4, count4, t3_cdf)
         CALL check_err (iret)
      ENDIF
!
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO k = zz1,zz2
            DO j = yy1,yy2
               DO i = xx1,xx2
#ifdef SPMD
                  IF (vit_global (i,j,k) > 0.5) THEN
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= axmon_io (i,j,k,NNN)/ (nmonth (mon0))
#else
                  IF (vit (i,j,k) > 0.5) THEN
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= axmon (i,j,k,NNN)/ (nmonth (mon0))
#endif
                  ELSE
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= spval
                  END IF
               END DO
            END DO
         END DO
 
         iret = nf_put_vara_real (ncid,ax_id, start4, count4, t3_cdf)
         CALL check_err (iret)
! 
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO k = zz1,zz2
            DO j = yy1,yy2
               DO i = xx1,xx2
#ifdef SPMD
                  IF (vit_global (i,j,k) > 0.5) THEN
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= aymon_io (i,j,k,NNN)/ (nmonth (mon0))
#else
                  IF (vit (i,j,k) > 0.5) THEN
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= aymon (i,j,k,NNN)/ (nmonth (mon0))
#endif
                  ELSE
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= spval
                  END IF
               END DO
            END DO
         END DO
 
         iret = nf_put_vara_real (ncid,ay_id, start4, count4, t3_cdf)
         CALL check_err (iret)
! 
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO k = zz1,zz2
            DO j = yy1,yy2
               DO i = xx1,xx2
#ifdef SPMD
                  IF (vit_global (i,j,k) > 0.5) THEN
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= azmon_io (i,j,k,NNN)/ (nmonth (mon0))
#else
                  IF (vit (i,j,k) > 0.5) THEN
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= azmon (i,j,k,NNN)/ (nmonth (mon0))
#endif
                  ELSE
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= spval
                  END IF
               END DO
            END DO
         END DO
 
         iret = nf_put_vara_real (ncid,az_id, start4, count4, t3_cdf)
         CALL check_err (iret)
! 
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO k = zz1,zz2
            DO j = yy1,yy2
               DO i = xx1,xx2
#ifdef SPMD
                  IF (vit_global (i,j,k) > 0.5) THEN
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= dxmon_io (i,j,k,NNN)/ (nmonth (mon0))
#else
                  IF (vit (i,j,k) > 0.5) THEN
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= dxmon (i,j,k,NNN)/ (nmonth (mon0))
#endif
                  ELSE
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= spval
                  END IF
               END DO
            END DO
         END DO
 
         iret = nf_put_vara_real (ncid,dx_id, start4, count4, t3_cdf)
         CALL check_err (iret)
! 
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO k = zz1,zz2
            DO j = yy1,yy2
               DO i = xx1,xx2
#ifdef SPMD
                  IF (vit_global (i,j,k) > 0.5) THEN
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= dymon_io (i,j,k,NNN)/ (nmonth (mon0))
#else
                  IF (vit (i,j,k) > 0.5) THEN
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= dymon (i,j,k,NNN)/ (nmonth (mon0))
#endif
                  ELSE
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= spval
                  END IF
               END DO
            END DO
         END DO
 
         iret = nf_put_vara_real (ncid,dy_id, start4, count4, t3_cdf)
         CALL check_err (iret)
! 
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO k = zz1,zz2
            DO j = yy1,yy2
               DO i = xx1,xx2
#ifdef SPMD
                  IF (vit_global (i,j,k) > 0.5) THEN
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= dzmon_io (i,j,k,NNN)/ (nmonth (mon0))
#else
                  IF (vit (i,j,k) > 0.5) THEN
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= dzmon (i,j,k,NNN)/ (nmonth (mon0))
#endif
                  ELSE
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= spval
                  END IF
               END DO
            END DO
         END DO
 
         iret = nf_put_vara_real (ncid,dz_id, start4, count4, t3_cdf)
         CALL check_err (iret)
! 
#ifdef ISO
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO k = zz1,zz2
            DO j = yy1,yy2
               DO i = xx1,xx2
#ifdef SPMD
                  IF (vit_global (i,j,k) > 0.5) THEN
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= axmon_iso_io (i,j,k,NNN)/ (nmonth (mon0))
#else
                  IF (vit (i,j,k) > 0.5) THEN
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= axmon_iso (i,j,k,NNN)/ (nmonth (mon0))
#endif
                  ELSE
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= spval
                  END IF
               END DO
            END DO
         END DO
 
         iret = nf_put_vara_real (ncid,ax_iso_id, start4, count4, t3_cdf)
         CALL check_err (iret)
! 
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO k = zz1,zz2
            DO j = yy1,yy2
               DO i = xx1,xx2
#ifdef SPMD
                  IF (vit_global (i,j,k) > 0.5) THEN
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= aymon_iso_io (i,j,k,NNN)/ (nmonth (mon0))
#else
                  IF (vit (i,j,k) > 0.5) THEN
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= aymon_iso (i,j,k,NNN)/ (nmonth (mon0))
#endif
                  ELSE
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= spval
                  END IF
               END DO
            END DO
         END DO
 
         iret = nf_put_vara_real (ncid,ay_iso_id, start4, count4, t3_cdf)
         CALL check_err (iret)
! 
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO k = zz1,zz2
            DO j = yy1,yy2
               DO i = xx1,xx2
#ifdef SPMD
                  IF (vit_global (i,j,k) > 0.5) THEN
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= azmon_iso_io (i,j,k,NNN)/ (nmonth (mon0))
#else
                  IF (vit (i,j,k) > 0.5) THEN
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= azmon_iso (i,j,k,NNN)/ (nmonth (mon0))
#endif
                  ELSE
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= spval
                  END IF
               END DO
            END DO
         END DO
 
         iret = nf_put_vara_real (ncid,az_iso_id, start4, count4, t3_cdf)
         CALL check_err (iret)
! 
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO k = zz1,zz2
            DO j = yy1,yy2
               DO i = xx1,xx2
#ifdef SPMD
                  IF (vit_global (i,j,k) > 0.5) THEN
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= dxmon_iso_io (i,j,k,NNN)/ (nmonth (mon0))
#else
                  IF (vit (i,j,k) > 0.5) THEN
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= dxmon_iso (i,j,k,NNN)/ (nmonth (mon0))
#endif
                  ELSE
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= spval
                  END IF
               END DO
            END DO
         END DO
 
         iret = nf_put_vara_real (ncid,dx_iso_id, start4, count4, t3_cdf)
         CALL check_err (iret)
! 
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO k = zz1,zz2
            DO j = yy1,yy2
               DO i = xx1,xx2
#ifdef SPMD
                  IF (vit_global (i,j,k) > 0.5) THEN
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= dymon_iso_io (i,j,k,NNN)/ (nmonth (mon0))
#else
                  IF (vit (i,j,k) > 0.5) THEN
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= dymon_iso (i,j,k,NNN)/ (nmonth (mon0))
#endif
                  ELSE
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= spval
                  END IF
               END DO
            END DO
         END DO
 
         iret = nf_put_vara_real (ncid,dy_iso_id, start4, count4, t3_cdf)
         CALL check_err (iret)
! 
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO k = zz1,zz2
            DO j = yy1,yy2
               DO i = xx1,xx2
#ifdef SPMD
                  IF (vit_global (i,j,k) > 0.5) THEN
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= dzmon_iso_io (i,j,k,NNN)/ (nmonth (mon0))
#else
                  IF (vit (i,j,k) > 0.5) THEN
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= dzmon_iso (i,j,k,NNN)/ (nmonth (mon0))
#endif
                  ELSE
                     t3_cdf (i-xx1+1,j-yy1+1,k-zz1+1,1)= spval
                  END IF
               END DO
            END DO
         END DO
 
         iret = nf_put_vara_real (ncid,dz_iso_id, start4, count4, t3_cdf)
         CALL check_err (iret)
#endif
! 
         iret = nf_CLOSE (ncid)
         CALL check_err (iret)
!
      RETURN
      END SUBROUTINE diag_tracer
!
!
!========================================
      SUBROUTINE diag_heat_transport(NNN)
!========================================
!
      implicit none
      integer :: NNN,NIN
      real,dimension(jmt_global,km) :: dpoux
      real,dimension(imt,jmt_global,km) :: ttpp
      real,dimension(imt,jmt_global) :: vvv
!     allocate (work_1(imt,jmt_global,1),&
!             work_2(imt,jmt_global,1),work_3(imt,jmt_global,1))
!
      if (diag_mth.eqv..false.) return
!
!$OMP PARALLEL DO PRIVATE (K,J)
      do k=1,km
      do j=1,jmt_global
#ifdef SPMD
      dpoux(j,k)=DZP(k)/OUX_global(J)
#else
      dpoux(j,k)=DZP(k)/OUX(J)
#endif
      enddo
      enddo

      ttpp=0.0

! vertical averaged

!$OMP PARALLEL DO PRIVATE (K,J,I)
      do k=1,km
      do j=2,jmm_global
      do i=2,imm
#ifdef SPMD
      if (NNN==1) then
      ttpp(i,j,k)=vsmon_io(i,j,k)/nmonth(mon0)*viv_global(i,j,k)*0.25*&
             (tsmon_io(i  ,j+1,k)+tsmon_io(i  ,j,k)+&
              tsmon_io(i-1,j+1,k)+tsmon_io(i-1,j,k))/nmonth(mon0)
      else
      ttpp(i,j,k)=vsmon_io(i,j,k)/nmonth(mon0)*viv_global(i,j,k)*0.25*&
              (ssmon_io(i  ,j+1,k)+ssmon_io(i  ,j,k)+&
               ssmon_io(i-1,j+1,k)+ssmon_io(i-1,j,k))/nmonth(mon0)
      endif
#else
      if (NNN==1) then
      ttpp(i,j,k)=vsmon(i,j,k)/nmonth(mon0)*viv(i,j,k)*0.25*&
              (tsmon(i  ,j+1,k)+tsmon(i  ,j,k)+&
               tsmon(i-1,j+1,k)+tsmon(i-1,j,k))/nmonth(mon0)
      else
      ttpp(i,j,k)=vsmon(i,j,k)/nmonth(mon0)*viv(i,j,k)*0.25*&
              (ssmon(i  ,j+1,k)+ssmon(i  ,j,k)+&
               ssmon(i-1,j+1,k)+ssmon(i-1,j,k))/nmonth(mon0)
      endif
#endif
      enddo
      enddo
      enddo
!
      do j=1,jmt_global
      do i=1,imt
      work_1(i,j,1)=0.0
      do k=1,km
      work_1(i,j,1)=work_1(i,j,1)+ttpp(i,j,k)*dpoux(j,k)
      enddo
      enddo
      enddo

      do j=1,jmt_global
      do i=1,imt
      work_2(i,j,1)=0.0
      do k=1,km
#ifdef ISO
#ifdef SPMD
      work_2(i,j,1)=work_2(i,j,1)+ddymon_iso_io(i,j,k,NNN)/nmonth(mon0)*dpoux(j,k)
#else
      work_2(i,j,1)=work_2(i,j,1)+ddymon_iso(i,j,k,NNN)/nmonth(mon0)*dpoux(j,k)
#endif
#else
#ifdef SPMD
      work_2(i,j,1)=work_2(i,j,1)+ddymon_io(i,j,k,NNN)/nmonth(mon0)*dpoux(j,k)
#else
      work_2(i,j,1)=work_2(i,j,1)+ddymon(i,j,k,NNN)/nmonth(mon0)*dpoux(j,k)
#endif
#endif
      enddo
      enddo
      enddo

#ifdef ISO
      do j=1,jmt_global
      do i=1,imt
      work_3(i,j,1)=0.0
      do k=1,km
#ifdef SPMD
      work_3(i,j,1)=work_3(i,j,1)+aaymon_iso_io(i,j,k,NNN)/nmonth(mon0)*dpoux(j,k)
#else
      work_3(i,j,1)=work_3(i,j,1)+aaymon_iso(i,j,k,NNN)/nmonth(mon0)*dpoux(j,k)
#endif
      enddo
      enddo
      enddo
#endif

      DO NIN=1,2

      IF(NIN==1) THEN
!$OMP PARALLEL DO PRIVATE (J,I)
      do j=1,jmt_global
      do i=1,imt
#ifdef SPMD
      vvv(i,j)=vit_global(i,j,1)
#else
      vvv(i,j)=vit(i,j,1)
#endif
      enddo
      enddo
!
      ELSE
!$OMP PARALLEL DO PRIVATE (K,J,I)
      do j=1,jmt_global
      do i=1,imt
#ifdef SPMD
         if (basin(i,j)==1.or.basin(i,j)==2) then
            vvv(i,j)=1.0*vit_global(i,j,1)
         else
            vvv(i,j)=0.0
         end if
#else
         if (basin(i,j)==1.or.basin(i,j)==2) then
            vvv(i,j)=1.0*vit(i,j,1)
         else
            vvv(i,j)=0.0
         end if
#endif
      enddo
      enddo
      ENDIF
!
      do j=1,jmm_global
      do i=2,imt
       ttpp(I,J,1)= vvv(I -1,J)*vvv(I -1,J +1)*vvv(I,J)*vvv(I,J +1)
      enddo
       ttpp(1,j,1)=ttpp(imm,j,1)
      enddo
      do i=1,imt
       ttpp(i,jmt_global,1)=0.0
      enddo
!
!$OMP PARALLEL DO PRIVATE (J,I)
      do j=1,jmt_global
      mth_adv(j,NIN,NNN)=0.0
      mth_dif(j,NIN,NNN)=0.0
#ifdef ISO
      mth_adv_iso(j,NIN,NNN)=0.0
#endif
      do i=2,imm
      mth_adv(j,NIN,NNN)=mth_adv(j,NIN,NNN)-work_1(i,j,1)*ttpp(i,j,1)
      mth_dif(j,NIN,NNN)=mth_dif(j,NIN,NNN)+work_2(i,j,1)*ttpp(i,j,1)
#ifdef ISO
      mth_adv_iso(j,NIN,NNN)=mth_adv_iso(j,NIN,NNN)-work_3(i,j,1)*ttpp(i,j,1)
#endif
      enddo
      enddo

      IF (NNN==1) then
!$OMP PARALLEL DO PRIVATE (J)
      do j=1,jmt_global
      mth_adv(j,NIN,NNN)=mth_adv(j,NIN,NNN)*D0*CP*1.0E-15
      mth_dif(j,NIN,NNN)=mth_dif(j,NIN,NNN)*D0*CP*1.0E-15
#ifdef ISO
      mth_adv_iso(j,NIN,NNN)=mth_adv_iso(j,NIN,NNN)*D0*CP*1.0E-15
#endif
      mth(j,NIN,NNN)=mth_adv(j,NIN,NNN)+mth_dif(j,NIN,NNN)+mth_adv_iso(j,NIN,NNN)
      enddo
      ELSE
!$OMP PARALLEL DO PRIVATE (J)
      do j=1,jmt_global
      mth_adv(j,NIN,NNN)=(mth_adv(j,NIN,NNN)*1000.+35)/35.
      mth_dif(j,NIN,NNN)=(mth_dif(j,NIN,NNN)*1000.+35)/35.
#ifdef ISO
      mth_adv_iso(j,NIN,NNN)=(mth_adv_iso(j,NIN,NNN)*1000.+35)/35.
#endif
      mth(j,NIN,NNN)=mth_adv(j,NIN,NNN)+mth_dif(j,NIN,NNN)+mth_adv_iso(j,NIN,NNN)
      enddo
      ENDIF
!
      ENDDO

!     deallocate (work_1,work_2,work_3)
      RETURN
      END SUBROUTINE diag_heat_transport

end module diag_mod
