!  CVS: $Id: inirun.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      SUBROUTINE INIRUN
!     =================
!     INITIALIZING FOR ALL PHYSICAL FIELDS

#include <def-undef.h>
use param_mod
use pconst_mod
use dyn_mod
use tracer_mod
use forc_mod
use work_mod
#ifdef SPMD
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
#endif
#ifdef COUP
use shr_msg_mod
use shr_sys_mod
use buf_mod
use control_mod
#endif

      IMPLICIT NONE
#include <netcdf.inc>
!
  !----- local  ------
  integer            :: fid    ! nc domain file ID
  integer            :: dimid  ! nc dimension id
  integer            :: vid    ! nc variable ID
  integer            :: rcode  ! nc return code
  integer            :: ntim   ! temporary
!     Define Variables.
      integer*4   :: ncid, iret
      integer*4,  dimension(4) :: start(4)
      integer*4,  dimension(4) :: count(4)
      real(r8)  :: ddx !,ddy
      real(r8), allocatable :: yv(:,:,:)
      INTEGER :: NMFF,jg0,jg,jg1,jr

#ifdef COUP
!
     ncpl=1
!
     ! obtain grid variables
			nx = imt-num_overlap
			ny=  jmt-num_overlap
			nxg= imt - num_overlap
			nyg= jmt_global
     if (mytid ==0 ) then
        ny=jmt-num_overlap/2
				jj_start=1
		 		jj_end=ny
      else if (mytid == n_proc -1 ) then
        ny = j_loop-num_overlap/2
				jj_start=2
				jj_end=ny+1
			else
				jj_start=2
				jj_end=ny+1
     end if

!     allocate(xc(nx,ny))
!     allocate(yc(nx,ny))
!     allocate(yv(4,nx,ny))
!     allocate(mask(nx,ny))
!     allocate(area(nx,ny))
     allocate(xc(imt,jmt))
     allocate(yc(imt,jmt))
     allocate(yv(4,imt,jmt))
     allocate(mask(imt,jmt))
     allocate(area(imt,jmt))   
    
     xc = 0.0D0
     yc = 0.0D0
     yv = 0.0D0
     mask = 0.0D0
     area = 0.0D0  
     where(vit(:,:,1) > 0.0) mask = 1 
    do j=1,jmt !ny
    do i=1,imt !nx
       xc(i,j) = 0.0D0 + dlam*dble(i-1)
       yc(i,j) = 90.-wkj(j)/torad
    end do
    end do

        do j=1,jmt !ny
         do i=1,imt !nx
         	jg=j_global(j)
         	if(jg>=jmt_global)then
          	yv(1,i,j)=90.- WKJ_global(jmt_global)/torad-0.5
         	else
          	yv(1,i,j)=90.-(WKJ_global(jg)+WKJ_global(jg-(-1)))*0.5/torad
        	endif
          yv(2,i,j)=yv(1,i,j)
          yv(3,i,j)=90.-(WKJ_global(jg)+WKJ_global(jg+(-1)))*0.5/torad
          yv(4,i,j)=yv(3,i,j)
         enddo
         enddo

      do j =1,jmt !ny
       do i= 1,imt !nx
          area (i,j) = (sin(yv(3,i,j)*torad)-sin(yv(1,i,j)*torad))* torad !4.0* ATAN (1.0)/180.0 !
    end do
    end do


!    if (mytid == 0 ) then
!       do i= 1,nx
!         ddx = 1.0D0 - sinu(jmt-j)
!         ddy = dyt(jmt)/radius*0.5D0
!         area (i,ny) = ddx*ddy
!          area (i,ny) = (1.0D0 - cosu(jmt-ny))*dyt(jmt)/radius*0.5D0
!       end do
!    end if
!

  allocate (t_cpl(imt,jmt))
  allocate (s_cpl(imt,jmt))
  allocate (u_cpl(imt,jmt))
  allocate (v_cpl(imt,jmt))
  allocate (dhdx(imt,jmt))
  allocate (dhdy(imt,jmt))
  allocate (Q   (imt,jmt))

  allocate (taux (imt,jmt))
  allocate (tauy (imt,jmt))
  allocate (netsw(imt,jmt))
  allocate (lat1 (imt,jmt))
  allocate (sen  (imt,jmt))
  allocate (lwup (imt,jmt))
  allocate (lwdn (imt,jmt))
  allocate (melth(imt,jmt))
  allocate (salt (imt,jmt))
  allocate (prec (imt,jmt))
  allocate (evap (imt,jmt))
  allocate (meltw(imt,jmt))
  allocate (roff (imt,jmt))
  allocate (ifrac(imt,jmt))
  allocate (patm (imt,jmt))
  allocate (duu10n(imt,jmt))
#ifdef USE_OCN_CARBON
  allocate(co2_cpl(imt,jmt))
  allocate(pco2(imt,jmt))
#endif    
!
     write(6,*) 'Begin msg_pass(init)'
     call msg_pass('init')
     write(6,*) 'End msg_pass(init)'
    
    write(6,*) 'Begin mpi_barrier(init)'
     call mpi_barrier(mpi_comm_ocn,ierr)
    write(6,*) 'End mpi_barrier(init)'

#endif

!$OMP PARALLEL DO PRIVATE (J,I)
      DO J = 1,JMT
         DO I = 1,IMT
            UB (I,J)= 0.0D0
            VB (I,J)= 0.0D0
            H0 (I,J)= 0.0D0
            UBP (I,J)= 0.0D0
            VBP (I,J)= 0.0D0
            H0P (I,J)= 0.0D0
         END DO
      END DO


!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = 1,JMT
            DO I = 1,IMT
               U (I,J,K)= 0.0D0
               V (I,J,K)= 0.0D0
               UP (I,J,K)= 0.0D0
               VP (I,J,K)= 0.0D0
            END DO
         END DO
      END DO


!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KMP1
         DO J = 1,JMT
            DO I = 1,IMT
               WS (I,J,K)= 0.0D0
            END DO
         END DO
      END DO


!$OMP PARALLEL DO PRIVATE (J,I)
      DO J = 1,JMT
         DO I = 1,IMT
            H0L (I,J)= 0.0D0
            H0F (I,J)= 0.0D0
            H0BL (I,J)= 0.0D0
            H0BF (I,J)= 0.0D0
         END DO
      END DO

!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = 1,JMT
            DO I = 1,IMT
               UTL (I,J,K)= 0.0D0
               UTF (I,J,K)= 0.0D0
               VTL (I,J,K)= 0.0D0
               VTF (I,J,K)= 0.0D0
            END DO
         END DO
      END DO
!         DO J = 1,JMT
!            DO I = 1,IMT
!               UTT (I,J)= 0.0
!               VTT (I,J)= 0.0
!            END DO
!         END DO

!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,NTRA
      DO J = 1,JMT
         DO I = 1,IMT
            NET (I,J,K)= 0.0D0
         END DO
      END DO
      END DO

!$OMP PARALLEL DO PRIVATE (J,I)
      DO J = 1,JMT
         DO I = 1,IMT
            ITICE (I,J)= 0D0
            ALEAD (I,J)= 0.0D0
            TLEAD (I,J)= 0.0D0
            HI (I,J)= 0.0D0
         END DO
      END DO


!$OMP PARALLEL DO PRIVATE (J,I)
      DO J = 1,JMT
         DO I = 1,IMT
            PXB (I,J)= 0.0D0
            PYB (I,J)= 0.0D0
            PAX (I,J)= 0.0D0
            PAY (I,J)= 0.0D0
            WHX (I,J)= 0.0D0
            WHY (I,J)= 0.0D0
            WGP (I,J)= 0.0D0
         END DO
      END DO

!     ------------------------------------------------------------------
!     Output Arrays
!     ------------------------------------------------------------------
      CALL YY00
!

      MONTH = 1

      IF (NSTART == 1) THEN

!     ------------------------------------------------------------------
!     READ LEVITUS ANNUAL MEAN TEMPERATURE AND SALINITY
!     ------------------------------------------------------------------
#ifdef SPMD
         if (mytid==0) then
#if (!defined CDFIN)
         OPEN (81,FILE ='TSinitial',STATUS ='OLD',FORM ='UNFORMATTED')
         READ (81) AT_io
         CLOSE (81)
#else
!----------------------------------------------------
! Open netCDF file.
!----------------------------------------------------
      iret=nf_open('TSinitial',nf_nowrite,ncid)
      call check_err (iret)

!----------------------------------------------------
!   Retrieve data
!----------------------------------------------------
      start(1)=1 ; count(1)=imt
      start(2)=1 ; count(2)=jmt_global
      start(3)=1 ; count(3)=km
      start(4)=1 ; count(4)=1

      iret=nf_get_vara_double(ncid,   5,start,count, at_io(1,1,1,1))
      call check_err (iret)

      start(1)=1 ; count(1)=imt
      start(2)=1 ; count(2)=jmt_global
      start(3)=1 ; count(3)=km
      start(4)=1 ; count(4)=1
      iret=nf_get_vara_double(ncid,   6,start,count, at_io(1,1,1,2))
      call check_err (iret)
!
      iret = nf_close (ncid)
      call check_err (iret)



#endif
!
         end if
!Yu
         call global_to_local_4d(at_io,at,km,2)
!	write(*,'(i4,11f8.2)') mytid,((at(i,j,1,1),i=190,200),j=10,20)
#else
#if (!defined CDFIN)
         OPEN (81,FILE ='TSinitial',STATUS ='OLD',FORM ='UNFORMATTED')
         READ (81) AT_io
         CLOSE (81)
#else
!----------------------------------------------------
! Open netCDF file.
!----------------------------------------------------
      iret=nf_open('TSinitial',nf_nowrite,ncid)
      call check_err (iret)

!----------------------------------------------------
!   Retrieve data
!----------------------------------------------------
      start(1)=1 ; count(1)=imt
      start(2)=1 ; count(2)=jmt_global
      start(3)=1 ; count(3)=km
      start(4)=1 ; count(4)=1

      iret=nf_get_vara_double(ncid,   5,start,count, at_io(1,1,1,1))
      call check_err (iret)

      start(1)=1 ; count(1)=imt
      start(2)=1 ; count(2)=jmt_global
      start(3)=1 ; count(3)=km
      start(4)=1 ; count(4)=1
      iret=nf_get_vara_double(ncid,   6,start,count, at_io(1,1,1,2))
      call check_err (iret)

      iret = nf_close (ncid)
      call check_err (iret)


#endif
     at=at_io
#endif

!----------------------------------------------------
!   assign 0 to land grids of TSinital
!----------------------------------------------------
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO K = 1,KM
            DO J = 1,JMT
               DO I = 1,IMT
                  AT (I,J,K,1) = AT (I,J,K,1)*VIT(I,J,K)
                  AT (I,J,K,2) = (AT (I,J,K,2)- 35.0D0)*0.001D0*VIT(I,J,K)
               END DO
            END DO
         END DO
!
!!$OMP PARALLEL DO PRIVATE (K,J,I)
!         DO K = 1,KM
!            DO J = 1,JMT
!               DO I = 1,IMT
!                  AT (I,J,K,2) = (AT (I,J,K,2) - 35.0)*0.001
!               END DO
!            END DO
!         END DO

         DO N = 1,NTRA
!$OMP PARALLEL DO PRIVATE (K,J,I)
            DO K = 1,KM
               DO J = 1,JMT
                  DO I = 1,IMT
                     ATB (I,J,K,N) = AT (I,J,K,N)
#if (defined BOUNDARY)
                     RESTORE (I,J,K,N) = AT (I,J,K,N)
#endif
                  END DO
               END DO
            END DO


!nick
!$OMP PARALLEL DO PRIVATE (J,I)
            DO J = 1,JMT
               DO I = 1,IMT
                  ATB (I,J,0,N) = 0.0D0
               END DO
            END DO
!nick
         END DO
      ELSE

!     ------------------------------------------------------------------
!     READ INTERMEDIATE RESULTS (fort.22/fort.21)
!     ------------------------------------------------------------------

#if (defined BOUNDARY)
#ifdef SPMD
         if (mytid==0) then
#if (!defined CDFIN)
         OPEN (81,FILE ='TSinitial',STATUS ='OLD',FORM ='UNFORMATTED')
         READ (81) RESTORE_in_io
         CLOSE (81)
#else
!----------------------------------------------------
! Open netCDF file.
!----------------------------------------------------
      iret=nf_open('TSinitial',nf_nowrite,ncid)
      call check_err (iret)

!----------------------------------------------------
!   Retrieve data
!----------------------------------------------------
      start(1)=1 ; count(1)=imt
      start(2)=1 ; count(2)=jmt_global
      start(3)=1 ; count(3)=km
      start(4)=1 ; count(4)=1

      iret=nf_get_vara_double(ncid,   5,start,count, restore_io(1,1,1,1))
      call check_err (iret)

      start(1)=1 ; count(1)=imt
      start(2)=1 ; count(2)=jmt_global
      start(3)=1 ; count(3)=km
      start(4)=1 ; count(4)=1
      iret=nf_get_vara_double(ncid,   6,start,count, restore_io(1,1,1,2))
      call check_err (iret)
!
      iret = nf_close (ncid)
      call check_err (iret)
#endif
         end if
!Yu
      call global_to_local_4d(restore_io,restore,km,2)
#else
#if (!defined CDFIN)
         OPEN (81,FILE ='TSinitial',STATUS ='OLD',FORM ='UNFORMATTED')
         READ (81) RESTORE_io
         CLOSE (81)
#else
!----------------------------------------------------
! Open netCDF file.
!----------------------------------------------------
      iret=nf_open('TSinitial',nf_nowrite,ncid)
      call check_err (iret)

!----------------------------------------------------
!   Retrieve data
!----------------------------------------------------
      start(1)=1 ; count(1)=imt
      start(2)=1 ; count(2)=jmt_global
      start(3)=1 ; count(3)=km
      start(4)=1 ; count(4)=1

      iret=nf_get_vara_double(ncid,   5,start,count, restore_io(1,1,1,1))
      call check_err (iret)

      start(1)=1 ; count(1)=imt
      start(2)=1 ; count(2)=jmt_global
      start(3)=1 ; count(3)=km
      start(4)=1 ; count(4)=1
      iret=nf_get_vara_double(ncid,   6,start,count, restore_io(1,1,1,2))
      call check_err (iret)
!
      iret = nf_close (ncid)
      call check_err (iret)
#endif
      restore=restore_io
#endif

!----------------------------------------------------
!   assign 0 to land grids of TSinital
!----------------------------------------------------
!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO K = 1,KM
            DO J = 1,JMT
               DO I = 1,IMT
                  RESTORE (I,J,K,1) = RESTORE (I,J,K,1)*VIT(I,J,K)
                  RESTORE (I,J,K,2) = (RESTORE (I,J,K,2) - 35.0)*0.001*VIT(I,J,K)
               END DO
            END DO
         END DO
!
!!$OMP PARALLEL DO PRIVATE (K,J,I)
!         DO K = 1,KM
!            DO J = 1,JMT
!               DO I = 1,IMT
!                  RESTORE (I,J,K,2) = (RESTORE (I,J,K,2) - 35.0)*0.001
!               END DO
!            END DO
!         END DO
#endif

!
#ifdef SPMD
         if (mytid==0) then
         open(22,file='fort.22',form='unformatted')
         READ (22)H0_io,U_io,V_io,AT_io,HI_io,ITICE_io,ALEAD_io,MONTH
          !write(483,'(362f10.3)')(at_io(1,j,1,1),j=1,jmt_global)
          !write(454,'(362f10.3)')(u_io(1,j,1),j=1,jmt_global)
         end if
!Yu
      call mpi_bcast(month,1,mpi_integer,0,mpi_comm_ocn,ierr)
      call global_to_local_4d(h0_io,h0,1,1)
      call global_to_local_4d(hi_io,hi,1,1)
      call global_to_local_4d(itice_io,itice,1,1)
      call global_to_local_4d(alead_io,alead,1,1)
      call global_to_local_4d(u_io,u,km,1)
      call global_to_local_4d(v_io,v,km,1)
      call global_to_local_4d(at_io,at,km,2)

#ifdef COUP
         if (nstart==2) then
            month=(cdate/10000-1)*12+mod(cdate,10000)/100
            do j=1,jmt !ny
            do i=1,imt !nx
               t_cpl (i,j)  = 273.15D0+at(i,j,1,1)
               s_cpl (i,j)  = at(i,j,1,2)*1000.0D0+35.0D0
               q     (i,j)  = 0.0D0
               u_cpl (i,j)  = 0.0D0
               v_cpl (i,j)  = 0.0D0
               dhdx  (i,j)  = 0.0D0
               dhdy  (i,j)  = 0.0D0
            end do
            end do
         else
            if (mytid == 0) then
               read(22)t_cpl_io,s_cpl_io,u_cpl_io,v_cpl_io,dhdx_io,dhdy_io,q_io
               cdate  =  10000*((month-1)/12+1)+100*(mod(month-1,12)+1)+1
            end if
!
            call global_to_local_4d(t_cpl_io,work,1,1)
            do j=1,jmt !ny
            do i=1,imt !nx
               t_cpl(i,j)=work(i,j)
            end do
            end do
!
            call global_to_local_4d(s_cpl_io,work,1,1)
            do j=1,jmt !ny
            do i=1,imt !nx
               s_cpl(i,j)=work(i,j)
            end do
            end do
!
            call global_to_local_4d(u_cpl_io,work,1,1)
            do j=1,jmt !ny
            do i=1,imt !nx
               u_cpl(i,j)=work(i,j)
            end do
            end do
!
            call global_to_local_4d(v_cpl_io,work,1,1)
            do j=1,jmt !ny
            do i=1,imt !nx
               v_cpl(i,j)=work(i,j)
            end do
            end do
!
            call global_to_local_4d(dhdx_io,work,1,1)
            do j=1,jmt !ny
            do i=1,imt !nx
               dhdx(i,j)=work(i,j)
            end do
            end do
!
            call global_to_local_4d(dhdy_io,work,1,1)
            do j=1,jmt !ny
            do i=1,imt !nx
               dhdy(i,j)=work(i,j)
            end do
            end do
!
            call global_to_local_4d(q_io,work,1,1)
            do j=1,jmt !ny
            do i=1,imt !nx
               q(i,j)=work(i,j)
            end do
            end do
!
         end if
#endif
#else
         open(22,file='fort.22',form='unformatted')
         READ (22)H0_io,U_io,V_io,AT_io,HI_io,ITICE_io,ALEAD_io,MONTH
#ifdef COUP
         if (nstart==2) then
            month=(cdate/10000-1)*12+mod(cdate,10000)/100
            do j=1,jmt !ny
            do i=1,imt !nx
               t_cpl (i,j)  = 273.15D0+at(i,j,1,1)
               s_cpl (i,j)  = at(i,j,1,2)*1000.0D0+35.0D0
               q     (i,j)  = 0.0D0
               u_cpl (i,j)  = 0.0D0
               v_cpl (i,j)  = 0.0D0
               dhdx  (i,j)  = 0.0D0
               dhdy  (i,j)  = 0.0D0
            end do
            end do

         else
            read(22)t_cpl_io,s_cpl_io,u_cpl_io,v_cpl_io,dhdx_io,dhdy_io,q_io
            cdate  =  10000*((month-1)/12+1)+100*(mod(month-1,12)+1)+1
            do j=1,jmt !ny
            do i=1,imt !nx
               t_cpl (i,j)  = t_cpl_io(i,j)
               s_cpl (i,j)  = s_cpl_io(i,j)
               q     (i,j)  = q_io(i,j)
               u_cpl (i,j)  = u_cpl_io(i,j)
               v_cpl (i,j)  = v_cpl_io(i,j)
               dhdx  (i,j)  = dhdx(i,j)
               dhdy  (i,j)  = djdy(i,j)
            end do
            end do


         end if
#endif
         CLOSE(22)
         h0=H0_io
         U=U_io
         V=V_io
         AT=AT_io
         hi=HI_io
         ITICE=ITICE_io
         ALEAD=ALEAD_io
#endif


         NMFF = MOD (MONTH -1,12)

         CALL VINTEG (U,UB)
         CALL VINTEG (V,VB)

!$OMP PARALLEL DO PRIVATE (K,J,I)
         DO K = 1,KM
!YU         DO J = 2,JMM
            DO J = jst,jmt
               DO I = 1,IMT
                  UP (I,J,K) = U (I,J,K)
                  VP (I,J,K) = V (I,J,K)
                  UTF (I,J,K) = U (I,J,K)
                  VTF (I,J,K) = V (I,J,K)
                  ATB (I,J,K,1) = AT (I,J,K,1)
                  ATB (I,J,K,2) = AT (I,J,K,2)
               END DO
            END DO
         END DO
!            DO J = jst,jmt
!               DO I = 1,IMT
!                  UTT (I,J) = U (I,J,1)
!                  VTT (I,J) = V (I,J,1)
!               END DO
!            END DO


!$OMP PARALLEL DO PRIVATE (J,I)
!YU      DO J = 2,JMM
         DO J = jst,jmt
            DO I = 1,IMT
               H0P (I,J)= H0 (I,J)
               UBP (I,J)= UB (I,J)
               VBP (I,J)= VB (I,J)
               H0F (I,J)= H0 (I,J)
               H0BF (I,J)= H0 (I,J)
            END DO
         END DO
      END IF
!
      if (mytid==0)then
          write(6,*)"END-----------INIRUN !"
#ifdef COUP
          call shr_sys_flush(6)
#endif
      endif

      RETURN
      END SUBROUTINE INIRUN


