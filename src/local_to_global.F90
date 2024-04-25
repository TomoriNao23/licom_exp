!  CVS: $Id: local_to_global.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     ================
      subroutine local_to_global(hist_output,rest_output)
!     ================
!     To transfer global 1-d data to local processor.
 
#include <def-undef.h>
#ifdef SPMD
use precision_mod
use param_mod
use pconst_mod
use output_mod
use dyn_mod
#ifdef COUP
use buf_mod
#endif
use tracer_mod
use work_mod
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn

      IMPLICIT NONE
      integer :: itag,j_start
      logical :: hist_output,rest_output 

     if (hist_output) then
         call local_to_global_4d(z0mon,z0mon_io,1,1)
         call local_to_global_4d(himon,himon_io,1,1)
         call local_to_global_4d(hdmon,hdmon_io,1,1)
         call local_to_global_4d(netmon,netmon_io,2,1)
         call local_to_global_4d(tsmon,tsmon_io,km,1)
         call local_to_global_4d(ssmon,ssmon_io,km,1)
         call local_to_global_4d(usmon,usmon_io,km,1)
         call local_to_global_4d(vsmon,vsmon_io,km,1)
         call local_to_global_4d(wsmon,wsmon_io,km,1)
         call local_to_global_4d(icmon,icmon_io,2,1)
#if (defined SMAG_OUT)
         call local_to_global_4d(am3mon,am3mon_io,km,1)
#endif
!
         call local_to_global_4d(axmon,axmon_io,km,2)
         call local_to_global_4d(aymon,aymon_io,km,2)
         call local_to_global_4d(azmon,azmon_io,km,2)
         call local_to_global_4d(dxmon,dxmon_io,km,2)
         call local_to_global_4d(dymon,dymon_io,km,2)
         call local_to_global_4d(dzmon,dzmon_io,km,2)
         call local_to_global_4d(ddymon,ddymon_io,km,2)
#ifdef ISO
         call local_to_global_4d(axmon_iso,axmon_iso_io,km,2)
         call local_to_global_4d(aymon_iso,aymon_iso_io,km,2)
         call local_to_global_4d(azmon_iso,azmon_iso_io,km,2)
         call local_to_global_4d(dxmon_iso,dxmon_iso_io,km,2)
         call local_to_global_4d(dymon_iso,dymon_iso_io,km,2)
         call local_to_global_4d(dzmon_iso,dzmon_iso_io,km,2)
!
         call local_to_global_4d(aaymon_iso,aaymon_iso_io,km,2)
         call local_to_global_4d(ddymon_iso,ddymon_iso_io,km,2)
#endif
         call local_to_global_4d(trendmon,trendmon_io,km,2)
         call local_to_global_4d(penmon,penmon_io,km,1)
         call local_to_global_4d(mldmon,mldmon_io,1,1)
         call local_to_global_4d(akmmon,akmmon_io,km,1)
         call local_to_global_4d(aktmon,aktmon_io,km,1)
         call local_to_global_4d(aksmon,aksmon_io,km,1)
!
!linpf091126
!         
         call local_to_global_4d(sumon,sumon_io,1,1)
         call local_to_global_4d(svmon,svmon_io,1,1)
         call local_to_global_4d(lthfmon,lthfmon_io,1,1)
         call local_to_global_4d(sshfmon,sshfmon_io,1,1)
         call local_to_global_4d(lwvmon,lwvmon_io,1,1)
         call local_to_global_4d(swvmon,swvmon_io,1,1)
!
!linpf091126
!
     end if
!
     if (rest_output) then
         call local_to_global_4d_2double(h0,h0_io,1,1)
         call local_to_global_4d_2double(hi,hi_io,1,1)
         call local_to_global_4d_2double(itice,itice_io,1,1)
         call local_to_global_4d_2double(alead,alead_io,1,1)
         call local_to_global_4d_2double(u,u_io,km,1)
         call local_to_global_4d_2double(v,v_io,km,1)
         call local_to_global_4d_2double(at,at_io,km,2)
#ifdef COUP
!
!         do j=2,jmt-1
!         do i=1,imt
!            work(i,j)=t_cpl(i,jmt-j+1)
!         end do
!         end do
          work=t_cpl
         call local_to_global_4d_2double(work,t_cpl_io,1,1)
!
!         do j=2,jmt-1
!         do i=1,imt
!            work(i,j)=s_cpl(i,jmt-j+1)
!         end do
!         end do
          work=s_cpl
         call local_to_global_4d_2double(work,s_cpl_io,1,1)
!
!         do j=2,jmt-1
!         do i=1,imt
!            work(i,j)=u_cpl(i,jmt-j+1)
!         end do
!         end do
          work=u_cpl
         call local_to_global_4d_2double(work,u_cpl_io,1,1)
!
!         do j=2,jmt-1
!         do i=1,imt
!            work(i,j)=v_cpl(i,jmt-j+1)
!         end do
!         end do
         work=v_cpl
         call local_to_global_4d_2double(work,v_cpl_io,1,1)
!
!         do j=2,jmt-1
!         do i=1,imt
!            work(i,j)=dhdx(i,jmt-j+1)
!         end do
!         end do
          work=dhdx
         call local_to_global_4d_2double(work,dhdx_io,1,1)
!
!         do j=2,jmt-1
!         do i=1,imt
!            work(i,j)=dhdy(i,jmt-j+1)
!         end do
!         end do
           work=dhdy
         call local_to_global_4d_2double(work,dhdy_io,1,1)
!
!         do j=2,jmt-1
!         do i=1,imt
!            work(i,j)=q(i,jmt-j+1)
!         end do
!         end do
          work=q
         call local_to_global_4d_2double(work,q_io,1,1)
!
#endif
     end if

!           end do
#endif
      return
      end subroutine local_to_global
! 
     subroutine local_to_global_4d(local,global,kk,mm)
!     ================
!    To transfer global 1-d data to local processor.

#include <def-undef.h>
#ifdef SPMD
use precision_mod
use param_mod
use pconst_mod
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn

      IMPLICIT NONE
      integer :: itag,j_start,kk,mm,nth
      real(r4)    :: local(imt,jmt,kk,mm),global(imt,jmt_global,kk,mm),work(imt,jmt,kk,mm)
!
      if (mytid == 0) then
 
!$OMP PARALLEL DO PRIVATE (I,J,M,K)
           do m=1,mm
           do k=1,kk
           do j=1,jmt
           do i=1,imt
              global(i,j_global(j),k,m)=local(i,j,k,m)
           end do
           end do
           end do
           end do
!
        do nth=1,nproc-1
!
           itag=nth+800
           call mpi_recv(j_start,1,mpi_integer,nth,itag,mpi_comm_ocn,status,ierr)
           itag=nth+801
           call mpi_recv(work,imt*jmt*kk*mm,mpi_real,nth,itag,mpi_comm_ocn,status,ierr)
!$OMP PARALLEL DO PRIVATE (I,J)
           do m=1,mm
           do k=1,kk
           do j=1,jmt
           do i=1,imt
              if ((j_start+j)<=jmt_global) then
                 global(i,j+j_start,k,m)=work(i,j,k,m)
              end if
           end do
           end do
           end do
           end do
         end do
      else
         j_start=j_global(1)-1
         itag=mytid+800
         call mpi_send(j_start,1,mpi_integer,0,itag,mpi_comm_ocn,ierr)
         itag=mytid+801
         call mpi_send(local,imt*jmt*kk*mm,mpi_real,0,itag,mpi_comm_ocn,ierr)
      end if
!
#endif
     return
     end subroutine local_to_global_4d
     subroutine local_to_global_4d_double(local,global,kk,mm)
!     ================
!    To transfer global 1-d data to local processor.

#include <def-undef.h>
#ifdef SPMD
use precision_mod
use param_mod
use pconst_mod
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn

      IMPLICIT NONE
      integer :: itag,j_start,kk,mm,nth
      real(r8):: local(imt,jmt,kk,mm),work(imt,jmt,kk,mm)
      real(r4):: global(imt,jmt_global,kk,mm)
!
      if (mytid == 0) then
 
!$OMP PARALLEL DO PRIVATE (I,J,M,K)
           do m=1,mm
           do k=1,kk
           do j=1,jmt
           do i=1,imt
              global(i,j_global(j),k,m)=local(i,j,k,m)
           end do
           end do
           end do
           end do
!
        do nth=1,nproc-1
!
           itag=nth+800
           call mpi_recv(j_start,1,mpi_integer,nth,itag,mpi_comm_ocn,status,ierr)
           itag=nth+801
           call mpi_recv(work,imt*jmt*kk*mm,MPI_PR,nth,itag,mpi_comm_ocn,status,ierr)
!$OMP PARALLEL DO PRIVATE (I,J)
           do m=1,mm
           do k=1,kk
           do j=1,jmt
           do i=1,imt
              if ((j_start+j)<=jmt_global) then
                 global(i,j+j_start,k,m)=work(i,j,k,m)
              end if
           end do
           end do
           end do
           end do
         end do
      else
         j_start=j_global(1)-1
         itag=mytid+800
         call mpi_send(j_start,1,mpi_integer,0,itag,mpi_comm_ocn,ierr)
         itag=mytid+801
         call mpi_send(local,imt*jmt*kk*mm,MPI_PR,0,itag,mpi_comm_ocn,ierr)
      end if
!
#endif
     return
     end subroutine local_to_global_4d_double

!add by linpf
    subroutine local_to_global_4d_2double(local,global,kk,mm)
!     ================
!    To transfer global 1-d data to local processor.

#include <def-undef.h>
#ifdef SPMD
use precision_mod
use param_mod
use pconst_mod
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn

      IMPLICIT NONE
      integer :: itag,j_start,kk,mm,nth
      real(r8):: local(imt,jmt,kk,mm),work(imt,jmt,kk,mm)
      real(r8):: global(imt,jmt_global,kk,mm)
!
      if (mytid == 0) then

!$OMP PARALLEL DO PRIVATE (I,J,M,K)
           do m=1,mm
           do k=1,kk
           do j=1,jmt
           do i=1,imt
              global(i,j_global(j),k,m)=local(i,j,k,m)
           end do
           end do
           end do
           end do


        do nth=1,nproc-1
!
           itag=nth+800
           call mpi_recv(j_start,1,mpi_integer,nth,itag,mpi_comm_ocn,status,ierr)
           itag=nth+801
           call mpi_recv(work,imt*jmt*kk*mm,MPI_PR,nth,itag,mpi_comm_ocn,status,ierr)
!$OMP PARALLEL DO PRIVATE (I,J)
           do m=1,mm
           do k=1,kk
           do j=1,jmt
           do i=1,imt
              if ((j_start+j)<=jmt_global) then
                 global(i,j+j_start,k,m)=work(i,j,k,m)
              end if
           end do
           end do
           end do
           end do
         end do
      else
         j_start=j_global(1)-1
         itag=mytid+800
         call mpi_send(j_start,1,mpi_integer,0,itag,mpi_comm_ocn,ierr)
         itag=mytid+801
         call mpi_send(local,imt*jmt*kk*mm,MPI_PR,0,itag,mpi_comm_ocn,ierr)
      end if
!
#endif
     return
     end subroutine local_to_global_4d_2double
!add by lych
subroutine local_to_global_4d_r8(local,global,kk,mm)
!     ================
!    To transfer global 1-d data to local processor.

#include <def-undef.h>
#ifdef SPMD
use precision_mod
use param_mod
use pconst_mod
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn

      IMPLICIT NONE
      integer :: itag,j_start,kk,mm,nth
      real(r8):: local(imt,jmt,kk,mm),work(imt,jmt,kk,mm)
      real(r8):: global(imt,jmt_global,kk,mm)
!
      if (mytid == 0) then

!$OMP PARALLEL DO PRIVATE (I,J,M,K)
           do m=1,mm
           do k=1,kk
           do j=1,jmt
           do i=1,imt
              global(i,j_global(j),k,m)=local(i,j,k,m)
           end do
           end do
           end do
           end do
!
        do nth=1,nproc-1
!
           itag=nth+800
           call mpi_recv(j_start,1,mpi_integer,nth,itag,mpi_comm_ocn,status,ierr)
           itag=nth+801
           call mpi_recv(work,imt*jmt*kk*mm,MPI_PR,nth,itag,mpi_comm_ocn,status,ierr)
!$OMP PARALLEL DO PRIVATE (I,J)
           do m=1,mm
           do k=1,kk
           do j=1,jmt
           do i=1,imt
              if ((j_start+j)<=jmt_global) then
                 global(i,j+j_start,k,m)=work(i,j,k,m)
              end if
           end do
           end do
           end do
           end do
end do
      else
         j_start=j_global(1)-1
         itag=mytid+800
         call mpi_send(j_start,1,mpi_integer,0,itag,mpi_comm_ocn,ierr)
         itag=mytid+801
         call mpi_send(local,imt*jmt*kk*mm,MPI_PR,0,itag,mpi_comm_ocn,ierr)
      end if
!
#endif
     return
     end subroutine local_to_global_4d_r8

