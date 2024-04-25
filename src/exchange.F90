
!  CVS: $Id: exchange.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     ================
      subroutine exchange_2d(aa)
!     ================
!     To compute bounary of subdomain for the each processor.

#include <def-undef.h>
#ifdef SPMD
use precision_mod
use param_mod
use pconst_mod
use msg_mod
      IMPLICIT NONE

      real(r8) :: aa(imt,jmt),bb(imt),cc(imt)
!
!$OMP PARALLEL DO PRIVATE (I)
     do i=1,imt
        bb(i)=aa(i,2)
        cc(i)=aa(i,jem)
     end do

      if (mytid< nproc-1) then
         call mpi_send(cc,imt,MPI_PR,mytid+1,tag_2d,mpi_comm_ocn,ierr)
      end if

      if (mytid>0) then
         call mpi_recv(aa(1,1),imt,MPI_PR,mytid-1,tag_2d,mpi_comm_ocn,status,ierr)
      end if

      if (mytid> 0) then
         call mpi_send(bb,imt,MPI_PR,mytid-1,tag_2d,mpi_comm_ocn,ierr)
      end if

      if (mytid< nproc-1) then
         call mpi_recv(aa(1,jmt),imt,MPI_PR,mytid+1,tag_2d,mpi_comm_ocn,status,ierr)
      end if

#endif
     end subroutine exchange_2d
!
!
      subroutine exchange_3d(aa,kk)
!     ================
!     To compute bounary of subdomain for the each processor.

#include <def-undef.h>
#ifdef SPMD
use precision_mod
use param_mod
use pconst_mod
use msg_mod
      IMPLICIT NONE

      integer :: kk
      real(r8)    :: aa(imt,jmt,kk),bb(imt,kk),cc(imt,kk),dd(imt,kk),ee(imt,kk)


!$OMP PARALLEL DO PRIVATE (K,I)
         do k=1,kk
         do i=1,imt
            bb(i,k)=aa(i,2,k)
            cc(i,k)=aa(i,jem,k)
         end do
         end do


      if (mytid< nproc-1) then
         call mpi_send(cc,imt*kk,MPI_PR,mytid+1,tag_3d,mpi_comm_ocn,ierr)
      end if

      if (mytid>0) then
         call mpi_recv(dd,imt*kk,MPI_PR,mytid-1,tag_3d,mpi_comm_ocn,status,ierr)
      end if

      if (mytid> 0) then
         call mpi_send(bb,imt*kk,MPI_PR,mytid-1,tag_3d,mpi_comm_ocn,ierr)
      end if

      if (mytid< nproc-1) then
         call mpi_recv(ee,imt*kk,MPI_PR,mytid+1,tag_3d,mpi_comm_ocn,status,ierr)
      end if



     if (mytid==0) then
!$OMP PARALLEL DO PRIVATE (K,I)
         do k=1,kk
         do i=1,imt
            aa(i,jmt,k)=ee(i,k)
         end do
         end do
     else if (mytid==nproc-1) then
!$OMP PARALLEL DO PRIVATE (K,I)
         do k=1,kk
         do i=1,imt
            aa(i,1,k)=dd(i,k)
         end do
         end do
     else
!$OMP PARALLEL DO PRIVATE (K,I)
         do k=1,kk
         do i=1,imt
            aa(i,1,k)=dd(i,k)
            aa(i,jmt,k)=ee(i,k)
         end do
         end do
      end if

#endif
      return
      end subroutine exchange_3d


      subroutine exchange_3d_iso(aa,kk)
!     ================
!     To compute bounary of subdomain for the each processor.

#include <def-undef.h>
#ifdef SPMD
use precision_mod
use param_mod
use pconst_mod
use msg_mod
      IMPLICIT NONE

      integer :: kk
      real(r8)    :: aa(imt,kk,jmt),bb(imt,kk),cc(imt,kk)


!$OMP PARALLEL DO PRIVATE (K,I)
         do k=1,kk
         do i=1,imt
            bb(i,k)=aa(i,k,2)
            cc(i,k)=aa(i,k,jem)
         end do
         end do

      if (mytid< nproc-1) then
         call mpi_send(cc,imt*kk,MPI_PR,mytid+1,tag_4d,mpi_comm_ocn,ierr)
      end if
      if (mytid > 0 ) then
         call mpi_recv(aa(1,1,1),imt*kk,MPI_PR,mytid-1,tag_4d,mpi_comm_ocn,status,ierr)
      end if


      if (mytid > 0 ) then
         call mpi_send(bb,imt*kk,MPI_PR,mytid-1,tag_4d,mpi_comm_ocn,ierr)
      end if
      if (mytid< nproc-1) then
         call mpi_recv(aa(1,1,jmt),imt*kk,MPI_PR,mytid+1,tag_4d,mpi_comm_ocn,status,ierr)
      end if

#endif
      return
      end subroutine exchange_3d_iso


      subroutine exchange_pack(aa,bb,cc)
!     ================
!     To compute bounary of subdomain for the each processor.

#include <def-undef.h>
#ifdef SPMD
use precision_mod
use param_mod
use pconst_mod
use msg_mod
      IMPLICIT NONE

      real(r8):: aa(imt,jmt),bb(imt,jmt),cc(imt,jmt)
      real(r8):: outbuf_1(imt*3),outbuf_2(imt*3)
      integer:: bsize,pos
      bsize=imt*3
!
!$omp parallel do private(i)
      do i=1,imt
         outbuf_2(i)=aa(i,jem)
         outbuf_2(imt+i)=bb(i,jem)
         outbuf_2(2*imt+i)=cc(i,jem)
      end do

      if (mytid< nproc-1) then
         call mpi_send(outbuf_2,imt*3,MPI_PR,mytid+1,tag_2d,mpi_comm_ocn,ierr)
      end if

      if (mytid>0) then
         call mpi_recv(outbuf_1,imt*3,MPI_PR,mytid-1,tag_2d,mpi_comm_ocn,status,ierr)
!$omp parallel do private(i)
        do i=1,imt
           aa(i,1)=outbuf_1(i)
           bb(i,1)=outbuf_1(imt+i)
           cc(i,1)=outbuf_1(2*imt+i)
        end do
      end if

!$omp parallel do private(i)
      do i=1,imt
         outbuf_1(i)=aa(i,2)
         outbuf_1(imt+i)=bb(i,2)
         outbuf_1(2*imt+i)=cc(i,2)
      end do

      if (mytid> 0) then
         call mpi_send(outbuf_1,imt*3,MPI_PR,mytid-1,tag_2d,mpi_comm_ocn,ierr)
      end if

      if (mytid< nproc-1) then
         call mpi_recv(outbuf_2,imt*3,MPI_PR,mytid+1,tag_2d,mpi_comm_ocn,status,ierr)
!$omp parallel do private(i)
        do i=1,imt
           aa(i,jmt)=outbuf_2(i)
           bb(i,jmt)=outbuf_2(imt+i)
           cc(i,jmt)=outbuf_2(2*imt+i)
        end do
      end if

#endif
     end subroutine exchange_pack

