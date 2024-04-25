!  CVS: $Id: global_to_local.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     ================
      subroutine global_to_local_1d(global,local)
!     ================
!     To transfer global 1-d data to local processor.
 
#include <def-undef.h>
#ifdef SPMD
use precision_mod
use param_mod
use pconst_mod
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
      IMPLICIT NONE
      real(r8) :: global(jmt_global),local(jmt)

!$OMP PARALLEL DO PRIVATE (J)
      do j=1,jmt
         local(j)=0.0
      end do

!$OMP PARALLEL DO PRIVATE (J)
      do j=1,jmt
         if (j_global(j)<=jmt_global) then
             local(j)=global(j_global(j))
         end if
      end do

#endif

      end subroutine global_to_local_1d
 
!--------------------------------------------------------------
      subroutine global_to_local_2d(global,local)
!     ================
!     To transfer global 2-d data to local processor.

#include <def-undef.h>
#ifdef SPMD
use precision_mod
use param_mod
use pconst_mod
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
      IMPLICIT NONE
      integer :: j_start,j_end,nth
      real(r8) :: global(imt,jmt_global),local(imt,jmt)

      if (mytid==0) then
         do nth=1,nproc-1
            call mpi_recv(j_start,1,mpi_integer,nth,tag_1d,mpi_comm_ocn,status,ierr)
!
!$OMP PARALLEL DO PRIVATE (J,I)
            do j=1,jmt
            do i=1,imt
               if ((j_start+j-1)<=jmt_global) then
                   local(i,j)=global(i,j_start+j-1)
               else
                   local(i,j)=0.0
               end if
            end do
            end do
            call mpi_send(local,imt*jmt,MPI_PR,nth,tag_1d,mpi_comm_ocn,ierr)
         end do
!
!$OMP PARALLEL DO PRIVATE (J,I)
         do j=1,jmt
         do i=1,imt
             local(i,j)=global(i,j_global(j))
         end do
         end do
      else
         do nth=1,nproc-1
            if (mytid==nth) then
               j_start=j_global(1)
               call mpi_send(j_start,1,mpi_integer,0,tag_1d,mpi_comm_ocn,ierr)
!
               call mpi_recv(local,imt*jmt,MPI_PR,0,tag_1d,mpi_comm_ocn,status,ierr)
            end if
         end do
      end if


!!$OMP PARALLEL DO PRIVATE (J,I)
!      do j=1,jmt
!      do i=1,imt
!         local(i,j)=0.0
!      end do
!      end do

!!$OMP PARALLEL DO PRIVATE (J,I)
!      do j=1,jmt
!         if (j_global(j)<=jmt_global) then
!             do i=1,imt
!             local(i,j)=global(i,j_global(j))
!             end do
!         end if
!      end do

#endif

      end subroutine global_to_local_2d

!--------------------------------------------------------------
      subroutine global_to_local_3d(global,local,kk)
!     ================
!     To transfer global 2-d data to local processor.

#include <def-undef.h>
#ifdef SPMD
use precision_mod
use param_mod
use pconst_mod
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
      IMPLICIT NONE
      integer :: kk
      real(r8) :: global(imt,jmt_global,kk),local(imt,jmt,kk)

!$OMP PARALLEL DO PRIVATE (K,I,J)
      do k=1,kk
      do j=1,jmt
      do i=1,imt
         local(i,j,k)=0.0
      end do
      end do
      end do

!$OMP PARALLEL DO PRIVATE (K,I,J)
      do k=1,kk
      do j=1,jmt
         if (j_global(j)<jmt_global) then
             do i=1,imt
                local(i,j,k)=global(i,j_global(j),k)
             end do
         end if
      end do
      end do

#endif

      end subroutine global_to_local_3d
 
!--------------------------------------------------------------
      subroutine global_to_local_4d(global,local,kk,mm)
!     ================
!     To transfer global 4-d data to local processor.

#include <def-undef.h>
#ifdef SPMD
use precision_mod
use param_mod
use pconst_mod
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
      IMPLICIT NONE
      integer :: kk,mm,j_start,j_end,nth
      real(r8) :: global(imt,jmt_global,kk,mm),local(imt,jmt,kk,mm)


      if (mytid==0) then
         do nth=1,nproc-1
            call mpi_recv(j_start,1,mpi_integer,nth,tag_1d,mpi_comm_ocn,status,ierr)
!
            do m=1,mm
            do k=1,kk
!$OMP PARALLEL DO PRIVATE (J,I)
            do j=1,jmt
            do i=1,imt
               if ((j_start+j-1)<=jmt_global) then
                   local(i,j,k,m)=global(i,j_start+j-1,k,m)
               else
                   local(i,j,k,m)=0.0
               end if
            end do
            end do
            end do
            end do
            call mpi_send(local,imt*jmt*kk*mm,MPI_PR,nth,tag_1d,mpi_comm_ocn,ierr)
         end do
!
         do m=1,mm
         do k=1,kk
!$OMP PARALLEL DO PRIVATE (J,I)
         do j=1,jmt
         do i=1,imt
             local(i,j,k,m)=global(i,j_global(j),k,m)
         end do
         end do
         end do
         end do
      else
         do nth=1,nproc-1
            if (mytid==nth) then
               j_start=j_global(1)
               call mpi_send(j_start,1,mpi_integer,0,tag_1d,mpi_comm_ocn,ierr)
!
               call mpi_recv(local,imt*jmt*kk*mm,MPI_PR,0,tag_1d,mpi_comm_ocn,status,ierr)
            end if
         end do
      end if

#endif

      end subroutine global_to_local_4d

!--------------------------------------------------------------
      subroutine global_to_local_4d_real(global,local,kk,mm)
!     ================
!     To transfer global 4-d data to local processor.

#include <def-undef.h>
#ifdef SPMD
use precision_mod
use param_mod
use pconst_mod
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
      IMPLICIT NONE
      integer :: kk,mm,j_start,j_end,nth
      real(r4) :: global(imt,jmt_global,kk,mm)
      real(r4) :: tmp(imt,jmt,kk,mm)
      real(r8) :: local(imt,jmt,kk,mm)


      if (mytid==0) then
         do nth=1,nproc-1
            call mpi_recv(j_start,1,mpi_integer,nth,tag_1d,mpi_comm_ocn,status,ierr)
!
            do m=1,mm
            do k=1,kk
!$OMP PARALLEL DO PRIVATE (J,I)
            do j=1,jmt
            do i=1,imt
               if ((j_start+j-1)<=jmt_global) then
                   tmp(i,j,k,m)=global(i,j_start+j-1,k,m)
               else
                   tmp(i,j,k,m)=0.0
               end if
            end do
            end do
            end do
            end do
            call mpi_send(tmp,imt*jmt*kk*mm,MPI_PR1,nth,tag_1d,mpi_comm_ocn,ierr)
         end do
!
         do m=1,mm
         do k=1,kk
!$OMP PARALLEL DO PRIVATE (J,I)
         do j=1,jmt
         do i=1,imt
             tmp(i,j,k,m)=global(i,j_global(j),k,m)
         end do
         end do
         end do
         end do
      else
         do nth=1,nproc-1
            if (mytid==nth) then
               j_start=j_global(1)
               call mpi_send(j_start,1,mpi_integer,0,tag_1d,mpi_comm_ocn,ierr)
!
               call mpi_recv(tmp,imt*jmt*kk*mm,MPI_PR1,0,tag_1d,mpi_comm_ocn,status,ierr)
            end if
         end do
      end if
!
      local=tmp
!
#endif

      end subroutine global_to_local_4d_real


