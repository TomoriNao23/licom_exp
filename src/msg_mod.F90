!  CVS: $Id: msg_mod.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
module msg_mod
#include <def-undef.h>
#if ( defined SPMD ) || (defined COUP)
#include <mpif.h>
!-------------------------------------------------------------------------------
!
! Author: Yongqiang YU  ( 12 Nov, 2002)
!
!-------------------------------------------------------------------------------
      integer, parameter :: tag_1d=10, tag_2d=100, tag_3d=200,tag_4d=300
      integer:: nproc
      integer :: status (MPI_STATUS_SIZE) ! Status of message
!      integer ,save            :: mpi_comm_ocn !by linpf 20100316
      integer :: mpi_comm_ocn
#endif
end module msg_mod
