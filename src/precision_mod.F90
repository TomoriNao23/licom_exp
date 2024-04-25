module precision_mod
#include <def-undef.h>
!
#ifdef SPMD
#include <mpif.h>
#endif
!-------------------------------------------------------------------------------
#ifdef D_PRECISION
integer, parameter :: r8 = selected_real_kind(12)
integer, parameter :: r4 = selected_real_kind(4)
integer, parameter :: i8 = selected_int_kind(13)

#ifdef SPMD
integer*4, parameter :: MPI_PR = MPI_DOUBLE_PRECISION
integer*4, parameter :: MPI_PR1 = MPI_REAL
#endif

#else
integer, parameter :: r8 = selected_real_kind(4)
integer, parameter :: r4 = selected_real_kind(4)
integer, parameter :: i8 = selected_int_kind(8)

#ifdef SPMD
integer*4, parameter :: MPI_PR = MPI_REAL
integer*4, parameter :: MPI_PR1 = MPI_REAL
#endif

#endif
end module precision_mod

