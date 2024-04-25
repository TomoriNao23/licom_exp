!  CVS: $Id: addps.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      SUBROUTINE ADDPS
!     ================
!     COMPENSATING THE LOSS OF GROSS MASS
 
#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use dyn_mod
#ifdef SPMD
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
#endif
      IMPLICIT NONE
!by linpf 20100309 
      real(r8), allocatable :: h0_en(:,:)     
!      real(r8),dimension(imt,jmt_global)::h0_en
!by linpf 20100309 
      REAL(r8)    :: ERROR,DH00,error0
      ERROR = 0.0D0
      
      
!by linpf 20100309       
!!$OMP PARALLEL DO PRIVATE (J,I),reduction(+:ERROR)      
#ifdef SPMD   
    allocate(h0_en(imt,jmt_global))
   call local_to_global_4d_2double(h0,h0_en,1,1)
         DO J = 1,jmt_global
         DO I = 1,imt
            ERROR = ERROR + DYT_global(J)* SINT_global (J)* H0_en(I,J)* vit_global(I,J,1)
         END DO
      END DO
#else
     
!!$OMP PARALLEL DO PRIVATE (J,I),reduction(+:ERROR)
      DO J = JSM,JMM
         DO I = 2,IMM
            ERROR = ERROR + DYT (J)* SINT (J)* H0 (I,J)* VIT (I,J,1)
         END DO
      END DO
      
#endif
!by linpf 20100309 
!      if(mytid==0)then
!       write(*,*)error
!      endif
      
#ifdef SPMD
!by linpf 20100309 
!       call mpi_reduce(error,error0,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr) 
!by linpf 20100309 
       error0=error
!       if(mytid==0)then
!        write(*,*)error0
!       endif
       call mpi_bcast(error0,1,MPI_PR,0,mpi_comm_ocn,ierr)
!        if(mytid==0)then
!        write(*,*)error0
!       endif
      dh00 = - error0/ asea
#else
      DH00 = - ERROR / ASEA
#endif
 
!$OMP PARALLEL DO PRIVATE (J,I)
      DO J = JST,JMT
         DO I = 1,IMT
            H0 (I,J)= (H0 (I,J) + DH00)* VIT (I,J,1)
         END DO
      END DO
 
      deallocate(h0_en)
 
      RETURN
      END SUBROUTINE ADDPS
 
 
