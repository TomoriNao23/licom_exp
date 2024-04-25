!  CVS: $Id: boundary.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      subroutine boundary
!     ================
!     To compute bounary of subdomain for the each processor.
 
#include <def-undef.h>
#ifdef SPMD
use param_mod
use pconst_mod
use msg_mod
      IMPLICIT NONE


      do j=jst,jmt
         if (mytid==0) then
            j_global(j)=jst_global+j-1
         else
            j_global(j)=jst_global-1+mytid*(jmt-num_overlap)+j
         end if
      end do
     
      if (j_global(jmt)<=jmt_global) then
         j_loop=jmt
      else
         j_loop=jmt-(j_global(jmt)-jmt_global)
      end if

      if (mytid==nproc-1) then
         if (j_global(jmt)<jmt_global.or.j_global(1)>jmt_global) then
           write(*,*) "ERROR in boundary! "
           write(*,*) "j_global(1),j_global(jmt),jmt_global=",j_global(1),j_global(jmt),jmt_global
           stop
         end if
      end if

      do k = 1,km
      do j = 1,jmm_global
         do i = 2,imt
            viv_global(i,j,k)= vit_global(i-1,j,k)*vit_global(i-1,j+1,k)* &
                         vit_global(i,j,k)* vit_global(i,j+1,k)
         end do
            viv_global(1,j,k)= viv_global (imm,j,k)
      end do
      do i = 1,imt
            viv_global(i,jmt_global,k)= 0.0
      end do
      end do


#endif

      end subroutine boundary
 
 
