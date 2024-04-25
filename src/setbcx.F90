!  CVS: $Id: setbcx.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
#include <def-undef.h>
 
#if (defined ISO)
!     ===================================
      SUBROUTINE setbcx (a, imt, jmtorkm)
!     ===================================
use precision_mod
      IMPLICIT NONE
      INTEGER :: imt,jmtorkm,k
      REAL(r8) :: a(imt,jmtorkm)
      DO k = 1,jmtorkm
         a (1,k) = a (imt -1,k)
 
         a (imt,k) = a (2,k)
      END DO
 
      RETURN
      END SUBROUTINE setbcx
 
#else
      SUBROUTINE setbcx ()
      RETURN
      END SUBROUTINE setbcx
#endif 
 
