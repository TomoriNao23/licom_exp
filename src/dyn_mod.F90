!  CVS: $Id: dyn_mod.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
module dyn_mod
#include <def-undef.h>
use precision_mod
use param_mod
!     ------------------------------------------------------------------
!     U V T S H0 W RHO
!     ------------------------------------------------------------------
      real(r8),dimension(imt,jmt)::ub,vb,h0,ubp,vbp,h0p
      real(r8),dimension(imt,jmt,km)::u,v,up,vp
      real(r8),dimension(imt,jmt,kmp1)::ws
      real(r8),dimension(imt,jmt)::h0l,h0f,h0bl,h0bf
      real(r8),dimension(imt,jmt,km)::utl,utf,vtl,vtf
!      real(r8),dimension(imt,jmt)::utt,vtt
!#ifdef SPMD
      real(r8),dimension(imt,jmt_global,km)::u_io,v_io
      real(r8),dimension(imt,jmt_global)::h0_io
#ifdef COUP
      real(r8),dimension(imt,jmt_global)::t_cpl_io,s_cpl_io,u_cpl_io,v_cpl_io,dhdx_io,dhdy_io,q_io
#endif
!#endif
!
!
!     ------------------------------------------------------------------
!     Pressure gradient
!     ------------------------------------------------------------------
      real(r8),dimension(:,:,:),allocatable::gg,dlu,dlv
      real(r8),dimension(:,:),allocatable::dlub,dlvb

#ifdef TIDE
!lhl20100801
!     ------------------------------------------------------------------
!     for tide
!     ------------------------------------------------------------------
      real(r8),dimension(imt,jmt) :: elevat,atm_elevat
!lhl20100801
#endif

end module dyn_mod

