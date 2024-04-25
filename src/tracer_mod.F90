!  CVS: $Id: tracer_mod.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
module tracer_mod
#include <def-undef.h>
!
use precision_mod
use param_mod
!     ------------------------------------------------------------------
!     U V T S H0 W RHO
!     ------------------------------------------------------------------
      real(r8),dimension(imt,jmt,km,NTRA)::at
      real(r8),dimension(imt,jmt,0:km,NTRA)::atb
      real(r8),dimension(imt,jmt,NTRA)::net
!
!lhl1204
!      real(r8),dimension(imt,jmt,km)::pdensity,pp,alpha,beta
!      real(r8),dimension(imt,jmt,km)::pdensityu,ppu,alphau,betau
!      real(r8),dimension(imt,jmt,km)::pdensity,alpha,beta
      real(r8),dimension(imt,jmt,km)::pdensity
      real(r8),dimension(imt,jmt)::amld
!lhl1204
!
      real(r8),dimension(imt,jmt,km,NTRA)::trend
      real(r8),dimension(imt,jmt,km,NTRA)::ax,ay,az
      real(r8),dimension(imt,jmt,km,NTRA)::dx,dy,dz
      real(r8),dimension(imt,jmt,km)::penetrate
!
      real(r8),dimension(imt,jmt,km,NTRA)::ddy
!
#ifdef ISO
      real(r8),dimension(imt,jmt,km,NTRA)::aay_iso,ddy_iso
      real(r8),dimension(imt,jmt,km,NTRA)::ax_iso,ay_iso,az_iso
      real(r8),dimension(imt,jmt,km,NTRA)::dx_iso,dy_iso,dz_iso
#endif
!
!#ifdef SPMD
      real(r8),dimension(imt,jmt_global,km,NTRA)::at_io
!#endif
!
!     ------------------------------------------------------------------
!     Sea Ice
!     ------------------------------------------------------------------
      real(r8),dimension(imt,jmt):: ITICE,ALEAD,TLEAD, HI
      real(r8),dimension(imt,jmt_global):: ITICE_io,ALEAD_io, HI_io
!#ifdef SPMD
!      real(r8),dimension(imt,jmt_global):: ITICE_io,ALEAD_io, HI_io
!#endif
!
#ifdef COUP
      real(r8),dimension(imt,jmt):: qice
#ifdef SPMD
      real(r8),dimension(imt,jmt_global):: qice_global
#endif
#endif
!
end module tracer_mod
