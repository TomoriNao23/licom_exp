#define LOGMSG()
!write(mytid+600,'(a,i4)')"LICOM",__LINE__
#define LOGMSGCarbon()
!write(600+mytid,'(a,i4)')"Licom CARBON",__LINE__

!  CVS: $Id: licom.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
 program licom
!    *.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*
!    +                                                               +
!    + ============================================================= +
!    +         LASG/IAP  CLIMATE OCEAN MODEL Version2.0 (LICOM2.0)
!    + ============================================================= +
!
!             A primitive-equation ocean model upgraded by the
!
!    +    State Key Laboratory of Numerical Modelling for            +
!    +    Atmospheric Sciences and Geophysical Fluid Dynamics (LASG) +
!    +    INSTITUTE OF ATMOSPHERIC PHYSICS (IAP)                     +
!    +    CHINESE ACADEMY OF SCIENCE (CAS)                           +
!    +    P.O. Box 9804, Beijing 100029, P.R.China                   +
!    +                                                               +
!    +                   December, 2002                              +
!    +                                                               +
!    +   AUTHORS: Hailong LIU  (lhl@lasg.iap.ac.cn)                  +
!    +            Xuehong ZHANG (zxh@lasg.iap.ac.cn)                  +
!    +            Yongqiang YU (yyq@lasg.iap.ac.cn)                  +
!    +                                                               +
!    +                                                               +
!
!
!       This model is based upon, but differs substantially from the
!                                     ~~~~~~~~~~~~~~~~~~~~~
!
!                        LASG/IAP L30T63 OGCM
!                   ==================================
!
!                        which was implemented by
!                      Xiangze Jin and Xuehong ZHANG
!
!
!             +------------------------------------------+
!             | DISTRIBUTION TERMS AND CONDITIONS NOTICE |
!             +------------------------------------------+
!
! (c) Copyright 2002 State Key Laboratory of Numerical Modelling for
! Atmospheric Sciences and Geophysical Fluid Dynamics (LASG) /
! Institute of Atmospheric Physics  (IAP)
!
! This software, the LASG/IAP Climate Ocean Model (LICOM), version 2.0 , was
! upgraded by State Key Laboratory of Numerical Modelling for
! Atmospheric Sciences and Geophysical Fluid Dynamics (LASG), Institute
! of Atmospheric Physics  (IAP), Chinese Academy of Sciences (CAS)
!
! Access and use of this software shall impose the following obligations
! and understandings on the user.  The user is granted the right,
! without any fee or cost, to use, copy, modify, alter, enhance and
! distribute this software, and any derivative works thereof, and its
! supporting documentation for any purpose whatsoever, except commercial
! sales, provided that this entire notice appears in all copies of the
! software, derivative works and supporting documentation.  Further, the
! user agrees to credit LASG/IAP in any publications that result
! from the use of this software or in any software package that includes
! this software.  The names LASG/IAP, however, may not be used in
! any advertising or publicity to endorse or promote any products or
! commercial entity unless specific written permission is obtained from
! LASG/IAP.
!
! THE LICOM MATERIALS ARE MADE AVAILABLE WITH THE UNDERSTANDING THAT
! UCAR/NCAR/CGD IS NOT OBLIGATED TO PROVIDE (AND WILL NOT PROVIDE) THE
! USER WITH ANY SUPPORT, CONSULTING, TRAINING, OR ASSISTANCE OF ANY KIND
! WITH REGARD TO THE USE, OPERATION AND PERFORMANCE OF THIS SOFTWARE, NOR
! TO PROVIDE THE USER WITH ANY UPDATES, REVISIONS, NEW VERSIONS, OR "BUG
! FIXES."
!
! THIS SOFTWARE IS PROVIDED BY LASG/IAP "AS IS" AND ANY EXPRESS OR
! IMPLIED WARRANTIES, INCLUDING BUT NOT LIMITED TO, THE IMPLIED
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED.  IN NO EVENT SHALL UCAR/NCAR/CGD BE LIABLE FOR ANY
! SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER,
! INCLUDING BUT NOT LIMITED TO CLAIMS ASSOCIATED WITH THE LOSS OF DATA
! OR PROFITS, WHICH MAY RESULT FROM AN ACTION IN CONTRACT, NEGLIGENCE OR
! OTHER TORTIOUS CLAIM THAT ARISES OUT OF OR IN CONNECTION WITH THE
! ACCESS, USE OR PERFORMANCE OF THIS SOFTWARE.
!
!-----------------------------------------------------------------------
!
! Purpose: Entry point for LICOM
!
!-----------------------------NOTICE------------------------------------
!
!
! Method: Call appropriate initialization, time-stepping, and finalization routines.
!
! Author: Hailong Liu and Yongqiang Yu, Dec. 31, 2002
!
!-----------------------------------------------------------------------

#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use dyn_mod
#ifdef COUP
use shr_msg_mod
use shr_sys_mod
use control_mod
#endif
#if ( defined SPMD ) || ( defined COUP)
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
#endif
use tracer_mod
use pmix_mod
#ifdef USE_OCN_CARBON
use carbon_mod,only:ahvp
use cforce_mod
#endif
      IMPLICIT NONE
#include <netcdf.inc>
#ifdef TIDE_OUT
      character*13 fname
      integer irec
      real(r4) h0_diurnal(imt,jmt_global)
      real(r4) ub_diurnal(imt,jmt_global),vb_diurnal(imt,jmt_global)
      real(r4) us_diurnal(imt,jmt_global,km),vs_diurnal(imt,jmt_global,km)
#endif
!
!---------------------------------------------------------------------
!     Initilizing Message Passing
!---------------------------------------------------------------------
mpi_comm_ocn=0
#ifdef COUP
     write(6,*)"Begin shr_msg_stdio"
     call shr_msg_stdio('ocn')
     write(6,*)"End shr_msg_stdio"
     write(6,*)"Begin msg_pass ('connect')"
     call msg_pass('connect')
     write(6,*)"End msg_pass ('connect')"
#endif
!
      mytid=0
!
	LOGMSG()
#ifdef SPMD
#if ( ! defined COUP )
      call mpi_init(ierr)
      call mpi_comm_dup(mpi_comm_world, mpi_comm_ocn, ierr)
      if (mytid==0) write(*,*) "COMM",mpi_comm_world, mpi_comm_ocn, ierr
#endif
	LOGMSG()
      call mpi_comm_rank (mpi_comm_ocn, mytid, ierr)
      call mpi_comm_size (mpi_comm_ocn, nproc, ierr)
      if (mytid==0) write(6,*) "Number of Processors is",nproc
#endif
	LOGMSG()


!---------------------------------------------------------------------
!     SET THE CONSTANTS USED IN THE MODEL
!---------------------------------------------------------------------
#ifdef COUP
      call shr_sys_flush(6)
#endif
	LOGMSG()
      if (mytid==0) write(6,*) 'Begin CONST'
      
#ifdef COUP
      call shr_sys_flush(6)
#endif      
      CALL CONST
      if (mytid==0) write(6,*) 'END CONST'
#ifdef COUP
      call shr_sys_flush(6)
#endif
	LOGMSG()
	
#ifdef COUP
      call shr_sys_flush(6)
#endif
#ifdef SHOW_TIME
      call run_time('CONST')
#endif

!---------------------------------------------------------------------
!     SET MODEL'S RESOLUTION,TOPOGRAPHY AND THE CONSTANT
!     PARAMETERS RELATED TO LATITUDES (J)
!---------------------------------------------------------------------
	LOGMSG()
      if (mytid==0) write(6,*) 'Begin GRIDS'
      CALL GRIDS
      if (mytid==0) write(6,*) 'END GRIDS'
#ifdef COUP
      call shr_sys_flush(6)
#endif
#ifdef SHOW_TIME
      call run_time('GRIDS')
#endif

!***********************************************************************
!          SET SOME CONSTANTS FOR THE BOGCM
!***********************************************************************
#ifdef USE_OCN_CARBON
      LOGMSGCarbon()
      CALL CTRLC
      LOGMSGCarbon()
#ifdef SHOW_TIME
      call run_time('CTRLC')
#endif
#endif

!---------------------------------------------------------------------
!     SET SURFACE FORCING FIELDS (1: Annual mean; 0: Seasonal cycle)
!---------------------------------------------------------------------
	LOGMSG()
      if (mytid==0) write(6,*) 'Begin RDRIVER'
      CALL RDRIVER
      if (mytid==0) write(6,*) 'END RDRIVER'
#ifdef SHOW_TIME
      call run_time('RDRIVER')
#endif

!***********************************************************************
!      SET FORCING DATA USED IN CARBON CYCLYE
!***********************************************************************
#ifdef USE_OCN_CARBON
      LOGMSGCarbon()
       CALL CFORCE
      LOGMSGCarbon()
#ifdef SHOW_TIME
       call run_time('CFORCE')
#endif
#endif

!---------------------------------------------------------------------
!     INITIALIZATION
!---------------------------------------------------------------------
	LOGMSG()
      if (mytid==0) write(6,*) 'Begin INIRUN'
      CALL INIRUN
      if (mytid==0) write(6,*) 'END INIRUN'
#ifdef SHOW_TIME
      call run_time('INRUN')
#endif

!**********************************************************************
!      INITIALIZATION CARBON MODEL
!**********************************************************************
#ifdef USE_OCN_CARBON
      LOGMSGCarbon()
      CALL INIRUN_PT
      LOGMSGCarbon()
#ifdef SHOW_TIME
      call run_time('INIRUN_PT')
#endif
#endif

!---------------------------------------------------------------------
!     INITIALIZATION OF ISOPYCNAL MIXING
!---------------------------------------------------------------------
	LOGMSG()
#ifdef ISO
      CALL ISOPYI
#ifdef SHOW_TIME
      call run_time('ISOPYI')
#endif
#endif

      MEND = MONTH + NUMBER

      if (mytid==0) then
         WRITE (6,FMT='(A,I5,A,I5)') 'NUMBER=',NUMBER,',MONTH=',MONTH
!*************************************************************************
      IF(NSTART==1) THEN
        WRITE (6,FMT='(A)') 'The physical model is initial run!'
      ELSE
        WRITE (6,FMT='(A)') 'The physical model is continuous run!'
      ENDIF
#ifdef USE_OCN_CARBON
      WRITE (6,FMT='(A,I5,A,I5)') 'MONTHR=',monthR, &
                 ',yearR=', yearStart + (monthR-1)/12

      IF(NSTARTC==1) THEN
        WRITE (6,FMT='(A)') 'The passive-tracer model is initial run!'
      ELSE
        WRITE (6,FMT='(A)') 'The passive-tracer model is continuous run!'
      ENDIF

#endif
!*****************************************************************************

#ifdef COUP
         call shr_sys_flush(6)
#endif
      endif

#ifdef TIDE
!lhl20100801
         time_tidal=0.0
!lhl20100801
#endif

      loop1 : do

      IY0 = (MONTH -1)/12
      IYFM = IY0+1
!     IYFM is the number of the current year

      MON0 = MONTH - IY0*12
      IMD = NMONTH (MON0)
!     IMD is the total day number of the current month

      ISB = 0
      ISC = 0
      IST = 0
!     ISB/C/T : SWITCH ON EULER FORWARD OR LEAP-FROG SCHEME FOR
!     BAROTROPIC, BAROCLINIC, AND THERMOHALINE PROCESSES RESPECTIVELY

!**********************************************************************
#ifdef USE_OCN_CARBON
      yearR=(monthR-1)/12+1
      isp=0
#endif
!**********************************************************************
#ifdef TIDE_OUT
      if (mytid.eq.0) then
      fname(1:3)='z0_'
      write(fname(4:7),'(i4.4)') IYFM
      write(fname(8:9),'(i2.2)') MON0
      fname(10:13)='.dat'
      open (118,file=fname,access="direct",form="unformatted",recl=imt*jmt_global*4)
      fname(1:3)='ub_'
      open (119,file=fname,access="direct",form="unformatted",recl=imt*jmt_global*4)
      fname(1:3)='vb_'
      open (120,file=fname,access="direct",form="unformatted",recl=imt*jmt_global*4)
      fname(1:3)='us_'
      open (121,file=fname,access="direct",form="unformatted",recl=imt*jmt_global*4*30)
      fname(1:3)='vs_'
      open (122,file=fname,access="direct",form="unformatted",recl=imt*jmt_global*4*30)
      endif
#endif
!---------------------------------------------------------------------
!     THE CYCLE OF THE CURRENT MONTH
!---------------------------------------------------------------------
      loop2 : DO IDAY = 1,IMD
#ifdef TIDE
      if (mytid.eq.0) print*,"time_tidal= ",time_tidal*IDTB
#endif
!       
#ifdef COUP
      if (nstart==2.and.iday==1) then
         nstart = 0
         cycle loop2
      end if
      cdate=((month-mon0)/12+1)*10000+mon0*100+iday
      sec=0
!
	LOGMSG()
      if (mytid==0) write(*,*) 'Begin to Send to cpl'
      call msg_pass('send')
      if (mytid==0) write(*,*) 'Send to cpl OK'
      
	LOGMSG()
      call msg_pass('recv')
      if (mytid==0)write(*,*) 'Send to cpl OK'
	LOGMSG()
!
      call post_cpl
	LOGMSG()
#ifdef SPMD
      call mpi_bcast(stop_now,1,mpi_integer,0,mpi_comm_ocn,ierr)
#endif
      if (stop_now==1) exit loop1
#else
      if (month==mend) exit loop1
#endif
	LOGMSG()
      if (mytid==0) then
         OPEN (56,FILE ='modeltime',FORM ='FORMATTED',STATUS ='UNKNOWN')
         WRITE (56,FMT='(A,A,I4,I6)') ' It is running on ', ABMON (MON0),IDAY,1+ IY0
!**********************************************************************
#ifdef USE_OCN_CARBON
         WRITE (56,FMT='(A,I4,A,I6)') ' The ptracer model is running on ', &
                                      IDAY, ABMON (monthR-(yearR-1)*12),yearStart+yearR-1
#endif
!*********************************************************************
         CLOSE (56)
      endif
	LOGMSG()
!---------------------------------------------------------------------
!     INTERPOLATE THE OBSERVED MONTHLY MEAN DATA TO CERTAIN DAY
!---------------------------------------------------------------------
	LOGMSG()
#if (!defined COUP)
      CALL INTFOR
!         CALL INTFOR (IYFM,MON0)
#endif
!
!*********************************************************************
!   INTERPOLATE THE FORCING DATA TO CERTAIN DAY
!*********************************************************************
#ifdef USE_OCN_CARBON
#if (!defined COUP)
      LOGMSGCarbon()
      CALL INTFORW
      LOGMSGCarbon()
#ifdef SHOW_TIME
      call run_time('INTFORW')
#endif
#endif
      LOGMSGCarbon()
      CALL INTFOR_PT
      LOGMSGCarbon()
#ifdef SHOW_TIME
      call run_time('INTFOR_PT')
#endif
#endif

	LOGMSG()
#if (defined SOLARCHLORO)
       call sw_absor
! compute the short wave pentration dependent on chlorophy
#endif

#ifdef  SHOW_TIME
         call run_time('INTFOR')
#endif

!---------------------------------------------------------------------
!     THERMAL CYCLE
!---------------------------------------------------------------------

         DO II = 1,NSS

!lhl090724
! daily forcing for CORE data
!
        ! CALL CORE_DAILY(II)
!lhl090724

!     COMPUTE DENSITY, BAROCLINIC PRESSURE AND THE RELAVANT VARIABLES
	LOGMSG()
!	      if (mytid==0) write(6,*) 'Begin READYT'
            CALL READYT
!	      if (mytid==0) write(6,*) 'END READYT'
#ifdef  SHOW_TIME
            call run_time('READYT')
#endif

!---------------------------------------------------------------------
!     BAROCLINIC & BAROTROPIC CYCLE
!---------------------------------------------------------------------
            DO JJ = 1,NCC
!     COMPUTE MOMENTUM ADVECTION, DIFFUSION & THEIR VERTICAL INTEGRALS
	LOGMSG()
!	      if (mytid==0) write(6,*) 'Begin READYC'
               CALL READYC
!	      if (mytid==0) write(6,*) 'END READYC'
 	LOGMSG()
#ifdef  SHOW_TIME
               call run_time('READYC')
#endif


!     PREDICTION OF BAROTROPIC MODE
	LOGMSG()
               CALL BAROTR
#ifdef  SHOW_TIME
               call run_time('BAROTR')
#endif

!     PREDICTION OF BAROCLINIC MODE
	LOGMSG()
!	      if (mytid==0) write(6,*) 'Begin Bclinc'
               CALL BCLINC
!             if (mytid==0) write(6,*) 'End Bclinc'
#ifdef  SHOW_TIME
               call run_time('BCLINC')
#endif
!lhl
#ifdef  DEBUG
	LOGMSG()
         CALL ENERGY
#ifdef  SHOW_TIME
               call run_time('ENERGY')
#endif
#endif
!lhl

            END DO

	LOGMSG()
!*******************************************************************
!     PREDICTION OF PASSIVE TRACER
!*******************************************************************
#ifdef USE_OCN_CARBON
#ifdef carbonDebug
      print*, "II in licom======================================================================",II
#endif
#ifdef printcall
#ifdef SPMD
      print*,"call ptracer in licom, mytid=",mytid
#else
      print*,"call ptracer in licom"
#endif
#endif
      LOGMSGCarbon()
      CALL NEWPTRACER
      LOGMSGCarbon()

#ifdef SHOW_TIME
      call run_time('PTRACER')
#endif
#endif
	LOGMSG()
!	if (mytid==0) write(6,*) 'Begin Tracer'
            CALL TRACER
!        if (mytid==0) write(6,*) 'END Tracer'

	LOGMSG()
#ifdef  SHOW_TIME
            call run_time('TRACER')
#endif

	LOGMSG()
            CALL ICESNOW
#ifdef  SHOW_TIME
            call run_time('ICESNOW')
#endif

!***********************************************************************
!     PERFORM CONVECTIVE ADJUSTMENT IF UNSTABLE STRATIFICATION OCCURS
!************************************************************************
#ifdef USE_OCN_CARBON
      LOGMSGCarbon()
            CALL CONVADJ_PT
      LOGMSGCarbon()
#else
!lhl1204
!#if (!defined CANUTO)
            CALL CONVADJ
!#endif
#endif

#ifdef  SHOW_TIME
            call run_time('CONVADJ')
#endif
!*************************************************************************
!     ACCUMULATE SOME VARIABLES FOR MONTHLY OUTPUT
        LOGMSG()
         CALL ACCUMM
#ifdef  SHOW_TIME
         call run_time('ACCUMM')
#endif
!
#ifdef TIDE_OUT
         call local_to_global_4d_double(h0,h0_diurnal,1,1)
         call local_to_global_4d_double(ub,ub_diurnal,1,1)
         call local_to_global_4d_double(vb,vb_diurnal,1,1)
         call local_to_global_4d_double(u,us_diurnal,1,30)
         call local_to_global_4d_double(v,vs_diurnal,1,30)
         if (mytid.eq.0) then
         irec=II+(IDAY-1)*NSS
         write(118,rec=irec) ((h0_diurnal(i,j),i=1,imt),j=1,jmt_global)
         write(119,rec=irec) ((ub_diurnal(i,j),i=1,imt),j=1,jmt_global)
         write(120,rec=irec) ((vb_diurnal(i,j),i=1,imt),j=1,jmt_global)
         write(121,rec=irec) (((us_diurnal(i,j,k),i=1,imt),j=1,jmt_global),k=1,km)
         write(122,rec=irec) (((vs_diurnal(i,j,k),i=1,imt),j=1,jmt_global),k=1,km)
         endif
#endif
!
         END DO


!     COMPENSATE THE LOSS OF GROSS MASS
	LOGMSG()
         CALL ADDPS
#ifdef  SHOW_TIME
         call run_time('ADDPS')
#endif
!
	LOGMSG()
         CALL ENERGY
#ifdef  SHOW_TIME
         call run_time('ENERGY')
#endif
!
!     MONITOR THE MODEL INTEGRATION

!     ACCUMULATE SOME VARIABLES FOR MONTHLY OUTPUT
!	LOGMSG()
!         CALL ACCUMM
!#ifdef  SHOW_TIME
!         call run_time('ACCUMM')
!#endif
!*********************************************************************
!      ACCUMULATE SOME VARIABLES IN CARBON MODEL
!*********************************************************************
#ifdef USE_OCN_CARBON
      LOGMSGCarbon()
      CALL ACCUMM_PT
      LOGMSGCarbon()
#ifdef SHOW_TIME
      call run_time('ACCUMM_PT')
#endif
#endif

#ifdef COUP
         call shr_sys_flush(6)
#endif
!
	LOGMSG()
#ifdef COUP
         call flux_cpl
#endif

      END DO loop2

#ifdef TIDE_OUT
         if (mytid.eq.0) then
         close(118)
         close(119)
         close(120)
         close(121)
         close(122)
         endif
#endif

	LOGMSG()
      MONTH = MONTH +1
#ifdef USE_OCN_CARBON
      MONTHR=MONTHR+1
#endif

!---------------------------------------------------------------------
!     SAVE, SMOOTH & RESET SOME ARRAYS
!---------------------------------------------------------------------
	LOGMSG()
      CALL SSAVECDF
#ifdef  SHOW_TIME
      call run_time('SSAVECDF')
#endif
!*********************************************************************
!     SAVE, SMOOTH & RESET SOME ARRAYS IN CARBON MODEL
!*********************************************************************
#ifdef USE_OCN_CARBON
#ifdef printcall
#ifdef SPMD
      print*,"call ssave_pt in licom, mytid=",mytid
#else
      print*,"call ssave_pt in licom"
#endif
#endif
      CALL SSAVE_PT
#ifdef SHOW_TIME
      call run_time('SSAVE_PT')
#endif
#endif
      end do loop1

#ifdef COUP
      call msg_pass('disconnect')
#else
#ifdef SPMD
      call mpi_finalize(ierr)
#endif
#endif
      close(6)
      STOP
      END  program LICOM
