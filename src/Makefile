#INCLDIR  = -I. -I/usr/include -I/usr/local/include -I /opt/intel/impi/4.1.1.036/intel64/include -I/public/yyqtest/tools/netcdf-4.3/include
#CPPFLAGS = -P -traditional
#CPP      = /usr/bin/cpp $(CPPFLAGS) $(INCLDIR)
#NLIB =   -L/opt/intel/impi/4.1.1.036/intel64/lib -L/public/yyqtest/tools/netcdf-4.3/lib -lnetcdff
#FFLAGS = -O3 -mcmodel=medium  -shared-intel -r8  -i4 -convert big_endian -assume byterecl -save -no-vec  -traceback -heap-arrays  64
#FC =  /opt/intel/impi/4.1.1.036/intel64/bin/mpiifort \$(FFLAGS)
INCLDIR  = -I. -I/usr/include -I/usr/local/include -I /opt/intel/impi/5.0.2.044/intel64/include -I/public/yyq/tools/netcdf3/include
#INCLDIR  = -I. -I/usr/include -I/usr/local/include -I/opt/intel/impi/5.0.2.044/intel64/include -I/usr/local/netcdf3/include
CPPFLAGS = -P -traditional
CPP      = /usr/bin/cpp $(CPPFLAGS) $(INCLDIR)
NLIB =   -L/opt/intel/impi/5.0.2.044/intel64/lib -L/public/yyq/tools/netcdf3/lib -lnetcdf
#NLIB =   -L/opt/intel/impi/5.0.2.044/intel64/lib -L/usr/local/netcdf3/lib -lnetcdf
#FFLAGS = -O3 -mcmodel=medium  -shared-intel -r8  -i4 -convert big_endian -assume byterecl -save -no-vec  -traceback -heap-arrays  64
FFLAGS =  -mcmodel=large -shared-intel -O0 -r8 -i4  -convert big_endian -assume byterecl -no-vec
#FC =  /opt/intel/impi/4.1.1.036/intel64/bin/mpiifort \$(FFLAGS)
FC = /opt/intel/impi/5.0.2.044/intel64/bin/mpiifort \$(FFLAGS)
OBJS   = precision_mod.o param_mod.o pconst_mod.o  pmix_mod.o isopyc_mod.o cdf_mod.o\
tracer_mod.o dyn_mod.o work_mod.o output_mod.o forc_mod.o msg_mod.o sw_mod.o\
accumm.o diag_mod.o  const.o    filter.o   invtri.o   k1_3.o     rdriver.o   vinteg.o\
addps.o    convadj.o  grids.o    isoadv.o   k2_3.o     readyc.o   smuvh.o    yy00.o\
barotr.o   dens.o     icesnow.o  isoflux.o  k3_123.o   readyt.o   setbcx.o\
bclinc.o   density.o  inirun.o   isopyc.o   licom.o     tracer.o\
energy.o   intfor.o   isopyi.o   mm00.o     upwell.o  \
ssave-cdf.o smag.o  boundary.o global_to_local.o \
exchange.o local_to_global.o clock0f.o run_time.o sw_absor.o thermalexpansion.o canuto_new.o core_daily.o
#icesnow1.o comqd.o diag_mod.o
###################################################################################
licom2: $(OBJS)
	$(FC)  -o ../exe/licom2  $(OBJS) $(NLIB)

.SUFFIXES : .F90
.F90.o:
	$(CPP) $< > $*.f90  ;  $(FC) -c  $*.f90  ;  mv  $*.f90 $*.i
.c.o:
	cc  -c $<
clean:
	rm -f *.o *.i *.il *.f90 *.f *.mod *.d licom2
