#!/bin/csh -f
echo '-----------------------------------------------------------------'
echo '     Set env variables available to model setup scripts (below): '
echo '-----------------------------------------------------------------'

setenv CASENAME wind0
setenv LICOMROOT /data06/yyq/mls/licom2
setenv SRCPATH  $LICOMROOT/src
setenv BLDPATH  $LICOMROOT/bld
setenv DATAPATH $LICOMROOT/data
setenv EXEROOT  $LICOMROOT/$CASENAME
setenv EXESRC   $EXEROOT/src
setenv EXEDIR   $EXEROOT/exe
setenv RUNTYPE  initial                   # run type:continue or initial
set HISTOUT = 1                           # model historic output, no use for this version !!!!!
set RESTOUT = 1                           # model restar file outpur every day
set NTASKS = 8
set NTHRDS  = 1
set LID = "`date +%y%m%d-%H%M%S`"    


echo '-----------------------------------------------------------------'
echo '     Copy the source codes                                       ' 
echo '-----------------------------------------------------------------'

mkdir -p  $EXEROOT
cd        $EXEROOT
mkdir -p  $EXESRC
mkdir -p  $EXEDIR
cp    -pf $SRCPATH/* $EXESRC/.
cp    -pf $BLDPATH/Makefile $EXESRC/Makefile
#
cd $EXESRC
#
echo '-----------------------------------------------------------------'
echo '     Produce Makefile                                            '
echo '-----------------------------------------------------------------'


echo '-----------------------------------------------------------------'
echo '     Produce the pre-compile file def-undef h                    '
echo '-----------------------------------------------------------------'

\cat >! def-undef.h << EOF
#define N_PROC  $NTASKS 
#define SPMD
#define  SYNCH
#undef  FRC_ANN
#define CDFIN
#undef  FRC_DAILY
#define SOLAR
#define  ACOS
#undef  BIHAR
#undef  SMAG_FZ
#undef  SMAG_OUT
#define NETCDF
#undef BOUNDARY
#define NODIAG
#undef  ICE
#undef SHOW_TIME
#undef DEBUG
#undef COUP
#define  ISO
#define D_PRECISION
#undef CANUTO
#undef SOLARCHLORO
#define LDD97
#undef TSPAS
#undef  SMAG
#define JMT_GLOBAL 115
#undef TIDE
#undef TIDE_OUT
EOF

echo '-----------------------------------------------------------------'
echo '     Compile and Link                                            '
echo '-----------------------------------------------------------------'

make clean 
make > makelog.$LID

if ( ! -e ../exe/licom2 ) then
echo "compile failure!"
exit 1
endif

echo '-----------------------------------------------------------------'
echo '     Produce the namelist file                                   '
echo '-----------------------------------------------------------------'

if      ($RUNTYPE == 'initial' ) then
  set NSTART = 1
else if ($RUNTYPE == 'continue') then
  set NSTART = 0
endif
#
cd $EXEDIR
#
\cat >! ocn.parm << EOF
 &namctl
  DLAM       = 2.0            !grid distance
  AM_TRO     = 15000
  AM_EXT     = 15000
  IDTB       = 120
  IDTC       = 2880
  IDTS       = 5760
  AFB1       = 0.20
  AFC1       = 0.43
  AFT1       = 0.43
  AMV        = 1.0E-3
  AHV        = 0.3E-4
  NUMBER     = 24
  NSTART     = $NSTART
  klv        = 30
  IO_HIST = $HISTOUT
  IO_REST = $RESTOUT
  diag_bsf = .true.
  diag_msf = .true.


 &end
EOF

echo '-----------------------------------------------------------------'
echo '     Link the data files to excutive directory                   '
echo '-----------------------------------------------------------------'

ln -s $DATAPATH/INDEX.DATA INDEX.DATA
ln -s $DATAPATH/TSinitial  TSinitial
ln -s $DATAPATH/MODEL.FRC  MODEL.FRC
ln -s $DATAPATH/dncoef.h1 dncoef.h1
 
echo '-----------------------------------------------------------------'
echo '     Run the model using run script                              '
echo '-----------------------------------------------------------------'
\cat >! run << EOF
mpirun -n $NTASKS ./licom2
#bsub -J control1 -a intelmpi -n $NTASKS  -q hpc_linux -o out1_1.lsf mpirun.lsf ./licom2
#bsub -R "span[ptile=4]" -J control1 -a intelmpi -n $NTASKS  -q hpc_linux -o out1_1.lsf mpirun.lsf ./licom2
EOF

chmod u+x run
