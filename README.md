Introduction
============
LICOM is LASG/IAP Climate system Ocean model
http://project.lasg.ac.cn/LFS/

Deployment steps
============
1.Granting executable permission to file case.sh 
```
chmod 777 ./bld/case.sh
```
2.MODIFY the parameters in case.sh
```
vim ./bld/case.sh
```
3.CREATE AN EXPERIENCE named $CASEAME at the location ./$CASENAME
```
./bld/case.sh
```
4.run the experimental
```
./$CASEAME/exe/run
```
or by the way of $ to submit to the backend (Screen output is redirected to ./$CASEAME/exe/screen)
```
cd ./$CASEAME/src/
make run 
```

Dependency library
============
1.cpp (PRECOMPILE .F90 TO .f90 FILES)
.e.g. GNU Compiler Collection

2.mpi_fortran (COMPILE .f90 FILES)
.e.g. Intel® oneAPI Base Toolkit + Intel® HPC Toolkit(High-Performance Computing)
.e.g. GNU Compiler Collection + Openmpi

3.netcdf for C and FORTRAN
.e.g. https://www.unidata.ucar.edu/software/netcdf/
.e.g. git clone https://github.com/dongli/starman
