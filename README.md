Introduction
============
LICOM is LASG/IAP Climate system Ocean model
http://project.lasg.ac.cn/LFS/

Deployment steps
============
## 1.Grant Executable Eermission to File case.sh 
```
chmod 777 ./bld/case.sh
```
## 2.Modify the Parameters in case.sh
```
vim ./bld/case.sh
```
## 3.CREATE an Experiment Named $CASEAME at the Location ./$CASENAME
```
./bld/case.sh
```
## 4.Run the Experiment
```
- ./$CASEAME/exe/run
```
or by the way of $ to submit to the backend (Screen output is redirected to ./$CASEAME/exe/screen)
```
- cd ./$CASEAME/src/
- make run 
```

Dependency library
============
## 1.cpp (PRECOMPILE .F90 TO .f90 FILES)   
- GNU Compiler Collection

## 2.mpi_fortran (COMPILE .f90 FILES)  
- Intel® oneAPI Base Toolkit + Intel® HPC Toolkit(High-Performance Computing)  
- GNU Compiler Collection + Openmpi

## 3.netcdf for C and FORTRAN  
- https://www.unidata.ucar.edu/software/netcdf/    
- git clone https://github.com/dongli/starman

## 4.Makefile

Future Work
============
- use CMAKE 
- use CUDAC to Implemente GPU parallelism with CPU parallelism
