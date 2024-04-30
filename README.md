Introduction
============
LICOM is LASG/IAP Climate system Ocean model
http://project.lasg.ac.cn/LFS/

Deployment steps
============
## 1.Grant Executable Permission to File case.sh 
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
- run script
```
./$CASEAME/exe/run
```
- by the way of $ to submit to the backend (Screen output is redirected to ./$CASEAME/exe/screen)
```
cd ./$CASEAME/src/
make run 
```

Dependency library
============
## cpp (PRECOMPILE .F90 TO .f90 FILES)   
- GNU Compiler Collection
## mpi_fortran (COMPILE .f90 FILES)  
- Intel® oneAPI Base Toolkit + Intel® HPC Toolkit(High-Performance Computing)  
- GNU Compiler Collection + Openmpi
## netcdf for C and FORTRAN  
- https://www.unidata.ucar.edu/software/netcdf/    
- starman (Another package manager for HPC warriors)
```
git clone https://github.com/dongli/starman
```
## Makefile

Future Work
============
- use CMAKE 
- use CUDAC to Implemente GPU parallelism with CPU parallelism
