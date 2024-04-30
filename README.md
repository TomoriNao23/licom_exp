Introduction
============
LICOM is LASG/IAP Climate system Ocean model http://project.lasg.ac.cn/LFS/  
This MODEL licom_exp is simplified by LICOM2.0  
## THIS MODEL IS ONLY FOR THE COURSE STUDYING OF <<NUMERICAL EXPERIMENTS FOR OCEANIC CIRCULATION AND AIR-SEA INTERACTION>>

Deployment steps
============
## 1.Grant Executable Permission 
to CSHELL SCRIPT case.sh 
```
chmod 777 ./bld/case.sh
```
## 2.Modify the Parameters 
in case.sh
```
vim ./bld/case.sh
```
## 3.Create an Experiment 
Named $CASEAME at the Location ./$CASENAME
```
./bld/case.sh
```
## 4.Run the Experiment
- SCEIPT run
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
