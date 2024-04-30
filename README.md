Introduction
============
LICOM is LASG/IAP Climate system Ocean model
http://project.lasg.ac.cn/LFS/

Deployment steps
============
1.Add file execution permissions 
```
chmod 777 ./bld/case.sh
```

2.MODIFY the parameters in ./bld/case.sh, and then CREATE AN EXPERIENCE named $CASEAME at the location ./$CASENAME
```
vim ./bld/case.sh
./bld/case.sh
```

3.run the experimental
```
./$CASEAME/exe/run
```
or by the way of $ to submit to the backend
```
cd ./$CASEAME/src/
make run (#Screen output is redirected to ./$CASEAME/exe/screen)
```

Dependency library
============
1.cpp (PRECOMPILE.F90 TO .f90)
.e.g. GNU Compiler Collection

2.mpi_fortran (Compile .f90 files)
.e.g. Intel® oneAPI Base Toolkit + Intel® HPC Toolkit(High-Performance Computing)
.e.g. GNU Compiler Collection + Openmpi

3.netcdf for C and FORTRAN
.e.g. https://www.unidata.ucar.edu/software/netcdf/
.e.g. git clone https://github.com/dongli/starman
