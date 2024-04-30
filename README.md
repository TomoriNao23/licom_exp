#THIS IS LICOM MODEL

#YOU CAN RUN THE MODEL BY 
chmod 777 ./bld/case.sh

#MODIFY the parameters in ./bld/case.sh, and then CREATE AN EXPERIENCE named $CASEAME at the location ./$CASENAME
vim ./bld/case.sh
./bld/case.sh

#run the experimental
./$CASEAME/exe/run 
#or by the way of $ to submit to the backend
cd ./$CASEAME/src/
make run (#Screen output is redirected to ./$CASEAME/exe/screen)


