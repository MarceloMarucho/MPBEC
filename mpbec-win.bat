dir mpbecini /s /b > tmpFile 
set /p mpbecdir= < tmpFile 
del tmpFile
chdir /d %mpbecdir%
matlab -nosplash -nodesktop -r "run mpbecini"
