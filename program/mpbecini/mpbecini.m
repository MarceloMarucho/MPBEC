cual=computer;
ese=pwd
addpath(genpath(ese));
if ((strcmp(cual,'PCWIN')==1) ||(strcmp(cual,'PCWIN64')==1))
cd ..
cd src
cd Main
else
cd ('/home')
mpbecdir=what('mpbecini');
cd (mpbecdir.path);
cd ..
cd src
cd Main
end
clc
run Welcome
movegui(Welcome,'center')
