% defining the name of the output files depending on the working grid.
if (dielw~=dielp)
    
if strcmp(bc, 'coarse')==1
    plotname='coarse_grid_MATLAB_pot';
    outputfile='coarse_grid_MATLAB_pot.dx';
    outputfile3='coarse_grid_dielx.dx';
    outputfile4='coarse_grid_diely.dx';
    outputfile5='coarse_grid_dielz.dx';
    outputfile6='coarse_grid_kappa.dx';
    outputfile2='coarse_grid_MATLAB_rho.dx';
% let's make a copy of the input file in the current directory. It is needed
% if the focus boundary condition is required.
suu = char( which(inputfile,'-all'));
 copyfile(suu,'copyinputfile.inm');
else
    outputfile='target_grid__MATLAB_pot.dx';
    outputfile2='target_grid_MATLAB_rho.dx';
    outputfile3='target_grid_dielx.dx';
    outputfile4='target_grid_diely.dx';
    outputfile5='target_grid_dielz.dx';
    outputfile6='target_grid_kappa.dx';
    
    plotname='target_grid_MATLAB_pot';
end

else
    
if strcmp(bc, 'coarse')==1
    plotname='ref_coarse_grid_MATLAB_pot';
    outputfile='ref_coarse_grid_MATLAB_pot.dx';
    outputfile3='ref_coarse_grid_dielx.dx';
    outputfile4='ref_coarse_grid_diely.dx';
    outputfile5='ref_coarse_grid_dielz.dx';
    outputfile6='ref_coarse_grid_kappa.dx';
    outputfile2='ref_coarse_grid_MATLAB_rho.dx';
% let's make a copy of the input file in the current directory. It is needed
% if the focus boundary condition is required.
suu = char( which(inputfile,'-all'));
 copyfile(suu,'copyinputfile.inm');
else
    outputfile='ref_target_grid__MATLAB_pot.dx';
    outputfile2='ref_target_grid_MATLAB_rho.dx';
    outputfile3='ref_target_grid_dielx.dx';
    outputfile4='ref_target_grid_diely.dx';
    outputfile5='ref_target_grid_dielz.dx';
    outputfile6='ref_target_grid_kappa.dx';
    
    plotname='ref_target_grid_MATLAB_pot';
end

end

% checking if the user is working on a pc o unix platarform
kt = strfind(nam_str, '/');

if numel(kt)>0 
    outpath=strcat(nam_str,'/',dirname);
else
    outpath=strcat(nam_str,'\',dirname);
end

% saving a copy of the input files in the output folder

if numgrid==1
 sus = char( which(inputfile));
 copyfile (sus, outpath);
 axx=strcat('program/src/Main/',inputfile);
 if (exist(axx, 'file')==2)
     delete(inputfile)
 end
 else
    sus = char( which('focusname.inm'));
    copyfile (sus, outpath);
     if (exist('program/src/Main/focusname.inm', 'file')==2)
     delete('focusname.inm')
     end
end
 
if (wout(4)==1.0)
    
disp('writing out the electrostatic potential solution')
% electrostatic potential solution

dxformat=MATLAB_pot;
namefile='POTENTIAL (KT/e)';
run dx_export
disp(['the file ' outputfile ' was generated'])
%disp('Done!....')
disp(' ')

% let's create a copy of the matlab solution which is required for the
% focus boundaary condition as the input file below
if numgrid ==1
copyfile(outputfile,'MATLAB_Solution.dx');
end
%lets make a copy if we have to generate plots using Jmol
if (wout(5)==1.0)
    copyfile(outputfile,'MEP.dx');
end

end % wout(4)

movefile (outputfile, outpath)

if (wout(1)==1.0)
    
% charge map
disp('writing out the charge map')
dxformat=charge;
namefile='CHARGE DENSITY (e/A^3)';
outputfile=outputfile2;
run dx_export
disp(['the file ' outputfile ' was generated'])
movefile (outputfile, outpath)
%disp('Done!')
disp(' ')

end

if (wout(3)==1.0)
    
rmin=[(xmin+0.5*h(1)) ymin zmin];

disp('writing out the dielx map')
dxformat=dielx;
namefile='dielx';
outputfile=outputfile3;
run dx_export
disp(['the file ' outputfile ' was generated'])
movefile (outputfile, outpath)
%disp('Done!')
disp(' ')

rmin=[xmin (ymin+0.5*h(2)) zmin];
disp('writing out the diely map')
dxformat=diely;
namefile='diely';
outputfile=outputfile4;
run dx_export
disp(['the file ' outputfile ' was generated'])
movefile (outputfile, outpath)
%disp('Done!')
disp(' ')

rmin=[xmin ymin (zmin+0.5*h(3))];
disp('writing out the dielz map')
dxformat=dielz;
namefile='dielz';
outputfile=outputfile5;
run dx_export
disp(['the file ' outputfile ' was generated'])
movefile (outputfile, outpath)
%disp('Done!')
disp(' ')

end 

if (wout(2)==1.0)
    
rmin=[xmin ymin zmin];
disp('writing out the ionic accessible map')
dxformat=kappa;
namefile='kappa';
outputfile=outputfile6;
run dx_export
disp(['the file ' outputfile ' was generated'])
movefile (outputfile, outpath)
%disp('Done!')
disp(' ')

end

if (wout(5)==1.0)
% generating plots
disp('Generating surface plots!....Please, close Jmol application to continue with the calculations')
disp(' ')
folderdx=pwd;
foldermain=folderdx;
folderpqr=in_nam_str;
path1=strcat(folderpqr,'/',pqr_str);
path2=strcat(foldermain,'/','mpbecinput.pqr');
copyfile(path1,path2);
status = system('java -jar Jmol.jar --silent -s jmolscript.st');
delete('mpbecinput.pqr');
delete('MEP.dx');
delete('jmolsurface.jvxl');
movefile('jmol_visualization_files.zip',outpath);


%disp(['the file ' outputfile ' was generated'])
%disp('Done!....')
%disp(' ')
 n=int16((dime(3)+1)/2);
 
 figure
 
 %plotting the electrostatic potential solutions
 name2= plotname;
 plot2=surf(MATLAB_pot(:,:,n),'FaceColor','interp','EdgeColor','none','FaceLighting','phong');
 % I do not want MATLAB to create a figure
 % get(0,'CurrentFigure')
 saveas(plot2,name2,'fig');
 disp(['the Matlab figure file ' name2 '.fig was generated'])
 movefile (strcat(name2,'.fig'), outpath)
 saveas(plot2,name2,'jpg');
 disp(['the picture file ' name2 '.jpg was generated'])
 movefile (strcat(name2,'.jpg'), outpath)
 disp(' ')
 disp('Generating surface plots!....Please, close Matlab Figure application to finish MAPBS')
disp(' ')
% 
%pause (3)
%close
disp('Done!....')
disp(' ')
end
