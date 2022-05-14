function [] = MPBEC(inputfile)
clc
global enerT outpath exitmpbec

%% MPBEC script to solve the linearized PB equation

% READ THE GUIDELINE DOCUMENT (PDF FILE) FOR MORE DETAILS.

% cd ..
% ese=pwd;
% addpath(genpath(ese));
% cd Main

%% Part 1.  Read the data
 
% read the inputfile.inm file and display it in the matlab cmd

 [dime, glen, wout, bulk, bc, digpres, selectsolver, dielx_str, diely_str, dielz_str, kappa_str, charge_str, pqr_str,pqr_cent_str,srfm, chgm, ener_str,in_nam_str, nam_str] = read_inm(inputfile);
 addpath(in_nam_str)

% creating the output folder

dirnametemp=strcat('outputfiles_',inputfile);
dirname = strtok(dirnametemp,'.');
mkdir(nam_str,dirname)

% creating the file which will contain all messages printed on the screen

logfile=strcat(inputfile,'.log');
diary (logfile)

disp('WELCOME!!!!!!....')
disp(' ')
disp('This code approximately solves the linear PB equation')
disp('for the electrostatic potential on a 3D grid.')
disp(' ')

% test on required memory and available RAM memory
disp('Checking available RAM memory....');
disp(' ');

run memory_check

% quit MPBEC ?
 
if exitmpbec==1 
    exitmpbec=0;
    close all force
    clear
    clc
    return
end

% if we use focus boundary condition, numgrid=1 evaluates the elect pot in
% the coarse grain and then numgrid=2 evaluates the boundary
% condition on the fine grid to obtain the solution in a subdomain. 
% If we use just one grid, then numgrid=1 refers to the fine
% grid and the numgrid=2 is not used.

oldbc=bc;
  
  
for numgrid=1:2

% if the bc=focusname.inm in the input file, I have to perfrom the calculation on
% the coarse grid and then the solution will be used to evalaute the boundary condition
% on the finest grid. Hence, In such case, I have to read the focusname.inm
% file.


if strcmp(bc, 'focusname.inm')==1
 %   read the inputfile
 
    [dime, glen, wout, bulk, bc, digpres, selectsolver, dielx_str, diely_str, dielz_str, kappa_str, charge_str, pqr_str,pqr_cent_str,srfm, chgm, ener_str,in_nam_str, nam_str] = read_inm(bc);

disp('Coarse grid calculation')
disp(' ')
addpath(in_nam_str)
end

% evaluating parameters such as zmagic, debye constant, kappa, thermal unit, etc.
run parameters

if (num2str(dielw)==num2str(dielp))
    disp('Calculations in the reference state !!!')
    disp(' ')
else
    disp('Calculations in the solvated state !!!')
    disp(' ')
end

disp(['Reading the input file ', inputfile])
disp(' ')
disp('GRID')
disp(' ')
disp(['number of points=   ',num2str(dime),',    box size (in A)=  ',num2str(glen), ',    Number of atoms=  ',num2str(atomN)])
disp(' ')
disp('PARAMETERS')
disp(' ')
disp(['Solvent dielectric coefficient = ',num2str(dielw),',   Solute dielectric coefficient = ',num2str(dielp),',  Temperature (K) = ',num2str(T)])
disp(' ')


%% Part 2.  solve A*pot=b

% initialize the timer
tic

disp('PRELIMINARLY CALCULATIONS')
disp(' ')
disp(['Calculating the ',bc,' boundary condition ....'])
disp(' ')

% quit MPBEC ?
 
if exitmpbec==1 
    exitmpbec=0;
    close all force
    clear
    clc
    return
end


if numgrid==1
run BoundaryCondition
%change the assigned value for bc in order to go to the fine grid
%calculation
if strcmp(oldbc, 'focusname.inm')==1
bc='coarse';
end

else 
% if numgrid=2 I dont have to evalaute bc using Dirichelet expressions.
% Instead, I have to use the solution on the coarse grid to evalaute the boundary condition on the finest target
% grid using linear interpolation
run fbc
end

% disp('Done!....')
disp(' ')

% quit MABPS ?
pause(1);
% run usedmemory
if exitmpbec==1 
    exitmpbec=0;
    close all force
    clear
    clc
    return
end

% evaluation / reading of the coeficient maps

if (strcmp(dielx_str,'nouploadedfiles')==1)

disp(['Generating the charge map using ',chgm, ' ....'])
disp(' ')
run chargemap

% disp('Done!....')
disp(' ')

disp(['Generating the dielectric and ionic maps using ',srfm, ' ....'])
disp(' ')
disp(['molecular solvent radius (A) = ',num2str(srad),',    Surface density (1/A^2) = ',num2str(surf_density),'    and Swin (A) = ',num2str(splineWin)])
disp(['Ionic Strength = ',num2str(bulkIonicStrength),',        Maximum Ionic Radius (A) = ',num2str(maxradprobe)])

% kappa and diel map calculations depend on the surface
% method (MOL, SMOL, SPL2)

% quit MABPS ?
pause(1);
% run usedmemory
if exitmpbec==1 
    exitmpbec=0;
    close all force
    clear
    clc
    return
end

if (strcmp(srfm,'MOL') ==1)
    
% calculate dielectric maps
    
run dielmap % fillcoCoefMOlDielNoSmooth

% quit MPBEC ?
pause(1);
% run usedmemory
if exitmpbec==1 
    exitmpbec=0;
    close all force
    clear
    clc
    return
end

% calculate kappa maps

run kappamap % fillcoCoefIon

% quit MPBEC ?
pause(1);
% run usedmemory
if exitmpbec==1 
    exitmpbec=0;
    close all force
    clear
    clc
    return
end

end

if (strcmp(srfm,'SMOL') ==1)
    
% calculate dielectric maps
    
run dielmap

% quit MPBEC ?
pause(1);
% run usedmemory
if exitmpbec==1 
    exitmpbec=0;
    close all force
    clear
    clc
    return
end

% calculate kappa maps

run kappamap % fillcoCoefIon

% quit MPBEC ?
pause(1);
% run usedmemory
if exitmpbec==1 
    exitmpbec=0;
    close all force
    clear
    clc
    return
end

if (dielw~=dielp) %(solvated state)
run harmonicsmoothing % fillcoCoefMolDielSmooth
end

end

% quit MPBEC ?
pause(1);
% run usedmemory
if exitmpbec==1 
    exitmpbec=0;
    close all force
    clear
    clc
    return
end

if (strcmp(srfm,'SPL2') ==1)
%calculate dielectric and kappa using swin window (Roux's approach)

run cubicspline % fillcoCoefSpline
end

% disp('Done!....')
disp(' ')

else
    
disp('reading the maps...')

% read the .dx files
dielx=data_parse(dielx_str, dime);
diely=data_parse(diely_str, dime);
dielz=data_parse(dielz_str, dime);
kappa=data_parse(kappa_str, dime);
charge=data_parse(charge_str, dime);
disp('Done!....')
disp(' ')
end

%temp=toc
%return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%temporary lines to test maps
% run outputtestmaps
% return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('APPROXIMATE SOLUTION FOR THE ELECTROSTATIC POTENTIAL')
disp(' ')

% Prepare the Laplacian operator A and b

disp('Constructing the sparse matrix A....')

% quit MPBEC ?
pause(1);
% run usedmemory
if exitmpbec==1 
    exitmpbec=0;
    close all force
    clear
    clc
    return
end

run BuildA

disp('Done!....')
disp(' ')

% quit MPBEC ?
pause(1);
% run usedmemory
if exitmpbec==1 
    exitmpbec=0;
    close all force
    clear
    clc
    return
end

% disp('Performing the LU decomposition....')
%  setup.type='nofill';
%  setup.milu='row';
% % 
%  [L U]=ilu(A,setup);
% 
% disp('Done!....')
% disp(' ')
% 
% disp('Solving the linearized PB equation')
%  setup.milu='off';
% [LL UU]=ilu(A,setup);
%  potein = UU\(LL\bb);
% disp('Done!....')
% disp(' ')

accuracy=10^-(digpres);
max_iteration=5800;

% Solve A*pot=b to obtain an approximate solution for "pot"
 

 run linearsolver

 % quit MPBEC 
 pause(1);
 % run usedmemory
if exitmpbec==1 
    exitmpbec=0;
    close all force
    clear
    clc
    return
end

% Stop the calculations if the linear solver coud not obtain the solution
% for "pot"
% 
if fail==1 
return
end

%%

% Add Boundary to Solution

run addboundcond


% 
disp(['Residual error=  ',num2str(errorh),',     Iteration number=  ',num2str(iteration_number)])
disp(' ')

%
if strcmp(ener_str, 'calceneryes')==1
    disp('Calculating the energy...')
run energy 

end
disp('END OF THE CALCULATIONS')
disp(' ')

toc
disp(' ')
%% Part 3: Write out the electrostatic potential and maps in dx format and
%% generate plots (also calculate energy if needed)

%make sure the user required output files

if (sum(wout) > 0.0)

run writeout

disp(['These files were saved in the directory ', outpath,' as required'])
disp(' ')

else
    disp('The user did not require output files')
end


%% Part 4: checking for the second grid

% If bc=focusname.inm I already evalauted the elect pot
% in the coarse grid and now I have to read the input file inputfile.inm to evalaute the
% elec pot in the fine grid. Otherwise I am done and I must
% exit the numgrid loop.

if (strcmp(bc, 'sdh')==0 && strcmp(bc, 'periodic')==0 && strcmp(bc, 'mixed')==0 && strcmp(bc, 'mdh')==0&& strcmp(bc, 'from the coarse grid solution')==0)

pause (3)

disp('Calculating')
disp('the electrostatic potential in the target grid....')
disp(' ')
inputfile='copyinputfile.inm';

disp('reading the solution obtained previously in the coarse grid....')
disp(' ')
[rminn,cgdime,cgh]=gridinf('MATLAB_Solution.dx');
MATLAB_pot_coarse=data_parse('MATLAB_Solution.dx', cgdime);
disp(' ')
  [dime, glen, wout, bulk, bc, digpres, selectsolver, dielx_str, diely_str, dielz_str, kappa_str, charge_str, pqr_str,pqr_cent_str,srfm, chgm, ener_str,in_nam_str, nam_str] = read_inm(inputfile);
     bc='from the coarse grid solution'; 
 %    THE SAME Boundary Condition 'bc' I USED TO EVALUATE THE POT IN THE COARSE GRID is used to
 %    solve the LPBE in the target grid!!
 
delete('copyinputfile.inm');
%  sus = char( which('focusname.inm','-all'));
%  movefile (sus, outpath)
addpath(in_nam_str)
else
% exit the numgrid loop
    break
end
end % end numgrid loop

disp('The information displayed on this screen is saved in the file MATLAB_screen.io')
disp(' ')
disp('Thanks for using MPBEC!!!!....')
diary off
movefile (logfile, outpath)
propkaname=regexprep(pqr_str,'pqr','propka');

aaa=dir; 
bbb=struct2cell(aaa);
if any(ismember(bbb(1,:),propkaname))
movefile (propkaname, outpath)
movefile ('pKa.pqr', outpath)
movefile ('input.pqr', outpath)
end

% if exist(propkaname,'file')
% movefile (propkaname, outpath)
% movefile ('pKa.pqr', outpath)
% movefile ('input.pqr', outpath)
% else
%   warningMessage = sprintf('Warning: file does not exist:\n%s', propkaname);
%   uiwait(msgbox(warningMessage));
% end
delete('MATLAB_Solution.dx')
delete('*.dx');
delete('*.inm');
delete('*.log');
wwq=ls;

%clear all
clearvars -except enerT outpath
end
