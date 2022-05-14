function varargout = GUI_Beginner(varargin)
% GUI_Beginner M-file for GUI_Beginner.fig
%      GUI_Beginner, by itself, creates a new GUI_Beginner or raises the existing
%      singleton*.


global enerT

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_Beginner_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_Beginner_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before GUI_Beginner is made visible.
function GUI_Beginner_OpeningFcn(hObject, eventdata, handles, varargin)

global enerT exitmpbec
set(handles.figure1, 'Name', 'Matlab Program for Biomolecular Electrostatic Calculations');
warning('off')
%  default command line output for GUI_Beginner
handles.output = hObject;
handles.dime1=0;
handles.dime2=0;
handles.dime3=0;
handles.glen1=0;
handles.glen2=0;
handles.glen3=0;
handles.cglen1=0;
handles.cglen2=0;
handles.cglen3=0;
handles.grid='Automatic';
handles.temp=298.15;
handles.solvent=78.54;
handles.dielx='';
handles.diely='';
handles.dielz='';
handles.pathdielx='';
handles.kappa='';
handles.charge='';
handles.mol1='';
handles.mol2='';
handles.mol3='';
handles.mol4='';
%handles.cout='';
handles.memorypath='.';
handles.epe='enerepeno';
handles.sfe='enersfeno';
handles.bin='enerbinno';
handles.pka='enerpkano';
handles.bc='mdh';
handles.digit=6;
handles.filename='inputfile.inm';
handles.filename2='focusname.inm';
handles.chgm='SPL2';
handles.srfm='MOL';
handles.solute=2;
handles.surfdens=10;
handles.inpufiles='';
handles.pdb2pqr='PQR';
handles.swin=0.3;
handles.srad=1.4;
handles.ion.charge=zeros(4,1);
handles.ion.charge(1)=1.0;
handles.ion.charge(2)=-1.0;
handles.ion.charge(3)=2.0;
handles.ion.charge(4)=-2.0;
handles.ion.conc=zeros(4,1);
handles.ion.conc(1)=0.15;
handles.ion.conc(2)=0.15;
handles.ion.conc(3)=0.0;
handles.ion.conc(4)=0.0;
handles.ion.radius=zeros(4,1);
handles.ion.radius(1)=2.0;
handles.ion.radius(2)=1.8;
handles.ion.radius(3)=2.0;
handles.ion.radius(4)=1.8;
handles.writeout=zeros(5,1);
handles.writeout(4)=1.0;
handles.writeout(5)=1.0;
handles.check=0;
handles.selectsolver='biconjgrad';
handles.focus='';
handles.ph=-2;
handles.complexsyst='';
handles.pdbmol1='';
handles.quit = 0;
% Update handles structure
cd ..
cd examples
cd solutions
handles.cout = pwd;
cd ..
cd ..
cd Main
set( findall(handles.uipanel7, '-property', 'Enable'), 'Enable', 'off');
set( findall(handles.edit7, '-property', 'Enable'), 'Enable', 'off');
set( findall(handles.uipanel13, '-property', 'Enable'), 'Enable', 'off');
set( findall(handles.uipanel9, '-property', 'Enable'), 'Enable', 'off');
set( findall(handles.uipanel10, '-property', 'Enable'), 'Enable', 'off');
set( findall(handles.uipanel37, '-property', 'Enable'), 'Enable', 'off');
set( findall(handles.edit19, '-property', 'Enable'), 'Enable', 'off');
set( findall(handles.text32, '-property', 'Enable'), 'Enable', 'off');
set( findall(handles.edit25, '-property', 'Enable'), 'Enable', 'off');
set( findall(handles.uipanel15, '-property', 'Enable'), 'Enable', 'off');
set( findall(handles.uipanel21, '-property', 'Enable'), 'Enable', 'off');
set( findall(handles.edit10, '-property', 'Enable'), 'Enable', 'off');
set( findall(handles.edit23, '-property', 'Enable'), 'Enable', 'off');
set( findall(handles.edit52, '-property', 'Enable'), 'Enable', 'off');
set( findall(handles.edit53, '-property', 'Enable'), 'Enable', 'off');
set( findall(handles.edit54, '-property', 'Enable'), 'Enable', 'off');
set( findall(handles.edit55, '-property', 'Enable'), 'Enable', 'off');
set( findall(handles.edit57, '-property', 'Enable'), 'Enable', 'off');
set( findall(handles.edit58, '-property', 'Enable'), 'Enable', 'off');
set( findall(handles.edit59, '-property', 'Enable'), 'Enable', 'off');
set( findall(handles.edit60, '-property', 'Enable'), 'Enable', 'off');
set( findall(handles.uipanel3, '-property', 'Enable'), 'Enable', 'off');
set( findall(handles.radiobutton98, '-property', 'Enable'), 'Enable', 'off');

guidata(hObject, handles);

function varargout = GUI_Beginner_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;

% read dielectric map along x 
function edit7_Callback(hObject, eventdata, handles)

if ((strcmp(get(hObject,'String'),'')==0 && strcmp(handles.inputfiles,'NO')==1))
    msgbox('To "upload the map files" please click "yes" on the above button and try again. Thanks.','!!WARNING!!');
else
filespec = {'*.dx', 'Open Format files (*.dx)'};
[handles.dielx, handles.pathdielx] = uigetfile(filespec, 'Select the Required File',handles.memorypath);

if (ischar(handles.dielx)==0)
    handles.dielx='';
elseif strcmp(handles.dielx, '')==0

handles.memorypath=handles.pathdielx;
end
set(hObject, 'String',handles.dielx);
guidata(hObject,handles);
end

function edit7_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% read temperature
function edit10_Callback(hObject, eventdata, handles)

handles.temp=str2double(get(hObject,'String'));
guidata(hObject, handles);

function edit10_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% read dielectric map along y
function edit12_Callback(hObject, eventdata, handles)

if ((strcmp(get(hObject,'String'),'')==0 && strcmp(handles.inputfiles,'NO')==1))
    msgbox('To "upload the map files" please click "yes" on the above button and try again. Thanks.','!!WARNING!!');
else
filespec = {'*.dx', 'Open Format files (*.dx)'};
[handles.diely, handles.pathdielx] = uigetfile(filespec, 'Select the Required File',handles.memorypath);

if (ischar(handles.diely)==0)
    handles.diely='';
elseif strcmp(handles.diely, '')==0

handles.memorypath=handles.pathdielx;
end
set(hObject, 'String',handles.diely);
guidata(hObject,handles);

end

function edit12_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% read dielectric map along z
function edit13_Callback(hObject, eventdata, handles)

if ((strcmp(get(hObject,'String'),'')==0 && strcmp(handles.inputfiles,'NO')==1))
    msgbox('To "upload the map files" please click "yes" on the above button and try again. Thanks.','!!WARNING!!');
else
filespec = {'*.dx', 'Open Format files (*.dx)'};
[handles.dielz, handles.pathdielx] = uigetfile(filespec, 'Select the Required File',handles.memorypath);

if (ischar(handles.dielz)==0)
    handles.dielz='';
elseif strcmp(handles.dielz, '')==0

handles.memorypath=handles.pathdielx;
end
set(hObject, 'String',handles.dielz);
guidata(hObject,handles);

end

function edit13_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% read kapp map
function edit14_Callback(hObject, eventdata, handles)

if ((strcmp(get(hObject,'String'),'')==0 && strcmp(handles.inputfiles,'NO')==1))
    msgbox('To "upload the map files" please click "yes" on the above button and try again. Thanks.','!!WARNING!!');
else
filespec = {'*.dx', 'Open Format files (*.dx)'};
[handles.kappa, handles.pathdielx] = uigetfile(filespec, 'Select the Required File',handles.memorypath);

if (ischar(handles.kappa)==0)
    handles.kappa='';
elseif strcmp(handles.kappa, '')==0

handles.memorypath=handles.pathdielx;
end
set(hObject, 'String',handles.kappa);
guidata(hObject,handles);

end

function edit14_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% read molecular structure molecule 1
function edit15_Callback(hObject, eventdata, handles)
if strcmp(handles.pdb2pqr, 'PQR')==1
filespec = {'*.pqr', 'Format File (*.pqr)'};
cd ..
cd examples
cd pqr
[handles.mol1, handles.pathdielx] = uigetfile(filespec, 'Select the Required File',handles.memorypath);
cd ..
cd ..
cd Main
end
if strcmp(handles.pdb2pqr, 'PDB')==1
    filespec = {'*.pdb', 'Format File (*.pdb)'};
cd ..
cd examples
cd pdb
[handles.mol1, handles.pathdielx] = uigetfile(filespec, 'Select the Required File');
if (ischar(handles.mol1)==0)
m1 = msgbox('Please select a file. Thanks.','!!WARNING!!'); 
 waitfor(m1)
 cd ..
cd ..
cd Main
else
cd ..
cd ..
% cd Main
% cd ..
cd pdb2pqr
h = msgbox('This could be a minute. Patience, MPBEC is working for you...','Converting pdb to pqr','help');
waitfor(h)
str=computer;
if (strcmp(str,'PCWIN')==1 ||strcmp(str,'PCWIN64')==1)
    cd pdb2pqrwin
    path1 = pwd;
    execut=strcat(path1,'\pdb2pqr');
else
    cd pdb2pqrlinux
    path1 = pwd;
    execut=strcat(path1,'/pdb2pqr');
end

cd ..
cd ..
cd Main
path6=strcat(handles.pathdielx,handles.mol1);
outputpqrfilename=regexprep(handles.mol1,'.pdb','.pqr');
path4=strcat(handles.pathdielx,outputpqrfilename);
warning('off','all')
file_id = fopen('pdb2pqrscript.st', 'rt');
forcefield=fgetl(file_id);
fclose(file_id);
argum=sprintf('%s --ff=%s %s %s',execut,forcefield, path6, path4);
status=system(argum);
handles.pdbmol1=handles.mol1;
handles.mol1=outputpqrfilename;
warning('on','all')
clc
end

end

if (ischar(handles.mol1)==0)
    handles.mol1='';
elseif strcmp(handles.mol1, '')==0

handles.memorypath=handles.pathdielx;
end
set(hObject, 'String',handles.mol1);
guidata(hObject,handles);

function edit15_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% read molecular structure molecule reference system
function edit16_Callback(hObject, eventdata, handles)
if strcmp(handles.pdb2pqr, 'PQR')==1
filespec = {'*.pqr', 'Format File (*.pqr)'};
[handles.mol4, handles.pathdielx] = uigetfile(filespec, 'Select the Required File',handles.memorypath);
else
    filespec = {'*.pdb', 'Format File (*.pdb)'};
[handles.mol4, handles.pathdielx] = uigetfile(filespec, 'Select the Required File',handles.memorypath);
if (ischar(handles.mol4)==0)
m1 = msgbox('Please select a file. Thanks.','!!WARNING!!'); 
 waitfor(m1)
 cd ..
cd ..
cd Main
else
cd ..
cd pdb2pqr
h = msgbox('This could be a minute. Patience, MPBEC is working for you...','Converting pdb to pqr','help');
waitfor(h)
str=computer;
if (strcmp(str,'PCWIN')==1 ||strcmp(str,'PCWIN64')==1)
    cd pdb2pqrwin
    path1 = pwd;
    execut=strcat(path1,'\pdb2pqr');
else
    cd pdb2pqrlinux
    path1 = pwd;
    execut=strcat(path1,'/pdb2pqr');
end

cd ..
cd ..
cd Main
path6=strcat(handles.pathdielx,handles.mol4);
outputpqrfilename=regexprep(handles.mol4,'.pdb','.pqr');
path4=strcat(handles.pathdielx,outputpqrfilename);
warning('off','all')
file_id = fopen('pdb2pqrscript.st', 'rt');
forcefield=fgetl(file_id);
fclose(file_id);
argum=sprintf('%s --ff=%s %s %s',execut,forcefield, path6, path4);
status=system(argum);
handles.mol4=outputpqrfilename;
warning('on','all')
clc
end
end
if (ischar(handles.mol4)==0)
    handles.mol4='';
elseif strcmp(handles.mol4, '')==0

handles.memorypath=handles.pathdielx;
end
set(hObject, 'String',handles.mol4);
guidata(hObject,handles);

function edit16_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% read digits of presicion linea solver
function edit17_Callback(hObject, eventdata, handles)

handles.digit=int8(str2double(get(hObject,'String')));
guidata(hObject, handles);

function edit17_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% solvent dielectric constant reading
function edit19_Callback(hObject, eventdata, handles)

handles.solvent=str2double(get(hObject,'String'));
guidata(hObject, handles);

function edit19_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% read output directory
function edit23_Callback(hObject, eventdata, handles)
cd ..
cd examples
cd solutions
handles.cout = uigetdir('please, select Directory to Open');
if (ischar(handles.cout)==0)
    handles.cout='';
end
set(hObject, 'String',handles.cout);
cd ..
cd ..
cd Main
guidata(hObject,handles);

% ????
function edit23_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Dirichlet boundary condition single charges (sdh)
function radiobutton1_Callback(hObject, eventdata, handles)

if (get(hObject,'Value') == get(hObject,'Max'))
handles.bc='sdh';
guidata(hObject, handles);
end

% need PDB files ? (open web browsser)
function pushbutton3_Callback(hObject, eventdata, handles)

url='http://www.rcsb.org';
 status1=0;
 status2=0;
 status6=0;
 status=0;
 status7=0;
 str=computer;
 if (strcmp(str,'PCWIN')==1 ||strcmp(str,'PCWIN64')==1)
     status6 =  dos('start http://www.rcsb.org');
      if status6~=0
     msgbox('Sorry. We cannot find a web browser application on your computer.','Sorry','none');
      end
 end
 if (strcmp(str,'GLNX86')==1 ||strcmp(str,'GLNXA64')==1)
 status = system('firefox http://www.rcsb.org');
 if status~=0
    status1 = system('opera http://www.rcsb.org');
 else
     status1=0;
 end
 if status1~=0
     status2 = system('google-chrome http://www.rcsb.org');
 else
     status2=0;
 end
     if status2~=0
     msgbox('Sorry. We cannot find a web browser application on your computer.','Sorry','none');
     end
 end
 if (strcmp(str,'MACI')==1 ||strcmp(str,'MACI64')==1)
     status7 = web(url,'-browser');
     if status7~=0
     msgbox('Sorry. We cannot find a web browser application on your computer.','Sorry','none');
     end
 end

 % ?????
function radiobutton5_Callback(hObject, eventdata, handles)

if (get(hObject,'Value') == get(hObject,'Max'))
handles.bc='periodic';
guidata(hObject, handles);
end

% ???
function radiobutton6_Callback(hObject, eventdata, handles)

if (get(hObject,'Value') == get(hObject,'Max'))
handles.bc='mixed';
guidata(hObject, handles);
end

% solute dielectric constant reading
function edit25_Callback(hObject, eventdata, handles)

handles.solute=str2double(get(hObject,'String'));
guidata(hObject, handles);

function edit25_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% read Surface density parameter in dielectric map approx
function edit26_Callback(hObject, eventdata, handles)

if ((strcmp(get(hObject,'String'),'')==0 && strcmp(handles.inputfiles,'YES')==1))
set(hObject,'String','');    
msgbox('This parameter will be ignored because you selected the option "upload the map files" instead. To assign a value to this parameter, please click "No" in the uploading map panel. Thanks.','!!WARNING!!');
else
handles.surfdens=str2double(get(hObject,'String'));
guidata(hObject, handles);    
end

function edit26_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% upload map files ? NO
function radiobutton8_Callback(hObject, eventdata, handles)

 if (get(hObject,'Value') == get(hObject,'Max'))

handles.inputfiles='NO';
guidata(hObject, handles);

 end
  if (get(hObject,'Value') == get(hObject,'Min'))

handles.inputfiles='YES';
guidata(hObject, handles);

  end

  % upload map files ? YES
function radiobutton9_Callback(hObject, eventdata, handles)

if (get(hObject,'Value') == get(hObject,'Max'))

handles.inputfiles='YES';
guidata(hObject, handles);

end
if (get(hObject,'Value') == get(hObject,'Min'))

handles.inputfiles='NO';
guidata(hObject, handles);

end
 
 % SPLO charge map approximation
function radiobutton16_Callback(hObject, eventdata, handles)


if (get(hObject,'Value') == get(hObject,'Max'))
handles.chgm='SPL0';
guidata(hObject, handles);
if (strcmp (handles.inputfiles, 'YES')==1 )
msgbox('This parameter will be ignored because you selected the option "upload the map files" instead. In order to use this method please click "No" in the uploading map panel. Thanks.','!!WARNING!!');
end
end

% SPL2 charge map approximation
function radiobutton17_Callback(hObject, eventdata, handles)


if (get(hObject,'Value') == get(hObject,'Max'))
handles.chgm='SPL2';
guidata(hObject, handles);
if (strcmp (handles.inputfiles, 'YES')==1 )
msgbox('This parameter will be ignored because you selected the option "upload the map files" instead. In order to use this method please click "No" in the uploading map panel. Thanks.','!!WARNING!!');
end
end

% SPL2 dielectric map approx
function radiobutton13_Callback(hObject, eventdata, handles)

if (get(hObject,'Value') == get(hObject,'Max'))
handles.srfm='SPL2';
guidata(hObject, handles);
if (strcmp (handles.inputfiles, 'YES')==1 )
msgbox('This parameter will be ignored because you selected the option "upload the map files" instead. In order to use this method please click "No" in the uploading map panel. Thanks.','!!WARNING!!');
end
end

% SMOL dielctric map approx
function radiobutton15_Callback(hObject, eventdata, handles)

if (get(hObject,'Value') == get(hObject,'Max'))
handles.srfm='SMOL';
guidata(hObject, handles);
if (strcmp (handles.inputfiles, 'YES')==1 )
msgbox('This parameter will be ignored because you selected the option "upload the map files" instead. In order to use this method please click "No" in the uploading map panel. Thanks.','!!WARNING!!');
end
end

% MOL dielectric map approx
function radiobutton12_Callback(hObject, eventdata, handles)

if (get(hObject,'Value') == get(hObject,'Max'))
handles.srfm='MOL';
guidata(hObject, handles);
if (strcmp (handles.inputfiles, 'YES')==1 )
msgbox('This parameter will be ignored because you selected the option "upload the map files" instead. In order to use this method please click "No" in the uploading map panel. Thanks.','!!WARNING!!');
end
end

% read charge map
function edit27_Callback(hObject, eventdata, handles)

if ((strcmp(get(hObject,'String'),'')==0 && strcmp(handles.inputfiles,'NO')==1))
    msgbox('To "upload the map files" please click "yes" in the left icon and then try it again. Thanks.','!!WARNING!!');
else
filespec = {'*.dx', 'Open Format files (*.dx)'};
[handles.charge, handles.pathdielx] = uigetfile(filespec, 'Select the Required File',handles.memorypath);

if (ischar(handles.charge)==0)
    handles.charge='';
elseif strcmp(handles.charge, '')==0

handles.memorypath=handles.pathdielx;
end
set(hObject, 'String',handles.charge);
guidata(hObject,handles);

end

function edit27_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function text32_CreateFcn(hObject, eventdata, handles)

s = sprintf('shifted on x-direction Dielectric function map file in dx format; this causes the pdie, sdie, srad, swin, and srfm parameters \n and the radii of the biomolecular atoms to be ignored when computing dielectric maps for the Poisson-Boltzmann equation. \n Note that the pdie and sdie values are still used for some boundary condition calculations. ');
set(hObject,'TooltipString',s)

function text20_CreateFcn(hObject, eventdata, handles)

s = sprintf('Specify the dielectric constant of the solvent. Bulk water at biologically-relevant \n temperatures is usually modeled with a dielectric constant of 78-80.');
set(hObject,'TooltipString',s)


function text28_CreateFcn(hObject, eventdata, handles)

s = sprintf('Specify the dielectric constant of the biomolecule. This is usually a value between 2 to 20, where\n  lower values consider only electronic polarization and higher values \n consider additional polarization due to intramolecular motion.');
set(hObject,'TooltipString',s)


function radiobutton17_CreateFcn(hObject, eventdata, handles)

s2 = sprintf('Cubic B-spline interpolation (basis spline). This option results in the location of the charge in \n the nearest and next-nearest grid points. Whereas this calculation requires \n higher computational cost than SPLO, it provides softer distribution of charges, \n which in turn generate mean electrostatic potential solutions less sensitive \n to the “Grid information�? parameters (number of grid points and box size). ');
set(hObject,'TooltipString',s2)

function radiobutton16_CreateFcn(hObject, eventdata, handles)

s1 = sprintf('Trilinear interpolation (linear splines). This method locates the charge of each atom in the nearest neighbor grid point. \n The resulting charge density map may not be a smooth function, having sharp corners at the data points. \n SPLO has been shown to be an effcient algorithm to generate charge density maps but the accuracy of the PB solution \n may depend on the mesh, the box size and the protein structure. ');
set(hObject,'TooltipString',s1)


function radiobutton13_CreateFcn(hObject, eventdata, handles)

s3 = sprintf('More suitable approaches (compared to MOL amd SMOL) are required in those calculations demanding the \n evaluation of the slope of the mean electrostatic potential (e.g. electric force). In this case, \n the dielectric and kappa maps are built up on a cubic-spline surface where an intermediate dielectric region \n is introduced at the solute-solvent boundary to smooth the transition. The parameter "Swin" specifies the rate of change \n (e.g. the window width) in the dielectric interface. The preassigned value for this parameter is 0.3Å.');
set(hObject,'TooltipString',s3)


function radiobutton15_CreateFcn(hObject, eventdata, handles)

s = sprintf('This option includes a harmonic interpolation that smooths the surface on the dielectric and kappa maps \n generated by MOL, which reduces the sensitivity to "Grid information" parameters" \n at the expense of increasing the computational cost.');
set(hObject,'TooltipString',s)

function radiobutton12_CreateFcn(hObject, eventdata, handles)

s = sprintf('The dielectric coefficients are assigned according to the definition of the solvent accessible surface. \n This surface is calculated using the Collony molecular surface approach. A solvent sphere is used \n to probe the molecule shape and create a surface defined by the union of the solvent-sized spheres, equidistant to the surface atoms in \n the molecule. The radius of the solvent probe "srad" is preassigned with the value 1.4Å, which represents the radius of a water molecule.\n The volume outside the solvent-accessible surface is assigned with "solvent dielectric coefficients" \n whereas the rest of the space is assigned with the "solute dielectric coeffcients". Setting "srad" equal to zero \n generates a typical van der Waals surface. The ion-accessibility map (kappa map) is built up on a modified van der Waals surface. \n The radius of each atom in the molecule is inflated with the radius of the ion species (See "Ionic species" panel). \n The space inside the modified surface is assigned with an ion-accessibility coefficient equal to zero whereas the volume outside \n of this surface is assigned with the corresponding kappa values. This method is very efficient, mainly for large proteins. However, the \n abrupt transition of the dielectric and kappa maps at the solute-solvent boundary may give rise to discontinuous \n variations of the electrostatic free energies.');
set(hObject,'TooltipString',s)


function radiobutton2_CreateFcn(hObject, eventdata, handles)

s = sprintf('"Focusing" boundary condition. Dirichlet condition where the potential \n at the boundary is set to the values computed by the previous (usually lower-resolution) PB calculation.\n This is used in sequential focusing performed manually in mg-manual calculations. \n All of the boundary points should lie within the domain of the previous calculation for best accuracy; \n if any boundary points lie outside, their values are computed using single \n Debye-Hückel boundary conditions (see above).');
set(hObject,'TooltipString',s)


function radiobutton1_CreateFcn(hObject, eventdata, handles)

s = sprintf('"Single Debye-Hückel" boundary condition. Dirichlet condition where the potential at the\n  boundary is set to the values prescribed by a Debye-Hückel model for a multiple,\n  non-interacting spheres with a point charges. The sphere radii are set to the atomic radii \n of the biomolecule and the sphere charges are set to the total charge of the protein. \n This condition works better than mdh for large biomolecules.');
set(hObject,'TooltipString',s)


function radiobutton8_CreateFcn(hObject, eventdata, handles)

s = sprintf('Answer "No" if you dont want to upload these files. Instead, MPBEC will calculate the corresponding maps. \n In this case you have to determine "chgm" and "sfrm". \n Remember, you always have to upload the pqr files.');
set(hObject,'TooltipString',s)

if (get(hObject,'Value') == get(hObject,'Max'))
handles.inputfiles='NO';
guidata(hObject, handles);

 end
  if (get(hObject,'Value') == get(hObject,'Min'))

handles.inputfiles='YES';
guidata(hObject, handles);

  end
 

function radiobutton9_CreateFcn(hObject, eventdata, handles)

s = sprintf('Answer "Yes" if you want to upload these files. In this case, you are not able to choose the "chgm" and "sfrm" methods. \n Remember, you always have to upload the pqr files.');
set(hObject,'TooltipString',s)
if (get(hObject,'Value') == get(hObject,'Max'))

handles.inputfiles='YES';
guidata(hObject, handles);

 end
  if (get(hObject,'Value') == get(hObject,'Min'))

handles.inputfiles='NO';
guidata(hObject, handles);

 end

% read Swin parameter in dielectric map approx
function edit28_Callback(hObject, eventdata, handles)

if ((strcmp(get(hObject,'String'),'')==0 && strcmp(handles.inputfiles,'YES')==1))
        set(hObject,'String','');
    msgbox('This parameter will be ignored because you selected the option "upload the map files" instead. To assign a value to this parameter, please click "No" in the uploading map panel. Thanks.','!!WARNING!!');
else
handles.swin=str2double(get(hObject,'String'));
guidata(hObject, handles);    
end


function edit28_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function text29_CreateFcn(hObject, eventdata, handles)

s = sprintf('Please provide the minimum number of points per area unit. Usually, it is equal to 10. Thanks.');
set(hObject,'TooltipString',s)


function text38_CreateFcn(hObject, eventdata, handles)

s = sprintf('This parameter will be used by the SPL2 method only. Otherwise, it will be ignored. \n It is a floating point number (usually 0.3 A) for the spline window width. For more details please \n read the reference Nina, Im, and Roux (doi:10.1016/S0301-4622(98)00236-1).');
set(hObject,'TooltipString',s)

% ????
function edit26_ButtonDownFcn(hObject, eventdata, handles)

% ????
function edit28_ButtonDownFcn(hObject, eventdata, handles)

% ????
function radiobutton5_CreateFcn(hObject, eventdata, handles)

s = sprintf('periodic boundary conditions over the xy plane and non periodic (Dirichlet) boundary condition along the z direction.');
set(hObject,'TooltipString',s)

% read Srad parameter in dielectric map approx
function edit29_Callback(hObject, eventdata, handles)

if ((strcmp(get(hObject,'String'),'')==0 && strcmp(handles.inputfiles,'YES')==1))
        set(hObject,'String','');
    msgbox('This parameter will be ignored because you selected the option "upload the map files" instead. To assign a value to this parameter, please click "No" in the uploading map panel. Thanks.','!!WARNING!!');
else
handles.srad=str2double(get(hObject,'String'));
guidata(hObject, handles);    
end


function edit29_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit29_ButtonDownFcn(hObject, eventdata, handles)



function text39_CreateFcn(hObject, eventdata, handles)

s = sprintf('Specify the radius of the solvent molecules; this parameter is used \n to define the dielectric function for probe-based dielectric definitions. \n This value is usually set to 1.4 Å for water. This keyword is ignored when any \n of the spline-based surfaces are used (e.g., SPL2), since they are not probe-based.');
set(hObject,'TooltipString',s)


% Dirichlet boundary condition multiple charges
function radiobutton20_Callback(hObject, eventdata, handles)

if (get(hObject,'Value') == get(hObject,'Max'))
handles.bc='mdh';
guidata(hObject, handles);
end

% read ion charge specie 1
function edit30_Callback(hObject, eventdata, handles)

handles.ion.charge(1)=str2double(get(hObject,'String'));
guidata(hObject, handles);

function edit30_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% read ion density specie 1
function edit31_Callback(hObject, eventdata, handles)

handles.ion.conc(1)=str2double(get(hObject,'String'));
guidata(hObject, handles);

function edit31_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% read ion radius specie 1
function edit32_Callback(hObject, eventdata, handles)

handles.ion.radius(1)=str2double(get(hObject,'String'));
guidata(hObject, handles);

function edit32_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% read ion charge specie 2
function edit33_Callback(hObject, eventdata, handles)

handles.ion.charge(2)=str2double(get(hObject,'String'));
guidata(hObject, handles);

function edit33_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% read ion density specie 2
function edit34_Callback(hObject, eventdata, handles)

handles.ion.conc(2)=str2double(get(hObject,'String'));
guidata(hObject, handles);

function edit34_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% read ion radius specie 2
function edit35_Callback(hObject, eventdata, handles)

handles.ion.radius(2)=str2double(get(hObject,'String'));
guidata(hObject, handles);

function edit35_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% read ion charge specie 3
function edit36_Callback(hObject, eventdata, handles)

handles.ion.charge(3)=str2double(get(hObject,'String'));
guidata(hObject, handles);

function edit36_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% read ion density specie 3
function edit37_Callback(hObject, eventdata, handles)

handles.ion,conc(3)=str2double(get(hObject,'String'));
guidata(hObject, handles);

function edit37_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% read ion radius specie 3
function edit38_Callback(hObject, eventdata, handles)

handles.ion.radius(3)=str2double(get(hObject,'String'));
guidata(hObject, handles);

function edit38_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% read ion charge specie 4
function edit39_Callback(hObject, eventdata, handles)

handles.ion.charge(4)=str2double(get(hObject,'String'));
guidata(hObject, handles);


function edit39_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% read ion density specie 4
function edit40_Callback(hObject, eventdata, handles)

handles.ion.conc(4)=str2double(get(hObject,'String'));
guidata(hObject, handles);

function edit40_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% read ion radius specie 4
function edit41_Callback(hObject, eventdata, handles)

handles.ion.radius(4)=str2double(get(hObject,'String'));
guidata(hObject, handles);

function edit41_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function radiobutton20_CreateFcn(hObject, eventdata, handles)

s = sprintf('"Multiple Debye-Hückel" boundary condition. Dirichlet condition where the potential at the\n  boundary is set to the values prescribed by a Debye-Hückel model for a multiple,\n  non-interacting spheres with a point charges. The sphere radii are set to the atomic radii \n of the biomolecule and the sphere charges are set to the total charge of the protein. \n This condition works better than sdh for closer boundaries but can be very slow for large biomolecules.');
set(hObject,'TooltipString',s)

function radiobutton21_CreateFcn(hObject, eventdata, handles)
s = sprintf('Biconjugated gradient stabilized method is used. If it fails, gmres is used then.');
set(hObject,'TooltipString',s)

% write uot charge map ?
function checkbox1_Callback(hObject, eventdata, handles)

if (get(hObject,'Value') == get(hObject,'Max'))
handles.writeout(1)=1.0;
else
  handles.writeout(1)=0.0;  
end
guidata(hObject, handles);

% write out Kappa map ?
function checkbox2_Callback(hObject, eventdata, handles)

if (get(hObject,'Value') == get(hObject,'Max'))
handles.writeout(2)=1.0;
else
  handles.writeout(2)=0.0;  
end
guidata(hObject, handles);

% write out dielectric maps ?
function checkbox3_Callback(hObject, eventdata, handles)

if (get(hObject,'Value') == get(hObject,'Max'))
handles.writeout(3)=1.0;
else
    handles.writeout(3)=0.0;
end
guidata(hObject, handles);

% write out MEP map ?
function checkbox4_Callback(hObject, eventdata, handles)

if (get(hObject,'Value') == get(hObject,'Max'))
handles.writeout(4)=1.0;
else
 handles.writeout(4)=0.0;   
end
guidata(hObject, handles);

% plot MEP on the molecular surface and isosurface ? 
function checkbox5_Callback(hObject, eventdata, handles)

if (get(hObject,'Value') == get(hObject,'Max'))
handles.writeout(5)=1.0;
else
    handles.writeout(5)=0.0;
end
guidata(hObject, handles);

% COMPUTE SOLUTION (red button)
function pushbutton2_Callback(hObject, eventdata, handles)

global enerT outpath

cd ..
ese=pwd;
addpath(genpath(ese));
cd Main

jFrame = get(handle(GUI_Beginner),'JavaFrame');
jFrame.setMinimized(true);
pause(1);
   if (exist(strcat('solv_',handles.filename),'file') == 2)
       argg=sprintf('%s %s %s %s %s','The input file name',handles.filename,'already exists.To prevent overwritting data, we will rename the file as',strcat('new',handles.filename),'Thanks! ');
       handles.filename=strcat('new',handles.filename);
       m2=  msgbox(argg,'Confirmation','warn','modal');  
     waitfor(m2)
   end
ionicstrenght=0.0;
for i=1:4
ionicstrenght=ionicstrenght+0.5*handles.ion.conc(i)*handles.ion.charge(i)^2;
end

if strcmp(handles.mol2, '')==1
    handles.mol4=handles.mol1;
end
% automatic grid calculation

           
  if (strcmp(handles.mol1, '')==0 && strcmp(handles.mol3, '')==1) 
           path4=strcat(handles.pathdielx,handles.mol1);
  end
  if ( strcmp(handles.mol3, '')==0)
           path4=strcat(handles.pathdielx,handles.mol3);
  end

file_id = fopen(path4, 'rt');
flag=[];
headerlinesIn=-1;
while isempty(flag)
    line=fgetl(file_id);
    flag=strfind(line,'ATOM');
    headerlinesIn=headerlinesIn+1;
end
fclose(file_id);

delimiterIn = ' ';

PQRdata=importdata(path4,delimiterIn, headerlinesIn); 

atomdata=PQRdata.data;
atomPx=atomdata(:,2);
atomPy=atomdata(:,3);
atomPz=atomdata(:,4);
lx=abs(max(atomPx)-min(atomPx));
ly=abs(max(atomPy)-min(atomPy));
lz=abs(max(atomPz)-min(atomPz));

  if ((strcmp(handles.grid, 'Automatic')==1))

dh=[lx ly lz]/255.0*3/2;
dhm=max(dh);
if (dhm>0.5)
    handles.cglen1=lx*3/2;
    handles.cglen2=ly*3/2;
    handles.cglen3=lz*3/2;
    handles.dime1=1+int32(lx*3/2/dhm);
    handles.dime2=1+int32(ly*3/2/dhm);
    handles.dime3=1+int32(lz*3/2/dhm);
    handles.focus='focusname.inm';
    handles.glen1=lx*5/4;
    handles.glen2=ly*5/4;
    handles.glen3=lz*5/4;
else
    handles.glen1=lx*3/2;
    handles.glen2=ly*3/2;
    handles.glen3=lz*3/2;
    if (handles.glen1<20)
        handles.glen1=19.5;
    end
    if (handles.glen2<20)
        handles.glen2=19.5;
    end
    if (handles.glen3<20)
        handles.glen3=19.5;
    end
    handles.dime1=1+int32(handles.glen1/0.5);
    handles.dime2=1+int32(handles.glen2/0.5);
    handles.dime3=1+int32(handles.glen2/0.5);
    handles.focus='';
end
  end
% end automatic grid calculation

% vector definition

maxionR=max(handles.ion.radius);
dime(1)=handles.dime1;
dime(2)=handles.dime2;
dime(3)=handles.dime3;
glen(1)=handles.glen1;
glen(2)=handles.glen2;
glen(3)=handles.glen3;
cglen(1)=handles.cglen1;
cglen(2)=handles.cglen2;
cglen(3)=handles.cglen3;
bulk(1)=ionicstrenght;
bulk(2)=maxionR;
bulk(3)=handles.solvent;
bulk(4)=handles.solute;
bulk(5)=handles.temp;
bulk(6)=handles.srad;
bulk(7)=handles.surfdens;
bulk(8)=handles.swin;
molecule=handles.mol1;
moleculeref=handles.mol4;

if strcmp(handles.focus, 'focusname.inm')==1
    boundcond=handles.focus;
else
    boundcond=handles.bc;
end

if (strcmp(handles.epe, 'enerepeyes')==1 || strcmp(handles.sfe, 'enersfeyes')==1 || strcmp(handles.bin, 'enerbinyes')==1 )
    energy='calceneryes';
else
    energy='calcenerno';
end

testnum=(norm(double(dime))+norm(double(glen))+norm(double(bulk))+...
    abs(handles.temp)+abs(handles.digit));

if (strcmp(handles.inputfiles,'NO')==1)
    handles.dielx='nouploadedfiles';
    handles.diely=handles.dielx;
    handles.dielz=handles.dielx;
    handles.kappa=handles.dielx;
    handles.charge=handles.dielx;
    
   qq1=norm(find(handles.ion.conc<0));
   qq2=norm(find(handles.ion.conc>5));
   qqq=qq1+qq2;
   rr1=norm(find(handles.ion.radius<0));
   rr2=norm(find(handles.ion.radius>5)) ;
   rrr=rr1+rr2;
    if (handles.surfdens<=0 || handles.surfdens>=50)
   m1=  msgbox('the surface density must be greater than 0 and smaller than 50. It is often set equal to 10. Thanks!!.','Confirmation','warn','modal');   
    waitfor(m1)
    handles.check=1;
 guidata(hObject, handles);

   elseif (handles.swin<=0 || handles.swin>=3)
   m1=  msgbox('the swin must be greater than 0 and smaller than 3. It is often set equal to 0.3. Thanks!!.','Confirmation','warn','modal');   
    waitfor(m1)
    handles.check=1;
 guidata(hObject, handles);

  elseif (handles.srad<0 || handles.srad>=3)
   m1=  msgbox('the molecular solvent radius must be greater than 0 and smaller than 3. It is often set equal to 1.4 for water. Thanks!!.','Confirmation','warn','modal');
    waitfor(m1)
    handles.check=1;
 guidata(hObject, handles);

  elseif (qqq>0)
  m1=   msgbox('the ionic concentration must be greater than or equal to 0 and smaller than 3 (M). Thanks!!.','Confirmation','warn','modal');  
   waitfor(m1)
   handles.check=1;
 guidata(hObject, handles);

  elseif (rrr>0)
  m1=   msgbox('the ionic radius must be greater than or equal to 0 and smaller than 5 (A). Thanks!!.','Confirmation','warn','modal');   
   waitfor(m1)
   handles.check=1;
 guidata(hObject, handles);
  else
     handles.check=0;
 guidata(hObject, handles);
   end

else
    handles.surfdens=0;
    handles.swin=handles.surfdens;
    handles.srad=handles.surfdens;
    handles.solute=handles.surfdens;
    handles.srfm='uploadedfiles';
    handles.chgm=handles.srfm;    
end

if handles.check==0
    
qa=double(int8(handles.digit)-handles.digit);

   if (isnan(testnum)==1 || testnum==0)
 m2=   msgbox('There is a typo. Please check the value of the numerical parameters. One of them is NaN. Thanks!!.','Confirmation','warn','modal');   
     waitfor(m2)

    elseif (strcmp(handles.mol1,'')==1)
   m2=  msgbox('Please upload the pqr file(s). Thanks!!.','Confirmation','warn','modal');  
    waitfor(m2)

   elseif (strcmp(handles.cout,'')==1)
   m2=  msgbox('Please enter the output files directory. Thanks!!.','Confirmation','warn','modal');   
     waitfor(m2)

   elseif (((strcmp(handles.grid,'Manual')==1) && handles.dime1<=0) || ((strcmp(handles.grid,'Manual')==1) && handles.dime2<=0) ||((strcmp(handles.grid,'Manual')==1) && handles.dime3<=0))
   m2=  msgbox('the number of points must be greater than 0. Thanks!!.','Confirmation','warn','modal');   
     waitfor(m2)
     
   elseif   (((strcmp(handles.focus,'focusname.inm')==1) && handles.glen1<lx) || ((strcmp(handles.focus,'focusname.inm')==1) && handles.glen2<ly) ||((strcmp(handles.focus,'focusname.inm')==1) && handles.glen3<lz))
       if (strcmp(handles.grid,'Manual')==1)
       argg=sprintf('%s %g %g %g %s','the fine box length along x, y, and z must be larger than the molecular size',lx,ly,lz,'Amgstrongs respectively. Please change box size and try again. Thanks! ');
       m2=  msgbox(argg,'Confirmation','warn','modal');  
     waitfor(m2)
       end
      
   elseif (((strcmp(handles.focus,'focusname.inm')==1) && handles.cglen1<3*lx/2) || ((strcmp(handles.focus,'focusname.inm')==1) && handles.cglen2<3*ly/2) ||((strcmp(handles.focus,'focusname.inm')==1) && handles.cglen3<3*lz/2))
      if (strcmp(handles.grid,'Manual')==1)
       argg=sprintf('%s %g %g %g %s','the focus (coarse) box length along x, y, and z must be at least 1.5 times larger than the molecular size',lx,ly,lz,'Amgstrongs respectively. Please change box size and try again. Thanks! ');
       m2=  msgbox(argg,'Confirmation','warn','modal');  
     waitfor(m2)
      end
      
   elseif (((strcmp(handles.focus,'focusname.inm')==0) && handles.glen1<3*lx/2) || ((strcmp(handles.focus,'focusname.inm')==0) && handles.glen2<3*ly/2) ||((strcmp(handles.focus,'focusname.inm')==0) && handles.glen3<3*lz/2))
      if (strcmp(handles.grid,'Manual')==1)
       argg=sprintf('%s %g %g %g %s','the fine box length along x, y, and z must be at least 1.5 times larger than the molecular size',lx,ly,lz,'Amgstrongs respectively. Please change box size and try again. Thanks! ');
       m2=  msgbox(argg,'Confirmation','warn','modal');  
     waitfor(m2)
      end
    
   elseif (handles.solvent<=0 || handles.solute<=0)
   m2=  msgbox('the dielectric coefficient must be greater than 0. Thanks!!.','Confirmation','warn','modal');   
     waitfor(m2)

   elseif ((handles.cglen1<handles.glen1 && strcmp(handles.focus,'focusname.inm')==1))
  m2=   msgbox('the coarse grained box must be greater than the fine (target) box. Thanks!!.','Confirmation','warn','modal');   
     waitfor(m2)

   elseif (handles.temp<=-273 || handles.temp> 600)
  m2=   msgbox('the temperature must be greater than -273 and smaller than 600. Thanks!!.','Confirmation','warn','modal');   
     waitfor(m2)
  

   elseif ((handles.cglen2<handles.glen2 && strcmp(handles.focus,'focusname.inm')==1) )
  m2=   msgbox('the coarse grained box must be greater than the fine (target) box. Thanks!!.','Confirmation','warn','modal');   
     waitfor(m2)
  

   elseif ((handles.cglen3<handles.glen3 && strcmp(handles.focus,'focusname.inm')==1))
  m2=   msgbox('the coarse grained box must be greater than the fine (target) box. Thanks!!.','Confirmation','warn','modal');   
     waitfor(m2)

   elseif (handles.digit<=0 || handles.digit>=14 || qa~=0.0)
  m2=   msgbox('the digits of precision must be an integer greater than 0 and smaller than 14. It is often set equal to 6. Thanks!!.','Confirmation','warn','modal');  
     waitfor(m2)
   elseif (strcmp(handles.mol2, '')==1  && strcmp(handles.bin,'enerbinyes')==1)
    m2= msgbox('To calculate the Binding Energy you must upload the Molecular structures (pqr file) for molecule 2 and complex and try again. Thanks.','!!WARNING!!');
    waitfor(m2)
     elseif (strcmp(handles.mol3, '')==1  && strcmp(handles.bin,'enerbinyes')==1)
    m2= msgbox('To calculate the Binding Energy you must upload the Molecular structures (pqr file) for molecule 2 and complex and try again. Thanks.','!!WARNING!!');
    waitfor(m2)
        elseif (strcmp(handles.mol3, '')==0  && strcmp(handles.mol4,'')==1)
    m2= msgbox('You must upload the Molecular structure reference system and try again. Thanks.','!!WARNING!!');
    waitfor(m2)
    elseif (strcmp(handles.sfe, 'enersfeyes')==1  && strcmp(handles.bin,'enerbinyes')==1)
    m2= msgbox('You cannot calculate both the Solvation Free Energy and Binding Energy simultanelously. Please unclick one of the options and try again. Thanks.','!!WARNING!!');
    waitfor(m2)
   elseif (strcmp(handles.pka, 'enerpkayes')==1 && (handles.ph<=0))
    m2=msgbox('Please assign a value for pH between 1 and 14 to perform the pKa calculations. Thanks.','!!WARNING!!');
     waitfor(m2);
     elseif (strcmp(handles.pka, 'enerpkayes')==1 && (handles.ph>=14))
    m2=msgbox('Please assign a value for pH between 1 and 14 to perform the pKa calculations. Thanks.','!!WARNING!!');
     waitfor(m2);
    elseif (strcmp(handles.pka, 'enerpkayes')==1 && strcmp(handles.pdb2pqr,'PDB')==0)
    m2=msgbox('Please upload the PDB file for molecule 1, not the PQR, to perform the pKa calculations. Thanks.','!!WARNING!!');
     waitfor(m2);
%   elseif (exist(strcat('solv_',handles.filename),'file') == 2)
%   argg=sprintf('%s %s %s','The input file name',handles.filename,'already exists. Please change the name and try again. Thanks! ');
%       m2=  msgbox(argg,'Confirmation','warn','modal'); 
%     waitfor(m2)
%       argg=sprintf('%s %s %s %s %s','The input file name',handles.filename,'already exists.To prevent overwritting data, we will rename the file as',strcat('1',handles.filename),'Thanks! ');
%       handles.filename=strcat('1',handles.filename);
%       m2=  msgbox(argg,'Confirmation','warn','modal');  
%     waitfor(m2)
     
   else
       
 % write out iinput files
 
axx=strcat('solv_',handles.filename);
fid = fopen(axx, 'wt');
fprintf(fid,'%g %g %g\n',dime);
fprintf(fid,'%f %f %f\n', glen);
fprintf(fid,'%f %f %f %f %f %f %f %f\n', bulk);
fprintf(fid,'%s\n',boundcond);
fprintf(fid,'%g %s\n',[handles.digit handles.selectsolver]);
fprintf(fid,'%s\n',handles.dielx);
fprintf(fid,'%s\n',handles.diely);
fprintf(fid,'%s\n',handles.dielz);
fprintf(fid,'%s\n',handles.kappa);
fprintf(fid,'%s\n',handles.charge);
fprintf(fid,'%s\n',molecule);
fprintf(fid,'%s\n',moleculeref);
fprintf(fid,'%s\n',handles.srfm);
fprintf(fid,'%s\n',handles.chgm);
fprintf(fid,'%s\n',energy);
fprintf(fid,'%s\n',handles.pathdielx);
fprintf(fid,'%s\n',handles.cout);
fprintf(fid,'%f %f %f %f %f\n', handles.writeout);
fclose(fid);

     if strcmp(handles.focus, 'focusname.inm')==1

fid = fopen(handles.filename2, 'wt');
fprintf(fid,'%g %g %g\n',dime);
fprintf(fid,'%f %f %f\n', cglen);
fprintf(fid,'%f %f %f %f %f %f %f %f\n', bulk);
fprintf(fid,'%s\n',handles.bc);
fprintf(fid,'%g %s\n',[handles.digit handles.selectsolver]);
fprintf(fid,'%s\n',handles.dielx);
fprintf(fid,'%s\n',handles.diely);
fprintf(fid,'%s\n',handles.dielz);
fprintf(fid,'%s\n',handles.kappa);
fprintf(fid,'%s\n',handles.charge);
fprintf(fid,'%s\n',molecule);
fprintf(fid,'%s\n',moleculeref);
fprintf(fid,'%s\n',handles.srfm);
fprintf(fid,'%s\n',handles.chgm);
fprintf(fid,'%s\n',energy);
fprintf(fid,'%s\n',handles.pathdielx);
fprintf(fid,'%s\n',handles.cout);
fprintf(fid,'%f %f %f %f %f\n', handles.writeout);
fclose(fid);

     end

     handles.check=3;
     guidata(hObject, handles);    

      if (strcmp(handles.pka, 'enerpkayes')==1)
cd ..
cd pdb2pqr
str=computer;
if (strcmp(str,'PCWIN')==1 ||strcmp(str,'PCWIN64')==1)
    cd pdb2pqrwin
    path1 = pwd;
    execut=strcat(path1,'\pdb2pqr');
else
    cd pdb2pqrlinux
    path1 = pwd;
    execut=strcat(path1,'/pdb2pqr');
end
cd ..
cd ..
cd Main
path4=pwd;
path6=strcat(handles.pathdielx,handles.pdbmol1);
outputpqrfilename=regexprep(handles.pdbmol1,'.pdb','.pqr');
path4=strcat(path4,'/',outputpqrfilename);
path2=strcat(handles.pathdielx,'/',outputpqrfilename);
inpqr=strcat(handles.pathdielx,outputpqrfilename);
copyfile (inpqr, 'input.pqr')
warning('off','all')
forcefield='PARSE';
pHflag=handles.ph;
argum=sprintf('%s --ff=%s --with-ph=%f %s %s',execut,forcefield, pHflag,path6, path4);
status=system(argum);
copyfile (path4, 'pKa.pqr')
movefile(outputpqrfilename,path2)
%delete(outputpqrfilename);
warning('on','all')

clc
      end
     
MPBEC(axx)


solvenermol1=enerT;

% sfe calculation

     if (strcmp(handles.sfe, 'enersfeyes')==1)

       bulk(3)=handles.solute;
       axx=strcat('ref_',handles.filename);

    fid = fopen(axx, 'wt');
fprintf(fid,'%g %g %g\n',dime);
fprintf(fid,'%f %f %f\n', glen);
fprintf(fid,'%f %f %f %f %f %f %f %f\n', bulk);
fprintf(fid,'%s\n',boundcond);
fprintf(fid,'%g %s\n',[handles.digit handles.selectsolver]);
fprintf(fid,'%s\n',handles.dielx);
fprintf(fid,'%s\n',handles.diely);
fprintf(fid,'%s\n',handles.dielz);
fprintf(fid,'%s\n',handles.kappa);
fprintf(fid,'%s\n',handles.charge);
fprintf(fid,'%s\n',molecule);
fprintf(fid,'%s\n',moleculeref);
fprintf(fid,'%s\n',handles.srfm);
fprintf(fid,'%s\n',handles.chgm);
fprintf(fid,'%s\n',energy);
fprintf(fid,'%s\n',handles.pathdielx);
fprintf(fid,'%s\n',handles.cout);
fprintf(fid,'%f %f %f %f %f\n', handles.writeout);
fclose(fid);

        if strcmp(handles.focus, 'focusname.inm')==1

fid = fopen(handles.filename2, 'wt');
fprintf(fid,'%g %g %g\n',dime);
fprintf(fid,'%f %f %f\n', cglen);
fprintf(fid,'%f %f %f %f %f %f %f %f\n', bulk);
fprintf(fid,'%s\n',handles.bc);
fprintf(fid,'%g %s\n',[handles.digit handles.selectsolver]);
fprintf(fid,'%s\n',handles.dielx);
fprintf(fid,'%s\n',handles.diely);
fprintf(fid,'%s\n',handles.dielz);
fprintf(fid,'%s\n',handles.kappa);
fprintf(fid,'%s\n',handles.charge);
fprintf(fid,'%s\n',molecule);
fprintf(fid,'%s\n',moleculeref);
fprintf(fid,'%s\n',handles.srfm);
fprintf(fid,'%s\n',handles.chgm);
fprintf(fid,'%s\n',energy);
fprintf(fid,'%s\n',handles.pathdielx);
fprintf(fid,'%s\n',handles.cout);
fprintf(fid,'%f %f %f %f %f\n', handles.writeout);
fclose(fid);

        end


     handles.check=3;
     guidata(hObject, handles);    
 
MPBEC(axx)
refenermol1=enerT;
sfe = solvenermol1-refenermol1;
disp(' ')
disp(['Solvation Free energy = ', num2str(sfe)])

fid = fopen('Solvation_Free_Energy.dat', 'wt');
fprintf(fid,'%s %g\n','Solvated state Energy (KJ/mol) =',solvenermol1);
fprintf(fid,'%s %g\n','Reference state Energy (KJ/mol) =',refenermol1);
fprintf(fid,'%s %g\n','Solvation Free Energy (KJ/mol) =',sfe);
fclose(fid);

movefile ('Solvation_Free_Energy.dat', outpath)

     end

% end sfe calculation

% binding energy calculation

     if (strcmp(handles.bin, 'enerbinyes')==1)

molecule=handles.mol2;
axx=strcat('mol2_',handles.filename);

    fid = fopen(axx, 'wt');
fprintf(fid,'%g %g %g\n',dime);
fprintf(fid,'%f %f %f\n', glen);
fprintf(fid,'%f %f %f %f %f %f %f %f\n', bulk);
fprintf(fid,'%s\n',boundcond);
fprintf(fid,'%g %s\n',[handles.digit handles.selectsolver]);
fprintf(fid,'%s\n',handles.dielx);
fprintf(fid,'%s\n',handles.diely);
fprintf(fid,'%s\n',handles.dielz);
fprintf(fid,'%s\n',handles.kappa);
fprintf(fid,'%s\n',handles.charge);
fprintf(fid,'%s\n',molecule);
fprintf(fid,'%s\n',moleculeref);
fprintf(fid,'%s\n',handles.srfm);
fprintf(fid,'%s\n',handles.chgm);
fprintf(fid,'%s\n',energy);
fprintf(fid,'%s\n',handles.pathdielx);
fprintf(fid,'%s\n',handles.cout);
fprintf(fid,'%f %f %f %f %f\n', handles.writeout);
fclose(fid);

        if strcmp(handles.focus, 'focusname.inm')==1

fid = fopen(handles.filename2, 'wt');
fprintf(fid,'%g %g %g\n',dime);
fprintf(fid,'%f %f %f\n', cglen);
fprintf(fid,'%f %f %f %f %f %f %f %f\n', bulk);
fprintf(fid,'%s\n',handles.bc);
fprintf(fid,'%g %s\n',[handles.digit handles.selectsolver]);
fprintf(fid,'%s\n',handles.dielx);
fprintf(fid,'%s\n',handles.diely);
fprintf(fid,'%s\n',handles.dielz);
fprintf(fid,'%s\n',handles.kappa);
fprintf(fid,'%s\n',handles.charge);
fprintf(fid,'%s\n',molecule);
fprintf(fid,'%s\n',moleculeref);
fprintf(fid,'%s\n',handles.srfm);
fprintf(fid,'%s\n',handles.chgm);
fprintf(fid,'%s\n',energy);
fprintf(fid,'%s\n',handles.pathdielx);
fprintf(fid,'%s\n',handles.cout);
fprintf(fid,'%f %f %f %f %f\n', handles.writeout);
fclose(fid);

        end


     handles.check=3;
     guidata(hObject, handles);    
 
MPBEC(axx)
solvenermol2=enerT;

molecule=handles.mol3;
axx=strcat('complex_',handles.filename);

    fid = fopen(axx, 'wt');
fprintf(fid,'%g %g %g\n',dime);
fprintf(fid,'%f %f %f\n', glen);
fprintf(fid,'%f %f %f %f %f %f %f %f\n', bulk);
fprintf(fid,'%s\n',boundcond);
fprintf(fid,'%g %s\n',[handles.digit handles.selectsolver]);
fprintf(fid,'%s\n',handles.dielx);
fprintf(fid,'%s\n',handles.diely);
fprintf(fid,'%s\n',handles.dielz);
fprintf(fid,'%s\n',handles.kappa);
fprintf(fid,'%s\n',handles.charge);
fprintf(fid,'%s\n',molecule);
fprintf(fid,'%s\n',moleculeref);
fprintf(fid,'%s\n',handles.srfm);
fprintf(fid,'%s\n',handles.chgm);
fprintf(fid,'%s\n',energy);
fprintf(fid,'%s\n',handles.pathdielx);
fprintf(fid,'%s\n',handles.cout);
fprintf(fid,'%f %f %f %f %f\n', handles.writeout);
fclose(fid);

        if strcmp(handles.focus, 'focusname.inm')==1

fid = fopen(handles.filename2, 'wt');
fprintf(fid,'%g %g %g\n',dime);
fprintf(fid,'%f %f %f\n', cglen);
fprintf(fid,'%f %f %f %f %f %f %f %f\n', bulk);
fprintf(fid,'%s\n',handles.bc);
fprintf(fid,'%g %s\n',[handles.digit handles.selectsolver]);
fprintf(fid,'%s\n',handles.dielx);
fprintf(fid,'%s\n',handles.diely);
fprintf(fid,'%s\n',handles.dielz);
fprintf(fid,'%s\n',handles.kappa);
fprintf(fid,'%s\n',handles.charge);
fprintf(fid,'%s\n',molecule);
fprintf(fid,'%s\n',moleculeref);
fprintf(fid,'%s\n',handles.srfm);
fprintf(fid,'%s\n',handles.chgm);
fprintf(fid,'%s\n',energy);
fprintf(fid,'%s\n',handles.pathdielx);
fprintf(fid,'%s\n',handles.cout);
fprintf(fid,'%f %f %f %f %f\n', handles.writeout);
fclose(fid);

        end


     handles.check=3;
     guidata(hObject, handles);    
 
MPBEC(axx)
solvenercomplex=enerT;

binener = solvenercomplex-solvenermol1-solvenermol2;
disp(' ')
disp(['Binding energy = ', num2str(binener)])

fid = fopen('Binding_Energy.dat', 'wt');
fprintf(fid,'%s %g\n','Solvated state Energy =',solvenermol1);
fprintf(fid,'%s %g\n','Reference state Energy =',solvenermol2);
fprintf(fid,'%s %g\n','SOlvation Free Energy =',solvenercomplex);
fprintf(fid,'%s %g\n','Binding Energy =',binener);
fclose(fid);

movefile ('Binding_Energy.dat', outpath)

     end
     
%       if (strcmp(handles.pka, 'enerpkayes')==1)
% cd ..
% cd pdb2pqr
% str=computer;
% if (strcmp(str,'PCWIN')==1 ||strcmp(str,'PCWIN64')==1)
%     cd pdb2pqrwin
%     path1 = pwd;
%     execut=strcat(path1,'\pdb2pqr');
% else
%     cd pdb2pqrlinux
%     path1 = pwd;
%     execut=strcat(path1,'/pdb2pqr');
% end
% cd ..
% cd ..
% cd Main
% path4=pwd;
% path6=strcat(handles.pathdielx,handles.pdbmol1);
% %execut=strcat(path1,'/pdb2pqr');
% outputpqrfilename=regexprep(handles.pdbmol1,'.pdb','.pqr');
% path4=strcat(path4,'/',outputpqrfilename);
% forcefield='PARSE';
% pHflag=handles.ph;
% argum=sprintf('%s --ff=%s --with-ph=%f %s %s',execut,forcefield, pHflag,path6, path4);
% status=system(argum);
% movefile('*.propka',outpath)
% disp('HHHHHHHHHHHHHHHHHHHHHHHHH')
% copyfile(path4,'jota.uw','f')
% delete(outputpqrfilename);
% 
% clc
%       end

% end binding energy calculation

% close(GUI_Beginner)

   end
end

% END COMPUTE SOLUTION 

% need user guide ? (open pdf file)
function pushbutton4_Callback(hObject, eventdata, handles)

status1=0;
status2=0;
status3=0;
status4=0;
status5=0;
%  status6=0;
status=0;
status7=0;
status8=0;
status9=0;
str=computer;
if (strcmp(str,'PCWIN')==1 ||strcmp(str,'PCWIN64')==1)
    winopen('userguide.pdf');

end
if (strcmp(str,'GLNX86')==1 ||strcmp(str,'GLNXA64')==1)
status = system('acroread userguide.pdf');
if status~=0
   status1 = system('evince userguide.pdf');
else
    status1=0;
end
if status1~=0
    status2 = system('xpdf userguide.pdf');
else
    status2=0;
end
   if status2~=0
    status3 = system('kpdf userguide.pdf');
   else
       status3=0;
   end
    if status3~=0
    status4 = system('okular userguide.pdf');
    else
        status4=0;
    end
    if status4~=0
    status5 = system('epdfview userguide.pdf');
    else
        status5=0;
    end
    if status5~=0
    msgbox('Sorry. We cannot find a pdf application on your computer.','Sorry','none');
    end
end
if (strcmp(str,'MACI')==1 ||strcmp(str,'MACI64')==1)
    status7 = system('open userguide.pdf');
if status7~=0
   status8 = system('open -a Safari userguide.pdf');
else
    status8=0;
end
    if status8~=0
         status9 = system('acroread userguide.pdf');
else
    status9=0;
    end
    if status9~=0
    msgbox('Sorry. We cannot find a pdf application on your computer.','Sorry','none');
    end
end
    
% map file panel (optional)
function uipanel7_SelectionChangeFcn(hObject, eventdata, handles)

if (get(hObject,'Value') == get(hObject,'Max'))

handles.inputfiles='YES';
guidata(hObject, handles);

end
if (get(hObject,'Value') == get(hObject,'Min'))

handles.inputfiles='NO';
guidata(hObject, handles);

end

 % biconjugate option linear solver
function radiobutton33_Callback(hObject, eventdata, handles)

if (get(hObject,'Value') == get(hObject,'Max'))
handles.selectsolver='biconjgrad';
guidata(hObject, handles);
end

% Gmres option linear solver
function radiobutton29_Callback(hObject, eventdata, handles)

if (get(hObject,'Value') == get(hObject,'Max'))
handles.selectsolver='gmres';
guidata(hObject, handles);
end

% Minres option linear solver
function radiobutton31_Callback(hObject, eventdata, handles)

if (get(hObject,'Value') == get(hObject,'Max'))
handles.selectsolver='minres';
guidata(hObject, handles);
end

function radiobutton33_CreateFcn(hObject, eventdata, handles)

s1 = sprintf('Biconjugate gradient solver stabilized by LU decomposition (droptol=0.000005). Max iterations = 400');
set(hObject,'TooltipString',s1)

function radiobutton29_CreateFcn(hObject, eventdata, handles)

s1 = sprintf('Generalized minimal residual solver stabilized by LU decomposition (drotol=0.000005). Restart = 10, max iterations = 400.');
set(hObject,'TooltipString',s1)

function radiobutton31_CreateFcn(hObject, eventdata, handles)

s1 = sprintf('Minimal residual solver stabilized by LU decomposition (drotol=0.000005). Max iterations = 400.');
set(hObject,'TooltipString',s1)

function radiobutton6_CreateFcn(hObject, eventdata, handles)

% read input file name
function edit42_Callback(hObject, eventdata, handles)

handles.filename=get(hObject,'String');

guidata(hObject, handles);

function edit42_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit43_Callback(hObject, eventdata, handles)
if ((strcmp(get(hObject,'String'),'')==0 && strcmp(handles.focus,'focusname.inm')==0))
set(hObject,'String','');    
msgbox('This parameter will be ignored. To assign a value to this parameter, please click the bottom "Focus" and try again. Thanks.','!!WARNING!!');
else
handles.cglen1=str2double(get(hObject,'String'));
guidata(hObject, handles); 
end

function edit43_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit44_Callback(hObject, eventdata, handles)
if ((strcmp(get(hObject,'String'),'')==0 && strcmp(handles.focus,'focusname.inm')==0))
set(hObject,'String','');    
msgbox('This parameter will be ignored. To assign a value to this parameter, please click the bottom "Focus" and try again. Thanks.','!!WARNING!!');
else
handles.cglen2=str2double(get(hObject,'String'));
guidata(hObject, handles); 
end

function edit44_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit45_Callback(hObject, eventdata, handles)
if ((strcmp(get(hObject,'String'),'')==0 && strcmp(handles.focus,'focusname.inm')==0))
set(hObject,'String','');    
msgbox('This parameter will be ignored. To assign a value to this parameter, please click the bottom "Focus" and try again. Thanks.','!!WARNING!!');
else
handles.cglen3=str2double(get(hObject,'String'));
guidata(hObject, handles); 
end

function edit45_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function radiobutton34_Callback(hObject, eventdata, handles)
if ((get(hObject,'Value') == get(hObject,'Max') && strcmp(handles.grid,'Manual')==1))
handles.focus='focusname.inm';
guidata(hObject, handles);
end
if ((get(hObject,'Value') == get(hObject,'Max') && strcmp(handles.grid ,'Automatic')==1))
msgbox('Unclick Automatic and click Manual option and try it again. Thanks.','!!WARNING!!');
end

% read molecular structure molecule 2
function edit47_Callback(hObject, eventdata, handles)
if ((strcmp(get(hObject,'String'),'')==0 && strcmp(handles.complexsyst,'yes')==0))
set(hObject,'String','');    
msgbox('To upload this file, please click the button "Optional Uploading" and try again. Thanks.','!!WARNING!!');
else
if strcmp(handles.pdb2pqr, 'PQR')==1    
filespec = {'*.pqr', 'Format File (*.pqr)'};
[handles.mol2, handles.pathdielx] = uigetfile(filespec, 'Select the Required File',handles.memorypath);
else
    filespec = {'*.pdb', 'Format File (*.pdb)'};
[handles.mol2, handles.pathdielx] = uigetfile(filespec, 'Select the Required File',handles.memorypath);
if (ischar(handles.mol2)==0)
m1 = msgbox('Please select a file. Thanks.','!!WARNING!!'); 
 waitfor(m1)
 cd ..
cd ..
cd Main
else
cd ..
cd pdb2pqr
h = msgbox('This could be a minute. Patience, MPBEC is working for you...','Converting pdb to pqr','help');
waitfor(h)
str=computer;
if (strcmp(str,'PCWIN')==1 ||strcmp(str,'PCWIN64')==1)
    cd pdb2pqrwin
    path1 = pwd;
    execut=strcat(path1,'\pdb2pqr');
else
    cd pdb2pqrlinux
    path1 = pwd;
    execut=strcat(path1,'/pdb2pqr');
end

cd ..
cd ..
cd Main
path6=strcat(handles.pathdielx,handles.mol2);
outputpqrfilename=regexprep(handles.mol2,'.pdb','.pqr');
path4=strcat(handles.pathdielx,outputpqrfilename);
warning('off','all')
file_id = fopen('pdb2pqrscript.st', 'rt');
forcefield=fgetl(file_id);
fclose(file_id);
argum=sprintf('%s --ff=%s %s %s',execut,forcefield, path6, path4);
status=system(argum);
handles.mol2=outputpqrfilename;
warning('on','all');
clc
end
end
if (ischar(handles.mol2)==0)
    handles.mol2='';
elseif strcmp(handles.mol2, '')==0

handles.memorypath=handles.pathdielx;
end
set(hObject, 'String',handles.mol2);
guidata(hObject,handles);

end

function edit47_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% read molecular structure molecule 1 and 2 (complex)
function edit48_Callback(hObject, eventdata, handles)
if ((strcmp(get(hObject,'String'),'')==0 && strcmp(handles.complexsyst,'yes')==0))
set(hObject,'String','');    
msgbox('To upload this file, please click the button "Optional Uploading" and try again. Thanks.','!!WARNING!!');
else
    if strcmp(handles.pdb2pqr, 'PQR')==1
filespec = {'*.pqr', 'Format File (*.pqr)'};
[handles.mol3, handles.pathdielx] = uigetfile(filespec, 'Select the Required File',handles.memorypath);
    else
        filespec = {'*.pdb', 'Format File (*.pdb)'};
[handles.mol3, handles.pathdielx] = uigetfile(filespec, 'Select the Required File',handles.memorypath);
if (ischar(handles.mol3)==0)
m1 = msgbox('Please select a file. Thanks.','!!WARNING!!'); 
 waitfor(m1)
 cd ..
cd ..
cd Main
else
cd ..
cd pdb2pqr
h = msgbox('This could be a minute. Patience, MPBEC is working for you...','Converting pdb to pqr','help');
waitfor(h)
str=computer;
if (strcmp(str,'PCWIN')==1 ||strcmp(str,'PCWIN64')==1)
    cd pdb2pqrwin
    path1 = pwd;
    execut=strcat(path1,'\pdb2pqr');
else
    cd pdb2pqrlinux
    path1 = pwd;
    execut=strcat(path1,'/pdb2pqr');
end

cd ..
cd ..
cd Main
path6=strcat(handles.pathdielx,handles.mol3);
outputpqrfilename=regexprep(handles.mol3,'.pdb','.pqr');
path4=strcat(handles.pathdielx,outputpqrfilename);
warning('off','all')
file_id = fopen('pdb2pqrscript.st', 'rt');
forcefield=fgetl(file_id);
fclose(file_id);
argum=sprintf('%s --ff=%s %s %s',execut,forcefield, path6, path4);
status=system(argum);
handles.mol3=outputpqrfilename;
warning('on','all');
clc
end
    end
if (ischar(handles.mol3)==0)
    handles.mol3='';
elseif strcmp(handles.mol3, '')==0

handles.memorypath=handles.pathdielx;
end
set(hObject, 'String',handles.mol3);
guidata(hObject,handles);

end

function edit48_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% upload molecular strucutre molecule 2 ? YES
function radiobutton35_Callback(hObject, eventdata, handles)
if (get(hObject,'Value') == get(hObject,'Max'))
handles.complexsyst='yes';
guidata(hObject, handles);
end

% Calculate solvation energy molecule 1 ?
function radiobutton86_Callback(hObject, eventdata, handles)
if (get(hObject,'Value') == get(hObject,'Max'))
handles.sfe='enersfeyes';
guidata(hObject, handles);
else
 handles.sfe='enersfeno';
guidata(hObject, handles);
end

% Calculate binding energy between molecule 1 and 2?
function radiobutton87_Callback(hObject, eventdata, handles)
if (get(hObject,'Value') == get(hObject,'Max'))
handles.bin='enerbinyes';
guidata(hObject, handles);
else
 handles.bin='enerbinno';
guidata(hObject, handles);
end

% Calculate pKa ?
function radiobutton88_Callback(hObject, eventdata, handles)
if (get(hObject,'Value') == get(hObject,'Max'))
handles.pka='enerpkayes';
guidata(hObject, handles);
else
 handles.pka='enerpkano';
guidata(hObject, handles);
end

% Calculate electrostatic potential molecule 1 ?
function radiobutton89_Callback(hObject, eventdata, handles)
if (get(hObject,'Value') == get(hObject,'Max'))
handles.epe='enerepeyes';
guidata(hObject, handles);
else
 handles.epe='enerepeno';
guidata(hObject, handles);
end

function radiobutton89_CreateFcn(hObject, eventdata, handles)
s1 = sprintf('It calculates the total electrostatic free energy stored in a single molecule. \n This energy basically represents the work that it takes to assemble the molecule. ');
set(hObject,'TooltipString',s1)


function radiobutton86_CreateFcn(hObject, eventdata, handles)
s1 = sprintf('The solvent effects on biomolecules are usually described by the solvation free energy. \n In this panel option, MPBEC calculates the polar contribution to the free energy released \n during the solvation of the biomolecule from vacuum to the solvent of interest. ');
set(hObject,'TooltipString',s1)


function radiobutton87_CreateFcn(hObject, eventdata, handles)
s1 = sprintf('It calculates the polar contribution to the binding free energy between two molecules (protein and ligand). \n Although this approach does not account for conformational changes of the complex after formation, \n it is a good approximation to study interactions between biomolecules.');
set(hObject,'TooltipString',s1)


function radiobutton88_CreateFcn(hObject, eventdata, handles)
s1 = sprintf('Protein active sites are able to protonate / deprotonate protons depending on the \n number of hydrogen ions (pH level) present in the aqueous electrolyte solution. This titration \n process may affect the function, stability and molecular recognition of proteins. MPBEC \n uses PROPKA to calculate the protonation (positioning of hydrogen ions) \n to a given molecular structure in PDB format for any pH level. PROPKA generates a new \n PDB file containing the renormalized protein charge distribution which is used to solve the linear PB equation. ');
set(hObject,'TooltipString',s1)

% ????
function edit49_Callback(hObject, eventdata, handles)

function edit49_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% ????
function radiobutton93_Callback(hObject, eventdata, handles)
if (get(hObject,'Value') == get(hObject,'Max'))

handles.pdb2pqr='PDB';
guidata(hObject, handles);

end
if (get(hObject,'Value') == get(hObject,'Min'))

handles.pdb2pqr='PQR';
guidata(hObject, handles);

end
 
% ?????
function radiobutton92_Callback(hObject, eventdata, handles)
 if (get(hObject,'Value') == get(hObject,'Max'))

handles.pdb2pqr='PQR';
guidata(hObject, handles);

 end
  if (get(hObject,'Value') == get(hObject,'Min'))

handles.pdb2pqr='PDB';
guidata(hObject, handles);

 end

% ??????
function radiobutton93_CreateFcn(hObject, eventdata, handles)

function radiobutton92_CreateFcn(hObject, eventdata, handles)

function edit50_Callback(hObject, eventdata, handles)

function edit50_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% read pH for pK a calculation
function edit51_Callback(hObject, eventdata, handles)
handles.ph=str2double(get(hObject,'String'));
guidata(hObject, handles);guidata(hObject, handles);

function edit51_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% visualize molecule 1 ?
function togglebutton2_Callback(hObject, eventdata, handles)
if (strcmp(handles.mol1,'')==1)
   
msgbox('please upload the corresponding pdb/pqr file and try again. Thanks.','!!WARNING!!');
else
    path6=strcat(handles.pathdielx,handles.mol1);
argum=sprintf('java -jar Jmol.jar --silent %s',path6);
status=system(argum);
clc
end

% visualize molecule 2?
function togglebutton3_Callback(hObject, eventdata, handles)
if ( strcmp(handles.mol2,'')==1)
set(hObject,'String','');    
msgbox('please upload the corresponding pdb/pqr file and try again. Thanks.','!!WARNING!!');
else
    path6=strcat(handles.pathdielx,handles.mol2);
argum=sprintf('java -jar Jmol.jar --silent %s',path6);
status=system(argum);
clc
end

% visualize molecule 1 and 2 ? (complex)
function togglebutton4_Callback(hObject, eventdata, handles)
if ( strcmp(handles.mol3,'')==1)
set(hObject,'String','');    
msgbox('please upload the corresponding pdb/pqr file and try again. Thanks.','!!WARNING!!');
else
    path6=strcat(handles.pathdielx,handles.mol3);
argum=sprintf('java -jar Jmol.jar --silent %s',path6);
status=system(argum);
clc
end

% PDB or PQR panel
function uipanel35_SelectionChangeFcn(hObject, eventdata, handles)
if (get(hObject,'Value') == get(hObject,'Min'))

handles.pdb2pqr='PDB';
guidata(hObject, handles);

end 
if (get(hObject,'Value') == get(hObject,'Max'))

    handles.pdb2pqr='PQR';
guidata(hObject, handles);

 end

function radiobutton97_CreateFcn(hObject, eventdata, handles)
if (get(hObject,'Value') == get(hObject,'Max'))

handles.pdb2pqr='PQR';

guidata(hObject, handles);
  else
handles.pdb2pqr='PDB';
guidata(hObject, handles);
end

% ??????
function radiobutton96_CreateFcn(hObject, eventdata, handles)
if (get(hObject,'Value') == get(hObject,'Min'))

handles.pdb2pqr='PDB';
guidata(hObject, handles);
else
      handles.pdb2pqr='PQR';

guidata(hObject, handles);

 end

% read PDB files? YES
function radiobutton97_Callback(hObject, eventdata, handles)
if (get(hObject,'Value') == get(hObject,'Max'))
handles.pdb2pqr='PDB';

guidata(hObject, handles);
  else

handles.pdb2pqr='PQR';
guidata(hObject, handles);
end

% ????
function radiobutton35_CreateFcn(hObject, eventdata, handles)

function radiobutton35_DeleteFcn(hObject, eventdata, handles)


% read grid points along x
function edit55_Callback(hObject, eventdata, handles)
handles.dime1=int16(str2double(get(hObject,'String')));
guidata(hObject,handles);

function edit55_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% read grid points along y
function edit56_Callback(hObject, eventdata, handles)
handles.dime2=int16(str2double(get(hObject,'String')));
guidata(hObject,handles);

function edit56_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% read grid points along z
function edit57_Callback(hObject, eventdata, handles)
handles.dime3=int16(str2double(get(hObject,'String')));
guidata(hObject,handles);

function edit57_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% read box length along x fine grid
function edit58_Callback(hObject, eventdata, handles)
handles.glen1=str2double(get(hObject,'String'));
guidata(hObject,handles);

function edit58_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% read box length along y fine grid
function edit59_Callback(hObject, eventdata, handles)
handles.glen2=str2double(get(hObject,'String'));
guidata(hObject,handles);


function edit59_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% read box length along z fine grid
function edit60_Callback(hObject, eventdata, handles)
handles.glen3=str2double(get(hObject,'String'));
guidata(hObject,handles);


function edit60_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% % read box length along x coarse grid
function edit52_Callback(hObject, eventdata, handles)
if ((strcmp(get(hObject,'String'),'')==0 && strcmp(handles.focus,'focusname.inm')==0))
set(hObject,'String','');    
msgbox('This parameter will be ignored. To assign a value to this parameter, please click the bottom "Focus" and try again. Thanks.','!!WARNING!!');
else
handles.cglen1=str2double(get(hObject,'String'));
guidata(hObject, handles); 
end

function edit52_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% read box length along y coarse grid
function edit53_Callback(hObject, eventdata, handles)
if ((strcmp(get(hObject,'String'),'')==0 && strcmp(handles.focus,'focusname.inm')==0))
set(hObject,'String','');    
msgbox('This parameter will be ignored. To assign a value to this parameter, please click the bottom "Focus" and try again. Thanks.','!!WARNING!!');
else
handles.cglen2=str2double(get(hObject,'String'));
guidata(hObject, handles); 
end

function edit53_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% read box length along z coarse grid
function edit54_Callback(hObject, eventdata, handles)
if ((strcmp(get(hObject,'String'),'')==0 && strcmp(handles.focus,'focusname.inm')==0))
set(hObject,'String','');    
msgbox('This parameter will be ignored. To assign a value to this parameter, please click the bottom "Focus" and try again. Thanks.','!!WARNING!!');
else
handles.cglen3=str2double(get(hObject,'String'));
guidata(hObject, handles); 
end

function edit54_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% focus calculation ?
function radiobutton98_Callback(hObject, eventdata, handles)
if ((get(hObject,'Value') == get(hObject,'Max') && strcmp(handles.grid,'Manual')==1))
handles.focus='focusname.inm';
guidata(hObject, handles);
end
if ((get(hObject,'Value') == get(hObject,'Max') && strcmp(handles.grid ,'Automatic')==1))
msgbox('Unclick Automatic and click Manual option and try it again. Thanks.','!!WARNING!!');
set(hObject,'Value',0.0);
guidata(hObject, handles);
end
if ((get(hObject,'Value') == get(hObject,'Min') && strcmp(handles.grid,'Manual')==1))
handles.focus='';
guidata(hObject, handles);
end

function uipanel37_SelectionChangeFcn(hObject, eventdata, handles)
if (get(hObject,'Value') == get(hObject,'Min'))

handles.grid='Automatic';

guidata(hObject, handles);
end 
if (get(hObject,'Value') == get(hObject,'Max'))
% disp('M')
    handles.grid='Manual';

guidata(hObject, handles);
 end

% Automatic grid calculation ? 
function radiobutton100_Callback(hObject, eventdata, handles)
if (get(hObject,'Value') == get(hObject,'Max'))

handles.grid='Automatic';

guidata(hObject, handles);
  else
handles.grid='Manual';
  
guidata(hObject, handles);
end

function radiobutton99_Callback(hObject, eventdata, handles)

function radiobutton99_CreateFcn(hObject, eventdata, handles)


function radiobutton100_CreateFcn(hObject, eventdata, handles)

% Need to quit the job ?
function pushbutton5_Callback(tyuhObject, eventdata, handles)
global exitmpbec
button= questdlg('Ready to quit?', ...
        'Exit Dialog','Yes','No','Yes');
switch button
            case 'Yes',
              disp('Exiting MPBEC calculation. This may take a few seconds...');
               handles.quit = 1;
              exitmpbec = handles.quit;
              guidata(tyuhObject, handles);
%              quit
            case 'No',
               guidata(tyuhObject, handles);
end


% --- Executes when selected object is changed in uipanel3.
function uipanel3_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel3 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
