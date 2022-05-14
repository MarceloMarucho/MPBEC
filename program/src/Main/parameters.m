% parameters

% define parameters
e_c = 4.803242384e-10;  %statcoulombs
k_B = 1.380662e-16;  %erg K-1
N_A = 6.022045e+23;
dielw=bulk(3); % solvent dielctric coeficient
diel=8.8541878e-12; %F m-1 or C^2 N-1 m-2
bulkIonicStrength=bulk(1);
dielp=bulk(4); % solute dielectric coeficient
maxionR=bulk(2); % maximum ionic radius in A
T=bulk(5); % temperature in Kelvin
srad=bulk(6); % solvent molecular radius
splineWin=bulk(8); % swin

vunit_ec=1.6021773e-19;
vunit_epso=8.8541878e-12;
vunit_kb=1.3806581e-23; 

%pre1=e_c^2/(dielw*k_B*T)*1e+8;
pre1=vunit_ec^2/(4.*pi*vunit_epso*dielw*vunit_kb*T)*1.0e+10;
xkappa = sqrt((bulkIonicStrength*1.0e-16*8*pi*N_A*e_c^2)/(1000*dielw*k_B*T));
zmagic=4*pi*e_c^2/(k_B*T)*10^8; % should be derived from unit conversions similar to zmagic in VPBE.c of APBS's source code
surf_density=bulk(7); % minimum number of points to be used on the generation of the reference sphere surface
%Squared_Debye_Huckel=25897036486.0*bulk(1)*bulk(3)/T; 
Squared_Debye_Huckel=xkappa^2;
VSMALL=0.0000001;

VPMGSMALL=0.00000001;

VLARGE=10000000.0;

VAPBS_DIM=3;
VCLIST_INFLATE=1.42;
empty=-1;
maxradprobe=0.5+maxionR;
sq2=sqrt(2.);

MAXIHTABLEDIM=75;

MAX_PROBERADIUS=5.0;

% atom list

% read the pqr data

file_id = fopen(pqr_str, 'rt');
flag=[];
headerlinesIn=-1;
while isempty(flag)
    line=fgetl(file_id);
    flag=strfind(line,'ATOM');
    headerlinesIn=headerlinesIn+1;
end
fclose(file_id);

delimiterIn = ' ';

PQRdata=importdata(pqr_str,delimiterIn, headerlinesIn); 

%PQRdata=importdata(pqr_str);

atomdata=PQRdata.data;
atomN=length(atomdata(:,1));
atomP=atomdata(:,2:4);
atomC=atomdata(:,5);
atomR=atomdata(:,6);

%center of the grid
% if the center of the grid is not based on the same molecule (pqr_str) for which I am calculating the
% potential, I have to read the position of the atoms of the corresponding
% molecule (pqr_cent_str)
if strcmp(pqr_cent_str, pqr_str)==0 
% read the pqr data

file_id = fopen(pqr_cent_str, 'rt');
flag=[];
headerlinesIn=-1;
while isempty(flag)
    line=fgetl(file_id);
    flag=strfind(line,'ATOM');
    headerlinesIn=headerlinesIn+1;
end
fclose(file_id);

delimiterIn = ' ';

PQRdata2=importdata(pqr_cent_str,delimiterIn, headerlinesIn); 

%PQRdata2=importdata(pqr_cent_str);
atomdatacent=PQRdata2.data;
atomPcent=atomdatacent(:,2:4);

else
atomPcent=atomP;
end

% let's evaluate the center of the grid with respect to the molecule
% (pqr_cent_str)

xcent=(min(atomPcent(:,1))+max(atomPcent(:,1)))/2.;
ycent=(min(atomPcent(:,2))+max(atomPcent(:,2)))/2.;
zcent=(min(atomPcent(:,3))+max(atomPcent(:,3)))/2.;

% find the spatial step sizes on each direction

h=zeros(3,1);

for dimension=1:3
  h(dimension)=glen(dimension)/(dime(dimension)-1);  
end

% let's define the min/max domains

xmin=xcent-glen(1)/2.;
ymin=ycent-glen(2)/2.;
zmin=zcent-glen(3)/2.;
xmax=xcent+glen(1)/2.;
ymax=ycent+glen(2)/2.;
zmax=zcent+glen(3)/2.;

% let's convert the atom position to the grid reference system

Pref(:,1)=atomP(:,1)-xmin;
Pref(:,2)=atomP(:,2)-ymin;
Pref(:,3)=atomP(:,3)-zmin;

%vector used by the export.m file

rmin=[xmin ymin zmin];

maxradatom=max(atomR);
