
irad=maxionR; % maximum ion radii

if (bulkIonicStrength > VPMGSMALL)
    ionmask=1.0;
else
    ionmask=0.0;
end

% preallocating arrays

kappa=ones(dime(1),dime(2),dime(3))*ionmask;

if (bulkIonicStrength > VPMGSMALL)
    
% loop through the atoms and set kappa=0 (inaccessible) if test point is
% inside the inflated van der waals radii
markval=0.0; % required by markSphere script
            infl=atomR+irad; %used by markSPhere
            infl2=infl.^2;
              idpos=atomP; %required by markSphere
kappa=markSpheref_vanderwalls(kappa,idpos,xmin,ymin,zmin,xmax,ymax,zmax,infl,infl2,h,markval,dime,atomR,VPMGSMALL);
end

clear markval arad irad idpos infl infl2