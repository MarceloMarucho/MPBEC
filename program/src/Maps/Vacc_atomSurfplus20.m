% Vacc_atomSurf

% imput srad (solvent radius) and the data of the specific atom of the molecule
% output the solvent accessibility surface area
% Vacc_atomSASA calls Vacc_SASA and this one calls Vacc_atomSurf which
% calls ivdwAccExlcus. The later check which points  will contribute to the
% solvent accessible surface. 
% the output is the number of points "npts", the corresponding points 
% xpts,ypts and zpts and and the value of the atomic surface area

% if (arad < VSMALL) return Vacc_Surf_ctor

% prad=maxradprobe; % is this correct? yes
% mrad=maxradatom+maxradprobe; % radii equal to the sum of the atomic van der waals radius and the probe radius
% % lets generate "refsphenpts" points "xpts,ypts, and zpts" for e reference
% % sphere of unit radius. Where refspherenpts is set up? I think it is the
% % same for all the atoms in the molecule. If so, it has to be run out of the atomic loop
% % Otherwise the number of points would dependon the specific atom.
% %rad=rad-0.004;
% %rad=3.52;
% refshenpts=ceil(4.0*pi*mrad^2*double(surf_density))
% run Vacc_SurfSphere % output xpts,ypts,zpts and refsphenpts
% % Lets determine which points will contribute

if (prad > MAX_PROBERADIUS) % check this!!!
    disp('probe radius exceeds the maximum allowed')
    return
end
 
area=zeros(atomN,1);
numm=zeros(atomN,1);
poson=zeros(atomN,refshenpts,3);

%memory_in_use = monitor_memory_whos
%if (memory_in_use>500)
%    disp('running out of memory')
%end

PII=4.*pi/double(refshenpts);
dist=(atomR+srad).^2;

fprintf('Calculating accessible surface')
for p=1:atomN
    if mod (p,60)==0 
        fprintf('.')
    end

%id=p;
    rad=atomR(p)+srad;
%j=1; % counter

npts=0; % counter

for i=1:refshenpts

    % lets define sphere points of radius rad centered at atom position apos
    
    pos(1)=rad*xpts(i)+atomP(p,1); 
    pos(2)=rad*ypts(i)+atomP(p,2);
    pos(3)=rad*zpts(i)+atomP(p,3);
    
% let check if this point is within the union of the spheres centered at
% the atomic centers with radiii rad 

% run vaccpoints % it returns acc=1 if accessible (outside the inflated van 
% der walls radius), otherwise acc=0.


cxp=floor((pos(1)-lower_cornerh(1))/rn(1))+1;
cyp=floor((pos(2)-lower_cornerh(2))/rn(2))+1;
czp=floor((pos(3)-lower_cornerh(3))/rn(3))+1;

acc=1;
%we can only test probes with radii less than the max probe radius specified

%if (head(cxp,cyp,czp)~=-1) % there is no atoms nearby and this area point "center" is accessible 
%      acc=1;
%      break
% % %    return
%  else
    % otherwise we have to determine if it is accessible or not.
    % let Check in the cell and the neiborgh
    if (cxp-1==0)
        lix=0;
    else
        lix=-1;
    end
    
        if (cyp-1==0)
        liy=0;
    else
        liy=-1;
        end
        
        if (czp-1==0)
        liz=0;
    else
        liz=-1;
        end
    
 for cxpne=lix:1
            if (acc==0)
           continue
            end
    for cypne=liy:1
            if (acc==0)
           continue
            end
       for czpne=liz:1
            if (acc==0)
           continue
            end
    % TODO for periodic boundary condition I have to pull back the outer
    % cells !!!
ih=head(cxp+cxpne,cyp+cypne,czp+czpne);
while (ih~=-1)
    distcheck2=(pos(1)-atomP(ih,1))^2+(pos(2)-atomP(ih,2))^2+...
        (pos(3)-atomP(ih,3))^2;
    if (distcheck2 < dist(ih)) 
        if (ih~=p) % (the two atoms should be differents)
            acc=0;
            break
        end
    end
    ih=ll(ih);
end
        end
    end
end

%end

if (acc==1)
    npts=npts+1;
   
    % assign the points
    
poson(p,npts,1)=pos(1);   
poson(p,npts,2)=pos(2);
poson(p,npts,3)=pos(3);

end

% bpts(i)=acc;
%end
end

area(p) = PII*rad^2*double(npts);
numm(p)=npts;

end
fprintf('\nDone!\n')
%run memory_check
clear pos dist
clear head ll