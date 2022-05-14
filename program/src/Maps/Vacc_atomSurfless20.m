% Vacc_atomSurf

if (prad > MAX_PROBERADIUS) % check this!!!
    disp('probe radius exceeds the maximum allowed')
    return
end

% allocation
area=zeros(atomN,1);
numm=zeros(atomN,1);
poson=zeros(atomN,refshenpts,3);
 fprintf('Calculating accessible surface')
for p=1:atomN
    if mod (p,5)==0 
        fprintf('.')
    end
id=p;
    rad=atomR(p)+srad;
j=1; % counter

npts=0; % counter
 
for i=1:refshenpts

    % lets define sphere points of radius rad centered at atom position
    % atomP
    
    pos(1)=rad*xpts(i)+atomP(p,1); 
    pos(2)=rad*ypts(i)+atomP(p,2);
    pos(3)=rad*zpts(i)+atomP(p,3);
    
% let check if this point is within the union of the spheres centered at
% the atomic centers with radiii rad 

% run vaccpoints % it returns acc=1 if accessible (outside the inflated van 
% der walls radius), otherwise acc=0.


    acc=1;
for ih=1:atomN
    distcheck2 = (pos(1)-atomP(ih,1))^2+(pos(2)-atomP(ih,2))^2+...
        (pos(3)-atomP(ih,3))^2;
     if (distcheck2 < (atomR(ih)+srad)^2) 
        if (ih~=id) % (the two atoms should be differents)
            acc=0;
            break
        end
    end
end



if (acc==1)
    npts=npts+1;
    
    % assign the points
    
poson(p,j,1)=pos(1);   
poson(p,j,2)=pos(2);
poson(p,j,3)=pos(3);
j=j+1;
end

end
area(p) = 4.*pi*rad^2*double(npts)/double(refshenpts);
numm(p)=npts;


end
fprintf('\nDone!\n')
clear pos
