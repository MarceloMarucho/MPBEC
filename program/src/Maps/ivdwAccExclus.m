function acc = ivdwAccExclus(cxp,cyp,czp,pos,atomP,srad,head,ll,Lc,atomR,id)

% ivdwAccExclus

% determine if a point (to test "pos") is within the union of thespheres
% centered at the atomic centers with radii equal to the sum of their van
% der waals radii and probe radius "radius". atomID is the ID of the atom
% to be ignored. 

% get the relevant cell that contain the position to test from the cell
% list

%we can only test probes with radii less than the max probe radius specified

%if (prad <= 5.) % check this!!!
% if (head(cxp,cyp,czp)==-1) % there is no atoms nearby and this area point "center" is accessible 
%     acc=1;
%     return
% else
    % otherwise we have to determine if it is accessible or not.
    % let Check in the cell and the neiborgh
for cxpne=cxp-1:cxp+1
    if (cxpne<1 || cxpne>Lc(1))
        continue
    end
    for cypne=cyp-1:cyp+1
        if (cypne<1 || cypne>Lc(2))
        continue
        end
        for czpne=czp-1:czp+1
            if (czpne<1 || czpne>Lc(3))
            continue
            end
    % TODO for periodic boundary condition I have to pull back the outer
    % cells!!!
ih=head(cxpne,cypne,czpne);
while (ih~=-1)
    distcheck2=(pos(1)-atomP(ih,1))^2+(pos(2)-atomP(ih,2))^2+...
        (pos(3)-atomP(ih,3))^2;
    if (distcheck2 < (atomR(ih)+srad)^2) % how to get this information?
        if (ih~=id) % (the two atoms should be differents)
            acc=0;
            return
        end
    end
    ih=ll(ih);
end
        end
    end
end
    acc=1; % if we are still here then the point is accessible
    return
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear head ll


    
    