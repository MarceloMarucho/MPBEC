% fillcoCoefSpline (smooth function based on Roux article)
if (strcmp(srfm,'SPL2')==1)

% % preallocating arrays
% 
kappa=ones(dime(1),dime(2),dime(3));
dielx=kappa;
diely=kappa;
dielz=kappa;

irad=maxionR; 

w2i=1.0/splineWin^2;
w3i=w2i/splineWin;
dime12=dime(1)*dime(2);

if (bulkIonicStrength > VPMGSMALL)
    ionmask=1.0;
else
    ionmask=0.0;
end

% loop through the atoms an do assign the dielectric
for p=1:atomN
% if we are on the grid
arad=atomR(p);
 if ((atomP(p,1) >=xmin) && (atomP(p,1)<= xmax)&& (atomP(p,2) >=ymin) && (atomP(p,2)...
            <= ymax) && (atomP(p,3) >=zmin) && (atomP(p,3) <= zmax))
        if ( arad > VPMGSMALL ) 
            
            % mark ion accessibility and dielectric values for later
            % assigment
            
            infl=irad+arad;
            itot=infl+splineWin; 
            itot2=itot^2;
            ictot=max(0,(infl-splineWin));
            ictot2=ictot^2;
            stot=arad+splineWin;
            stot2=stot^2;
            sctot=max(0,(arad-splineWin));
            sctot2=sctot^2;
            
            % we will search over the grid points which are in the greater
            % of these two radii
            
            rtot=max(itot,stot);
            rtot2=max(itot2,stot2);
            dx=rtot+0.5*h(1);
            dy=rtot+0.5*h(2);
            dz=rtot+0.5*h(3);
            
imin=max(0,floor((Pref(p,1)-dx)/h(1))); 
jmin=max(0,floor((Pref(p,2)-dy)/h(2)));
kmin=max(0,floor((Pref(p,3)-dz)/h(3)));

imax=min(dime(1)-1,ceil((Pref(p,1)+dx)/h(1))); % why here min goes with floor 
% and marksphere is with ceil
jmax=min(dime(2)-1,ceil((Pref(p,2)+dy)/h(2))); % dont I need to add a plus 1 
% like in discretitation script?
kmax=min(dime(3)-1,ceil((Pref(p,3)+dz)/h(3)));

for i=imin:imax
    ii=i+1;
        dx=Pref(p,1)-h(1)*i;
        dx2=dx^2;
        for j=jmin:jmax
            jj=j+1;
            dy=Pref(p,2)-h(2)*j;
            dy2=dy^2;
            for k=kmin:kmax
                kk=k+1;
                        dz=Pref(p,3)-h(3)*k;
                        dz2=dz^2;    
                        
                        % assign CCF

                        if (kappa(ii,jj,kk) > VPMGSMALL )
                            dist2=dx2+dy2+dz2;
                            if (dist2>= itot2)
                                kappa(ii,jj,kk)=1.0;
                            end
                            if (dist2<=ictot2)
                                kappa(ii,jj,kk)=0.0;
                            end
                            if ((dist2<itot2) && (dist2>ictot2))
                            dist=sqrt(dist2);
                            sm=dist-infl+splineWin;
                            sm2=sm^2;
                            kappa(ii,jj,kk)=0.75*sm2*w2i-0.25*sm*sm2*w3i;
                            end
                        end
                        
                        %assign A1CF
                        
                         if (dielx(ii,jj,kk) > VPMGSMALL )
                             dxx=Pref(p,1)-h(1)*(i+0.5);
                             dxx2=dxx^2;
                            dist2=dy2+dz2+dxx2;
                            if (dist2>= stot2)
                                dielx(ii,jj,kk)=1.0;
                            end
                            if (dist2<=sctot2)
                                dielx(ii,jj,kk)=0.0;
                            end
                            if ((dist2<stot2) && (dist2>sctot2))
                            dist=sqrt(dist2);
                            sm=dist-arad+splineWin;
                            sm2=sm^2;
                            dielx(ii,jj,kk)=0.75*sm2*w2i-0.25*sm*sm2*w3i;
                            end
                         end
                        
                          %assign A2CF
                        
                         if (diely(ii,jj,kk) > VPMGSMALL )
                             dyy=Pref(p,2)-h(2)*(j+0.5);
                             dyy2=dyy^2;
                            dist2=dyy2+dz2+dx2;
                            if (dist2>= stot2)
                                diely(ii,jj,kk)=1.0;
                            end
                            if (dist2<=sctot2)
                                diely(ii,jj,kk)=0.0;
                            end
                            if ((dist2<stot2) && (dist2>sctot2))
                            dist=sqrt(dist2);
                            sm=dist-arad+splineWin;
                            sm2=sm^2;
                            diely(ii,jj,kk)=0.75*sm2*w2i-0.25*sm*sm2*w3i;
                            end
                         end
                        
                          %assign A3CF
                        
                         if (dielz(ii,jj,kk) > VPMGSMALL )
                             dzz=Pref(p,3)-h(3)*(k+0.5);
                             dzz2=dzz^2;
                            dist2=dy2+dzz2+dx2;
                            if (dist2>= stot2)
                                dielz(ii,jj,kk)=1.0;
                            end
                            if (dist2<=sctot2)
                                dielz(ii,jj,kk)=0.0;
                            end
                            if ((dist2<stot2) && (dist2>sctot2))
                            dist=sqrt(dist2);
                            sm=dist-arad+splineWin;
                            sm2=sm^2;
                            dielz(ii,jj,kk)=0.75*sm2*w2i-0.25*sm*sm2*w3i;
                            end
                        end

            end
        end
end
        end
 end
end


% fill coefficients

dielx=dielx*(dielw-dielp)+dielp;
kappa=kappa*ionmask;
diely=diely*(dielw-dielp)+dielp;
dielz=dielz*(dielw-dielp)+dielp;

end