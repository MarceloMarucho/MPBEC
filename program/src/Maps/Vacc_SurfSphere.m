% Vacc_SurfSphere

% number of points over the surface nipts (input). The actual number of
% points is given by nactual. finally the real number of points
% npts=nactual

% set up constant
frac=double(refshenpts)/4.;
ntheta=int16(sqrt(pi*frac)); 
dtheta=pi/double(ntheta);
nphimax=2*ntheta;

% count the actual number of points

nactual=0;
for itheta=0:ntheta-1
    theta=dtheta*double(itheta);
    sintheta=sin(theta);
%    costheta=cos(theta);
    nphi=int16(sintheta*nphimax);
    nactual=nactual+nphi;
end
 refshenpts=nactual;

%preallocating array

xpts=zeros(nactual,1); 
ypts=zeros(nactual,1); 
zpts=zeros(nactual,1);
    
nactual=1;
for itheta=0:ntheta-1
    theta=dtheta*double(itheta);
    sintheta=sin(theta);
    costheta=cos(theta);
    nphi=int16(sintheta*nphimax);
    if (nphi~=0)
        dphi= 2.0*pi/double(nphi);
        for iphi=0:nphi-1
            phi=dphi*double(iphi);
            sinphi=sin(phi);
            cosphi=cos(phi);
            xpts(nactual)=cosphi*sintheta;
            ypts(nactual)=sinphi*sintheta;
            zpts(nactual)=costheta;
            nactual=nactual+1;
        end
    end
end


