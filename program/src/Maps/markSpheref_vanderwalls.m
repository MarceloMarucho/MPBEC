function eps= markSpheref_vanderwalls(eps,idpos,xmin,ymin,zmin,xmax,ymax,zmax,infl,infl2,h,markval,dime,atomR,VSMALL)


% convert to the grid reference frame
npts=size(idpos,1);

if markval~=0
    % dielectric map
for ii=1:npts
    if ((idpos(ii,1) >=xmin) && (idpos(ii,1)<= xmax)&& (idpos(ii,2) >=ymin) && (idpos(ii,2)...
            <= ymax) && (idpos(ii,3) >=zmin) && (idpos(ii,3) <= zmax))
        if (atomR(ii) > VSMALL)
posx=idpos(ii,1)-xmin;
posy=idpos(ii,2)-ymin;
posz=idpos(ii,3)-zmin;

xrange=infl(ii);
yrange=infl(ii);
zrange=infl(ii);


imin=max(0,ceil((posx-xrange)/h(1)));
jmin=max(0,ceil((posy-yrange)/h(2)));
kmin=max(0,ceil((posz-zrange)/h(3)));

imax=min(dime(1)-1,floor((posx+xrange)/h(1)));
jmax=min(dime(2)-1,floor((posy+yrange)/h(2)));
kmax=min(dime(3)-1,floor((posz+zrange)/h(3)));


 for k=kmin:kmax
             dz=posz-h(3)*k;
             dz2=dz^2;
for i=imin:imax
        dx=posx-h(1)*i;
        dx2=dx^2;
        for j=jmin:jmax
            dy=posy-h(2)*j;
            dy2=dy^2;
                if ((dx2+dy2+dz2) <= infl2(ii)) 
                 eps(i+1,j+1,k+1)=markval;
                 end
        end
end

 end
        end
   end
end

else
    % kappa map
    
for ii=1:npts
    if ((idpos(ii,1) >=(xmin-infl(ii))) && (idpos(ii,1)<=( xmax+infl(ii)))&& (idpos(ii,2) >=(ymin-infl(ii))) && (idpos(ii,2)...
            <= (ymax+infl(ii))) && (idpos(ii,3) >=(zmin-infl(ii))) && (idpos(ii,3) <= (zmax+infl(ii))))
        if (atomR(ii) > VSMALL)
posx=idpos(ii,1)-xmin;
posy=idpos(ii,2)-ymin;
posz=idpos(ii,3)-zmin;

xrange=infl(ii);
yrange=infl(ii);
zrange=infl(ii);


imin=max(0,ceil((posx-xrange)/h(1)));
jmin=max(0,ceil((posy-yrange)/h(2)));
kmin=max(0,ceil((posz-zrange)/h(3)));

imax=min(dime(1)-1,floor((posx+xrange)/h(1)));
jmax=min(dime(2)-1,floor((posy+yrange)/h(2)));
kmax=min(dime(3)-1,floor((posz+zrange)/h(3)));


 for k=kmin:kmax
             dz=posz-h(3)*k;
             dz2=dz^2;
for i=imin:imax
        dx=posx-h(1)*i;
        dx2=dx^2;
        for j=jmin:jmax
            dy=posy-h(2)*j;
            dy2=dy^2;
                if ((dx2+dy2+dz2) <= infl2(ii)) 
                 eps(i+1,j+1,k+1)=markval;
                 end
        end
end

 end
        end
   end
end

end
 
    clear imin imax jmin jmax kmin kmax posx posy posz



