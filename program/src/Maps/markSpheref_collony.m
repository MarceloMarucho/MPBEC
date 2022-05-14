function eps= markSpheref_collony(eps,idpos,xmin,ymin,zmin,infl,infl2,h,markval,dime,area,numm)

%auxiliary algorithms used previously

% markSphere algorithm

% convert to the grid reference frame
atomN=size(idpos,1);
for p=1:atomN

    if area(p)>0
        npts=numm(p);
        
for ii=1:npts
posx=idpos(p,ii,1)-xmin;
posy=idpos(p,ii,2)-ymin;
posz=idpos(p,ii,3)-zmin;

xrange=infl;
yrange=infl;
zrange=infl;

% box size

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
                if ((dx2+dy2+dz2) <= infl2) 
                 eps(i+1,j+1,k+1)=markval;
                 end
        end
end

 end
end
    end
end
 
    clear imin imax jmin jmax kmin kmax posx posy posz



