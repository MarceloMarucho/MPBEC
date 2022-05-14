
%charge density discretization 

% allocating arrays and vectors
charge=zeros(dime(1),dime(2),dime(3));
subcharge=zeros(prod(dime),1);
%
%memory_in_use = monitor_memory_whos
%if (memory_in_use>500)
%    disp('running out of memory')
%end

dime12=dime(1)*dime(2);

if (strcmp(chgm,'SPL0')==1)
    fprintf('Processing charge data')
for p=1:atomN
    if mod (p,50)==0 
        fprintf('.')
    end
    if ((atomP(p,1) >=xmin) && (atomP(p,1)<= xmax)&& (atomP(p,2) >=ymin) && (atomP(p,2)...
            <= ymax) && (atomP(p,3) >=zmin) && (atomP(p,3) <= zmax))

ifloat=Pref(p,1)/h(1)+1; % note there is a plus 1 that is not in the APBS source
jfloat=Pref(p,2)/h(2)+1;
kfloat=Pref(p,3)/h(3)+1;

ihi=ceil(ifloat);
ilo=floor(ifloat);
jhi=ceil(jfloat);
jlo=floor(jfloat);
khi=ceil(kfloat);
klo=floor(kfloat);

dx=ifloat-double(ilo);
dy=jfloat-double(jlo);
dz=kfloat-double(klo);


 partcharge=atomC(p)/prod(h);
     
      eme=(khi-1)*dime12+(jhi-1)*(dime(1))+ihi;
     subcharge(eme)=subcharge(eme)+dx*dy*dz*partcharge;
      eme=(khi-1)*dime12+(jlo-1)*(dime(1))+ihi;
     subcharge(eme)=subcharge(eme)+dx*(1.-dy)*dz*partcharge;
      eme=(klo-1)*dime12+(jlo-1)*(dime(1))+ihi;
     subcharge(eme)=subcharge(eme)+dx*(1.-dy)*(1.-dz)*partcharge;
      eme=(khi-1)*dime12+(jhi-1)*(dime(1))+ilo;
     subcharge(eme)=subcharge(eme)+(1.-dx)*dy*dz*partcharge;
      eme=(khi-1)*dime12+(jlo-1)*(dime(1))+ilo;
     subcharge(eme)=subcharge(eme)+(1.-dx)*(1.-dy)*dz*partcharge;
      eme=(klo-1)*dime12+(jhi-1)*(dime(1))+ihi;
     subcharge(eme)=subcharge(eme)+dx*dy*(1.-dz)*partcharge;
      eme=(klo-1)*dime12+(jhi-1)*(dime(1))+ilo;
     subcharge(eme)=subcharge(eme)+(1.-dx)*dy*(1.-dz)*partcharge;
      eme=(klo-1)*dime12+(jlo-1)*(dime(1))+ilo;
     subcharge(eme)=subcharge(eme)+(1.-dx)*(1.-dy)*(1.-dz)*partcharge;

    end
end
clear ihi ilo jhi jlo khi klo dx dy dz ifloat jfloat kfloat
end

if (strcmp(chgm,'SPL2')==1)  
    
    vfchi=@(x,y) 1.5+(double(x)-y);
%
 fprintf('Processing charge data')
for p=1:atomN
    if mod (p,150)==0 
        fprintf('.')
    end
    if ((atomP(p,1) >= (xmin-h(1))) && (atomP(p,1)<= (xmax+h(1)))&& ...
            (atomP(p,2) >=(ymin-h(2))) && (atomP(p,2)...
    <= (ymax+h(2))) && (atomP(p,3) >=(zmin-h(3))) && (atomP(p,3) <= (zmax+h(3))))

partcharge=atomC(p)/prod(h);

% let figure out which vertices we are next to

ifloat=Pref(p,1)/h(1)+1; % note there is a plus 1 that is not in the APBS source
jfloat=Pref(p,2)/h(2)+1;
kfloat=Pref(p,3)/h(3)+1;

ip1=ceil(ifloat);
ip2=ip1+1;
im1=floor(ifloat);
im2=im1-1;
jp1=ceil(jfloat);
jp2=jp1+1;
jm1=floor(jfloat);
jm2=jm1-1;
kp1=ceil(kfloat);
kp2=kp1+1;
km1=floor(kfloat);
km2=km1-1;


for ii=im2:ip2
    arg=vfchi(ii,ifloat);
    
    % lets evaluate the bispline
    m2m=0.0;
    m2=0.0;
    mx=0.0;
    if ((arg>=0.0)&&(arg<=2.0)) 
        m2m=1.0-abs(arg-1.0);
    else
        m2m=0.0;
    end
    if ((arg>=1.0)&&(arg<=3.0)) 
        m2=1.0-abs(arg-2.0);
    else
        m2=0.0;
    end
    if ((arg>=0.0)&&(arg<=3.0))
        mx=0.5*arg*m2m+0.5*(3.0-arg)*m2;
    else
        mx=0.0;
    end
    
    % done
    
    for jj=jm2:jp2
            arg=vfchi(jj,jfloat);
    
    % lets evaluate the bispline
    m2m=0.0;
    m2=0.0;
    my=0.0;
    if ((arg>=0.0)&&(arg<=2.0)) 
        m2m=1.0-abs(arg-1.0);
    else
        m2m=0.0;
    end
    if ((arg>=1.0)&&(arg<=3.0)) 
        m2=1.0-abs(arg-2.0);
    else
        m2=0.0;
    end
    if ((arg>=0.0)&&(arg<=3.0))
        my=0.5*arg*m2m+0.5*(3.0-arg)*m2;
    else
        my=0.0;
    end
    
    % done

        for kk=km2:kp2
         arg=vfchi(kk,kfloat);
    
    % lets evaluate the bispline
    m2m=0.0;
    m2=0.0;
    mz=0.0;
    if ((arg>=0.0)&&(arg<=2.0)) 
        m2m=1.0-abs(arg-1.0);
    else
        m2m=0.0;
    end
    if ((arg>=1.0)&&(arg<=3.0)) 
        m2=1.0-abs(arg-2.0);
    else
        m2=0.0;
    end
    if ((arg>=0.0)&&(arg<=3.0))
        mz=0.5*arg*m2m+0.5*(3.0-arg)*m2;
    else
        mz=0.0;
    end
    
    % done
    
    % lets evaluate the charge coefficients
    
            pepe=(kk-1)*dime12+(jj-1)*dime(1)+ii;
            subcharge(pepe)=subcharge(pepe)+partcharge*mx*my*mz;
        end
    end
end
    end
end
clear vfchi m2 m2m mx my mz arg
clear ifloat jfloat kfloat
end
fprintf('\nDone!\n')
% convert vector to array format
sdime12=dime(1)*dime(2);
for i=1:dime(1)
    for j=1:dime(2)
        for k=1:dime(3)
            pepe=(k-1)*sdime12+(j-1)*dime(1)+i;
            charge(i,j,k)=subcharge(pepe);
        end
    end
end
% run memory_check
clear subcharge 

