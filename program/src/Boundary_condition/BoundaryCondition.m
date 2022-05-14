% what boundary condition?

   if strcmp(bc, 'sdh')==1
   bx=1;
   by=1;
   bz=1;
   end
   if strcmp(bc, 'periodic')==1
   bx=2;
   by=2;
   bz=2;
   end
   if strcmp(bc, 'mixed')==1
   bx=2;
   by=2;
   bz=3;
   end
   if strcmp(bc, 'mdh')==1
   bx=3;
   by=3;
   bz=3;
   end
    % new boundary. Periodic along x and y and Neumann along z. (means
    % bz=4)
   if strcmp(bc, 'mixed2')==1
   bx=2;
   by=2;
   bz=4;
   end

% This routine evaluates the (Dirichlet) boundary condition


%generating the grid points

X=linspace(-glen(1)/2,glen(1)/2,dime(1));
Y=linspace(-glen(2)/2,glen(2)/2,dime(2));
Z=linspace(-glen(3)/2,glen(3)/2,dime(3));

%memory allocation for arrays

gxcf=zeros(dime(2),dime(3),2);
gycf=zeros(dime(1),dime(3),2);
gzcf=zeros(dime(1),dime(2),2);
potB=zeros(dime(1),dime(2),dime(3));
BC=[1 dime(1);1 dime(2); 1 dime(3)];

%memory_in_use = monitor_memory_whos
%if (memory_in_use>500)
%    disp('running out of memory')
%end

% mdh boundary condition
% xkappa=1.5*xkappa;
%% X Face Boundary
if (bx==3)
  fprintf('Processing boundary condition along x')
for j=1:dime(2)
    if mod (j,10)==0 
        fprintf('.')
    end
    for k=1:dime(3)
        for n=1:2
            i=BC(1,n);
            
            dist=sqrt((X(i)-atomP(:,1)).^2+(Y(j)-atomP(:,2)).^2+(Z(k)-atomP(:,3)).^2);
            val=pre1*atomC./dist.*exp(-xkappa*(dist-atomR))./(1+xkappa*atomR);               
            gxcf(j,k,n)=sum(val);
            potB(i,j,k)=sum(val);
        end
    end
end
fprintf('\nDone!\n')
end

%% Y Boundary
if (by==3)
    fprintf('Processing boundary condition along y')
for i=1:dime(1)
     if mod (i,10)==0 
        fprintf('.')
    end
    for k=1:dime(3)
        for n=1:2
            j=BC(2,n);
            dist=sqrt((X(i)-atomP(:,1)).^2+(Y(j)-atomP(:,2)).^2+(Z(k)-atomP(:,3)).^2);
            val=pre1*atomC./dist.*exp(-xkappa*(dist-atomR))./(1+xkappa*atomR);                
            gycf(i,k,n)=sum(val);
            potB(i,j,k)=sum(val);
        end
    end
end
fprintf('\nDone!\n')
end

%% Z Boundary
if (bz==3)
    fprintf('Processing boundary condition along z')
for i=1:dime(1)
     if mod (i,10)==0 
        fprintf('.')
    end
    for j=1:dime(2)
        for n=1:2
            k=BC(3,n);
            dist=sqrt((X(i)-atomP(:,1)).^2+(Y(j)-atomP(:,2)).^2+(Z(k)-atomP(:,3)).^2);
            val=pre1*atomC./dist.*exp(-xkappa*(dist-atomR))./(1+xkappa*atomR);                
            gzcf(i,j,n)=sum(val);
            potB(i,j,k)=sum(val);
        end
    end
end
fprintf('\nDone!\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%

% sdh boundary condition
solutecharge=sum(atomC);
solutex=xcent;
solutey=ycent;
solutez=zcent;


if (strcmp(bc, 'sdh')==1 || (atomN>19))

    % lets determine the charge, position and radius of the spherical solute
    
%solutecharge=sum(atomC);
%solutex=xcent;
%solutey=ycent;
%solutez=zcent;

soluterad = 0;
    xminsol = atomP(1,1);
    xmaxsol = atomP(1,1);
    yminsol = atomP(1,2);
    ymaxsol = atomP(1,2);
    zminsol = atomP(1,3);
    zmaxsol = atomP(1,3);
    
    for iatom=0: atomN-1
        atomRadius=atomR(iatom+1);
        xx = atomP(iatom+1,1);
        yy = atomP(iatom+1,2);
        zz = atomP(iatom+1,3);
        if ((xx+atomRadius) > xmaxsol) 
            xmaxsol = xx + atomRadius;
        end
        if ((xx-atomRadius) < xminsol) 
            xminsol = xx - atomRadius;
        end
        if ((yy+atomRadius) > ymaxsol) 
            ymaxsol = yy + atomRadius;
        end
        if ((yy-atomRadius) < yminsol) 
            yminsol = yy - atomRadius;
        end
        if ((zz+atomRadius) > zmaxsol) 
            zmaxsol = zz + atomRadius;
        end
        if ((zz-atomRadius) < zminsol) 
            zminsol = zz - atomRadius;
        end
        dispx = (xx - solutex);
        dispy = (yy - solutey);
        dispz = (zz - solutez);
        distsol = (dispx*dispx) + (dispy*dispy) + (dispz*dispz); 
        disti = sqrt(distsol) + atomRadius;
        if (disti > soluterad) 
            soluterad = disti;
        end
    end
    dix=(xmaxsol-xminsol);
    diy=(ymaxsol-yminsol);
    diz=(zmaxsol-zminsol);
    solutelen=[dix,diy,diz];
    clear disti distsol dispx dispy dispz xx yy zz atomRadius
    clear xminsol xmaxsol yminsol ymaxsol zminsol zmaxsol
end

%% X Face Boundary

sdhcharge = 0.0;
sdhdipole = zeros(3);
sdhquadrupole=zeros(9);
traced=zeros(9);

for iatom=0: atomN-1
                xr = atomP(iatom+1,1)-solutex;
                yr = atomP(iatom+1,2)-solutey;
                zr = atomP(iatom+1,3)-solutez;
                charge = atomC(iatom+1);
                sdhcharge = sdhcharge+charge;
                sdhdipole(1)= sdhdipole(1)+xr * charge;
                sdhdipole(2) = sdhdipole(2)+yr * charge;
                sdhdipole(3)= sdhdipole(3)+zr * charge;
                traced(1) = xr*xr*charge;
                traced(2) = xr*yr*charge;
                traced(3) = xr*zr*charge;
                traced(4) = yr*xr*charge;
                traced(5) = yr*yr*charge;
                traced(6) = yr*zr*charge;
                traced(7) = zr*xr*charge;
                traced(8) = zr*yr*charge;
                traced(9) = zr*zr*charge;
                qave = (traced(1) + traced(5) + traced(9)) / 3.0;
                sdhquadrupole(1) = sdhquadrupole(1)+1.5*(traced(1) - qave);
                sdhquadrupole(2) = sdhquadrupole(2)+1.5*(traced(2));
                sdhquadrupole(3) = sdhquadrupole(3)+1.5*(traced(3));
                sdhquadrupole(4) = sdhquadrupole(4)+1.5*(traced(4));
                sdhquadrupole(5) = sdhquadrupole(5)+1.5*(traced(5) - qave);
                sdhquadrupole(6)= sdhquadrupole(6)+1.5*(traced(6));
                sdhquadrupole(7) = sdhquadrupole(7)+ 1.5*(traced(7));
                sdhquadrupole(8) = sdhquadrupole(8)+1.5*(traced(8));
                sdhquadrupole(9) = sdhquadrupole(9)+1.5*(traced(9) - qave);
end

    ux = sdhdipole(1);
    uy = sdhdipole(2);
    uz = sdhdipole(3);
    qxx = sdhquadrupole(1) / 3.0;
    qxy = sdhquadrupole(2) / 3.0;
    qxz = sdhquadrupole(3) / 3.0;
    qyx = sdhquadrupole(4) / 3.0;
    qyy = sdhquadrupole(5) / 3.0;
    qyz = sdhquadrupole(6) / 3.0;
    qzx = sdhquadrupole(7) / 3.0;
    qzy = sdhquadrupole(8) / 3.0;
    qzz = sdhquadrupole(9) / 3.0;
    
    clear sdhquadrupole sdhdipole traced
    
   eps_r = dielw/dielp;
   tsr=zeros(3);
   
if bx==1
      fprintf('Processing boundary condition along x')
for j=1:dime(2)
    if mod (j,10)==0 
        fprintf('.')
    end
    for k=1:dime(3)
        for n=1:2
            i=BC(1,n);

    r=sqrt((X(i)-solutex)^2+(Y(j)-solutey)^2+(Z(k)-solutez)^2);         
    r2 = r*r;
    r3 = r2*r;
    r5 = r3*r2;
    tsr(1) = (1.0/dielw)/r;
    tsr(2) = (1.0/dielw)*(-1.0)/r3;
    tsr(3) = (1.0/dielw)*(3.0)/r5;
    if (xkappa < VSMALL) 
        tsr(2) = (3.0*eps_r)/(1.0 + 2.0*eps_r)*tsr(2);
        tsr(3) = (5.0*eps_r)/(2.0 + 3.0*eps_r)*tsr(3);
    else 
        ka = xkappa*soluterad;
        ka2 = ka*ka;
        ka3 = ka2*ka;
        kr = xkappa*r;
        kr2 = kr*kr;
        kr3 = kr2*kr;
        tsr(1) = exp(ka-kr) / (1.0 + ka) * tsr(1);
        tsr(2) = 3.0*eps_r*exp(ka-kr)*(1.0 + kr) * tsr(2);
        tsr(2) = tsr(2) / (1.0 + ka + eps_r*(2.0 + 2.0*ka + ka2));
        tsr(3) = 5.0*eps_r*exp(ka-kr)*(3.0 + 3.0*kr + kr2) * tsr(3);
        tsr(3) = tsr(3)/(6.0+6.0*ka+2.0*ka2+eps_r*(9.0+9.0*ka+4.0*ka2+ka3));
    end
                    val = pre1*sdhcharge*tsr(1);
                    val = val-pre1*ux*xr*tsr(2);
                    val= val- pre1*uy*yr*tsr(2);
                    val = val-pre1*uz*zr*tsr(2);
                    val = val+pre1*qxx*xr*xr*tsr(3);
                    val = val+pre1*qyy*yr*yr*tsr(3);
                    val = val+ pre1*qzz*zr*zr*tsr(3);
                    val = val+pre1*2.0*qxy*xr*yr*tsr(3);
                    val = val+pre1*2.0*qxz*xr*zr*tsr(3);
                    val = val+pre1*2.0*qyz*yr*zr*tsr(3);                           
           
            gxcf(j,k,n)=val;
            potB(i,j,k)=val;
        end
    end
end
fprintf('\nDone!\n')
end

%% Y Boundary
if by==1 
      fprintf('Processing boundary condition along y')
for i=1:dime(1)
    if mod (i,10)==0 
        fprintf('.')
    end
    for k=1:dime(3)
        for n=1:2
            j=BC(2,n);
            
    r=sqrt((X(i)-solutex)^2+(Y(j)-solutey)^2+(Z(k)-solutez)^2);         
    r2 = r*r;
    r3 = r2*r;
    r5 = r3*r2;
    tsr(1) = (1.0/dielw)/r;
    tsr(2) = (1.0/dielw)*(-1.0)/r3;
    tsr(3) = (1.0/dielw)*(3.0)/r5;
    if (xkappa < VSMALL) 
        tsr(2) = (3.0*eps_r)/(1.0 + 2.0*eps_r)*tsr(2);
        tsr(3) = (5.0*eps_r)/(2.0 + 3.0*eps_r)*tsr(3);
    else 
        ka = xkappa*soluterad;
        ka2 = ka*ka;
        ka3 = ka2*ka;
        kr = xkappa*r;
        kr2 = kr*kr;
        kr3 = kr2*kr;
        tsr(1) = exp(ka-kr) / (1.0 + ka) * tsr(1);
        tsr(2) = 3.0*eps_r*exp(ka-kr)*(1.0 + kr) * tsr(2);
        tsr(2) = tsr(2) / (1.0 + ka + eps_r*(2.0 + 2.0*ka + ka2));
        tsr(3) = 5.0*eps_r*exp(ka-kr)*(3.0 + 3.0*kr + kr2) * tsr(3);
        tsr(3) = tsr(3)/(6.0+6.0*ka+2.0*ka2+eps_r*(9.0+9.0*ka+4.0*ka2+ka3));
    end
                    val = pre1*sdhcharge*tsr(1);
                    val = val-pre1*ux*xr*tsr(2);
                    val= val- pre1*uy*yr*tsr(2);
                    val = val-pre1*uz*zr*tsr(2);
                    val = val+pre1*qxx*xr*xr*tsr(3);
                    val = val+pre1*qyy*yr*yr*tsr(3);
                    val = val+ pre1*qzz*zr*zr*tsr(3);
                    val = val+pre1*2.0*qxy*xr*yr*tsr(3);
                    val = val+pre1*2.0*qxz*xr*zr*tsr(3);
                    val = val+pre1*2.0*qyz*yr*zr*tsr(3);                           
                    
            gycf(i,k,n)=val;
            potB(i,j,k)=val;
        end
    end
end
fprintf('\nDone!\n')
end

%% Z Boundary
if bz==1
      fprintf('Processing boundary condition along z')
for i=1:dime(1)
    if mod (i,10)==0 
        fprintf('.')
    end
    for j=1:dime(2)
        for n=1:2
            k=BC(3,n);

    r=sqrt((X(i)-solutex)^2+(Y(j)-solutey)^2+(Z(k)-solutez)^2);         
    r2 = r*r;
    r3 = r2*r;
    r5 = r3*r2;
    tsr(1) = (1.0/dielw)/r;
    tsr(2) = (1.0/dielw)*(-1.0)/r3;
    tsr(3) = (1.0/dielw)*(3.0)/r5;
    if (xkappa < VSMALL) 
        tsr(2) = (3.0*eps_r)/(1.0 + 2.0*eps_r)*tsr(2);
        tsr(3) = (5.0*eps_r)/(2.0 + 3.0*eps_r)*tsr(3);
    else 
        ka = xkappa*soluterad;
        ka2 = ka*ka;
        ka3 = ka2*ka;
        kr = xkappa*r;
        kr2 = kr*kr;
        kr3 = kr2*kr;
        tsr(1) = exp(ka-kr) / (1.0 + ka) * tsr(1);
        tsr(2) = 3.0*eps_r*exp(ka-kr)*(1.0 + kr) * tsr(2);
        tsr(2) = tsr(2) / (1.0 + ka + eps_r*(2.0 + 2.0*ka + ka2));
        tsr(3) = 5.0*eps_r*exp(ka-kr)*(3.0 + 3.0*kr + kr2) * tsr(3);
        tsr(3) = tsr(3)/(6.0+6.0*ka+2.0*ka2+eps_r*(9.0+9.0*ka+4.0*ka2+ka3));
    end
                    val = pre1*sdhcharge*tsr(1);
                    val = val-pre1*ux*xr*tsr(2);
                    val= val- pre1*uy*yr*tsr(2);
                    val = val-pre1*uz*zr*tsr(2);
                    val = val+pre1*qxx*xr*xr*tsr(3);
                    val = val+pre1*qyy*yr*yr*tsr(3);
                    val = val+ pre1*qzz*zr*zr*tsr(3);
                    val = val+pre1*2.0*qxy*xr*yr*tsr(3);
                    val = val+pre1*2.0*qxz*xr*zr*tsr(3);
                    val = val+pre1*2.0*qyz*yr*zr*tsr(3);                           
                    
            gzcf(i,j,n)=val;
            potB(i,j,k)=val;
        end
    end
end
fprintf('\nDone!\n')
end

clear X Y Z tsr

