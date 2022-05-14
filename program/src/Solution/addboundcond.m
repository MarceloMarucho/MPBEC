potc=zeros(dime(1),dime(2),dime(3));

for i=2:dime(1)-1
    for j=2:dime(2)-1
        for k=2:dime(3)-1
             pe=(k-2)*(dime(1)-2)*(dime(2)-2)+(j-2)*(dime(1)-2)+i-1;
             potc(i,j,k)=pote(pe);
        end
    end
end
 
% solution

% adding Boundary condition
% if Dirichlet along the three directions
%if bx==1 && by==1 && bz==1
if (strcmp(bc, 'sdh')==1 ||strcmp(oldbc, 'focusname.inm')==1 || strcmp(bc, 'mdh')==1)
MATLAB_pot=potB;
end
% Periodic along x and y and Dirichlet along z
%if bx ==2 && by==2 && bz==1
if strcmp(bc, 'mixed')==1
%    MATLAB_pot(2:dime(1)-1,dime(2),dime(3))=potB(2:dime(1)-1,dime(2),dime(3));
 %   value of the nodes over the x-faces from periodicity symmetry
 for j=2: dime(2)-1
     for k=2:dime(3)-1
         MATLAB_pot(1,j,k)=potc(dime(1)-1,j,k);
         MATLAB_pot(dime(1),j,k)=potc(2,j,k);
     end
 end
 %   value of the nodes boundary over the y-faces from periodicity symmetry
 for i=2: dime(1)-1
     for k=2:dime(3)-1
         MATLAB_pot(i,1,k)=potc(i,dime(2)-1,k);
         MATLAB_pot(i,dime(2),k)=potc(i,2,k);
     end
 end
  %   value of the nodes over the z-faces from Hucke Debye values (including
  %   the corners and edges defining xz and yz interfaces)
 for i=1: dime(1)
     for j=1:dime(2)
         MATLAB_pot(i,j,1)=potB(i,j,1);
         MATLAB_pot(i,j,dime(3))=potB(i,j,dime(3));
     end
 end
 %  values of the nodes along the edges defined along z direction from continuity
   for k=2:dime(3)-1
    MATLAB_pot(1,1,k)=(MATLAB_pot(2,1,k)+MATLAB_pot(1,2,k))*0.5;
    MATLAB_pot(1,dime(2),k)=(MATLAB_pot(2,dime(2),k)+MATLAB_pot(1,dime(2)-1,k))*0.5; 
    MATLAB_pot(dime(1),1,k)=(MATLAB_pot(dime(1)-1,1,k)+MATLAB_pot(dime(1),2,k))*0.5;
    MATLAB_pot(dime(1),dime(2),k)=(MATLAB_pot(dime(1)-1,dime(2),k)+MATLAB_pot(dime(1),dime(2)-1,k))*0.5; 
   end
end

% Periodic along the three directions
%if bx ==2 && by==2 && bz==2
if strcmp(bc, 'periodic')==1
%    MATLAB_pot(2:dime(1)-1,dime(2),dime(3))=potB(2:dime(1)-1,dime(2),dime(3));
 %   value of the nodes over the x-faces from periodicity symmetry
 for j=2: dime(2)-1
     for k=2:dime(3)-1
         MATLAB_pot(1,j,k)=potc(dime(1)-1,j,k);
         MATLAB_pot(dime(1),j,k)=potc(2,j,k);
     end
 end
 %   value of the nodes boundary over the y-faces from periodicity symmetry
 for i=2: dime(1)-1
     for k=2:dime(3)-1
         MATLAB_pot(i,1,k)=potc(i,dime(2)-1,k);
         MATLAB_pot(i,dime(2),k)=potc(i,2,k);
     end
 end
  %   value of the nodes over the z-face from periodicity
 for i=1: dime(1)-1
     for j=1:dime(2)-1
         MATLAB_pot(i,j,1)=potc(i,j,dime(3)-1);
         MATLAB_pot(i,j,dime(3))=potc(i,j,2);
     end
 end
 %  values of the nodes along the edges defined along z direction from continuity
   for k=2:dime(3)-1
    MATLAB_pot(1,1,k)=(MATLAB_pot(2,1,k)+MATLAB_pot(1,2,k))*0.5;
    MATLAB_pot(1,dime(2),k)=(potc(2,dime(2),k)+potc(1,dime(2)-1,k))*0.5; 
    MATLAB_pot(dime(1),1,k)=(potc(dime(1)-1,1,k)+potc(dime(1),2,k))*0.5;
    MATLAB_pot(dime(1),dime(2),k)=(potc(dime(1)-1,dime(2),k)+potc(dime(1),dime(2)-1,k))*0.5; 
   end
%  values of the nodes along the edges defined along x direction from continuity
     for i=2:dime(1)-1
    MATLAB_pot(i,dime(2),1)=(potc(i,dime(2)-1,1)+potc(i,dime(2),2))*0.5;
    MATLAB_pot(i,dime(2),dime(3))=(potc(i,dime(2)-1,dime(3))+potc(i,dime(2),dime(3)))*0.5; 
    MATLAB_pot(i,1,1)=(potc(i,2,1)+potc(i,1,2))*0.5; 
    MATLAB_pot(i,1,dime(3))=(potc(i,2,dime(3))+potc(i,1,dime(3)-1))*0.5; 
     end
%  values of the nodes along the edges defined along y direction from continuity
     for j=2:dime(2)-1
    MATLAB_pot(1,j,dime(3))=(potc(2,j,dime(3))+potc(1,j,dime(3)-1))*0.5;
    MATLAB_pot(dime(1),j,dime(3))=(potc(dime(1)-1,j,dime(3))+potc(dime(1),j,dime(3)-1))*0.5; 
    MATLAB_pot(dime(1),j,1)=(potc(dime(1)-1,j,1)+potc(dime(1),j,2))*0.5; 
    MATLAB_pot(1,j,1)=(potc(2,j,1)+potc(1,j,2))*0.5; 
     end
     % values of nodes at the corners from coontinuity criteria
 MATLAB_pot(1,1,1)=(MATLAB_pot(2,1,1)+MATLAB_pot(1,2,1)+MATLAB_pot(1,1,2))/3.0;  
 MATLAB_pot(1,dime(2),1)=(MATLAB_pot(2,dime(2),1)+MATLAB_pot(1,dime(2)-1,1)+MATLAB_pot(1,dime(2),2))/3.0; 
 MATLAB_pot(dime(1),1,1)=(MATLAB_pot(dime(1)-1,1,1)+MATLAB_pot(dime(1),2,1)+MATLAB_pot(dime(1),1,2))/3.0; 
 MATLAB_pot(1,1,dime(3))=(MATLAB_pot(2,1,dime(3))+MATLAB_pot(1,2,dime(3))+MATLAB_pot(1,1,dime(3)-1))/3.0; 
 MATLAB_pot(1,dime(2),dime(3))=(MATLAB_pot(2,dime(2),dime(3))+MATLAB_pot(1,dime(2)-1,dime(3))+MATLAB_pot(1,dime(2),dime(3)-1))/3.0; 
 MATLAB_pot(dime(1),1,dime(3))=(MATLAB_pot(dime(1)-1,1,dime(3))+MATLAB_pot(dime(1),2,dime(3))+MATLAB_pot(dime(1),1,dime(3)-1))/3.0; 
 MATLAB_pot(dime(1),dime(2),1)=(MATLAB_pot(dime(1)-1,dime(2),1)+MATLAB_pot(dime(1),dime(2)-1,1)+MATLAB_pot(dime(1),dime(2),2))/3.0; 
 MATLAB_pot(dime(1),dime(2),dime(3))=(MATLAB_pot(dime(1)-1,dime(2),dime(3))+MATLAB_pot(dime(1),dime(2)-1,dime(3))+MATLAB_pot(dime(1),dime(2),dime(3)-1))/3.0; 
end
% adding the solution of the interior nodes from the solution of the linear
% equation
MATLAB_pot(2:dime(1)-1,2:dime(2)-1,2:dime(3)-1)=potc(2:dime(1)-1,2:dime(2)-1,2:dime(3)-1);