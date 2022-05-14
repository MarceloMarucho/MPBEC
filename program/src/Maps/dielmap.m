
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dielectric map calculation for surfmeth= MOl (discontinuos functions using
% collony and inflated van der Waals surface area definition)
% xmin, ymin, zmin, xmax, ymax, zmax are already calculated at discretization script
% nx. ny and nz are equal to dim vector from input file
% hx,hy and hzed are equat to h vector calcualted at MAPBS script

% allocate and assign solvent dielectric values to all the space

epsx=ones(dime(1),dime(2),dime(3))*dielw;
epsy=epsx;
epsz=epsx;

%memory_in_use = monitor_memory_whos
%if (memory_in_use>500)
%    disp('running out of memory')
%end

if (dielw==dielp) % reference state
    
 dielx=epsx;
 diely=epsy;
 dielz=epsz;

 return

end

% Because the space within the solute the dielectric constant is actually equal to the solute value
% the following loop trought the atoms of the molecule determines the corresponding set of points 
% which are inside the solvente-inflated van der waals radii and consequently where the
% dielectric value needs to be replaced by the solute value (markSphere).


  %SRAD IS THE SOLVENT MOLECULE RADIUS (OFTEN 1.4 a FOR WATER) AND SHOULD
  %BE READ FROM IMPUT FILE!!!! dielp is the solute dielectric value. 
  %TO EXTEND THE READ SCRIPT FOR THESE two new parameters
  % the inflated van der waals dadii reads
  
            infl=atomR+srad; %used by markSphere
            infl2=infl.^2;
            
            % position of the atom
            
            idpos=atomP; %used by markSphere
            
            % reset the solute dielectric value on each direction
            
            markval=dielp; %  Used by markSphere
            
            % mark x-shifted dielectric (shift the xmin)
                        
            xminold=xmin;
            xmin=xminold+0.5*h(1); %used by markSphere
            epsx=markSpheref_vanderwalls(epsx,idpos,xmin,ymin,zmin,xmax,ymax,zmax,infl,infl2,h,markval,dime,atomR,VSMALL);
            xmin=xminold;
            
            % mark y-shifted dielectric (shift the ymin)
            
            yminold=ymin;
            ymin=yminold+0.5*h(2); %used by markSphere
            epsy=markSpheref_vanderwalls(epsy,idpos,xmin,ymin,zmin,xmax,ymax,zmax,infl,infl2,h,markval,dime,atomR,VSMALL);
            ymin=yminold;

            % mark z-shifted dielectric (shift the zmin)
            
            zminold=zmin;
            zmin=zminold+0.5*h(3); %used by marksphere
            epsz=markSpheref_vanderwalls(epsz,idpos,xmin,ymin,zmin,xmax,ymax,zmax,infl,infl2,h,markval,dime,atomR,VSMALL);
            zmin=zminold;


clear idpos infl infl2
% for non-zero solvent radii the loops must also run over the solvent
% accessible surface points.

if (srad>VSMALL) 
    %lets check if we need the cell list approach
    
    if (atomN < 20)
        % we don't really need it
        ese3=which('Vacc_atomSurfless20');
        copyfile(ese3,'Vacc_atomSurf.m');
    else
        ese2=which('Vacc_atomSurfplus20');
    copyfile(ese2,'Vacc_atomSurf.m');
    
    run getcell % will provide the dimension of the cell and the list of  
%the position of the atoms "position" in the cell and the number of atoms
% "natomcell". It also provide the maximun probe RADIUS THAT THE CELL WAS BUILT WITH
    end
    
% overall calculation
% The atomic solvent accessible surface area (SASA) is evaluated for each
% atom in the molecule (Vacc_atomSurf). If the corresponding area in
% greater than zero, the set of points for this atom's SASA is determined
% and saved in a array (poson). For each of the points in such array the solvent
% accessibility is reset (markSphere).

prad=maxradprobe; 
mrad=maxradatom+maxradprobe; % radii equal to the sum of the atomic van der waals radius and the probe radius
% lets generate "refsphenpts" points "xpts,ypts, and zpts" for the reference
% sphere of unit radius. 
Maxarea=4.0*pi*mrad^2;
refshenpts=ceil(Maxarea*double(surf_density));
disp(' ')
run Vacc_SurfSphere % output xpts,ypts,zpts and refsphenpts
% Lets determine which points will contribute to the accessible surface
disp(['Constructing reference sphere with Max area =  ',num2str(Maxarea), '  and number of points  ',num2str(refshenpts)])
  
%    run Vacc_atomSurf  % this routine returns the "area", and "poson"
%    array and "numm" 
disp(' ')
run Vacc_atomSurf
    % calculate the number of points I have to correct
       
    % now the inflated radius is just the solvent radii
    
            infl=srad; %used by marksphere
            infl2=infl^2;
            
            %reset the solvent dielectric value on each direction
            
            markval=dielw; %used by marksphere
            
          
            % mark x-shifted dielectric (shift the xmin)
 
            xminold=xmin;
            xmin=xminold+0.5*h(1); %used by marksphere
            epsx=markSpheref_collony(epsx,poson,xmin,ymin,zmin,infl,infl2,h,markval,dime,area,numm);
            xmin=xminold;
 
            % mark y-shifted dielectric (shift the ymin)
 
            yminold=ymin;
            ymin=yminold+0.5*h(2); %used by marksphere
            epsy=markSpheref_collony(epsy,poson,xmin,ymin,zmin,infl,infl2,h,markval,dime,area,numm);
            ymin=yminold;
 
            % mark z-shifted dielectric (shift the zmin)
 
            zminold=zmin;
            zmin=zminold+0.5*h(3);  %used by marksphere
            epsz=markSpheref_collony(epsz,poson,xmin,ymin,zmin,infl,infl2,h,markval,dime,area,numm);
            zmin=zminold;
 
delete ('Vacc_atomSurf.m')
end

%end
clear idpos poson numm area

dielx=epsx;
diely=epsy;
dielz=epsz;


