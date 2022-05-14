% estimate maximum allocation memory required by MPBEC to perform the
% calculations

% read available RAM memory
availmemm=available_RAM_memory;
if (availmemm==-1) 
    return
end
% gradually increase the allocated memory

fc3=zeros(dime(1));
memory_in_use1 = monitor_MPBEC_used_memory;
memory_grid1=int16(memory_in_use1);
% first partial allocation
partmem=memory_grid1*16;
clear fc3
run parameters
% comparisoin between available and required RAM memory. IF running out of
% memory, alert the following message
if (partmem>0.75*availmemm)
    disp('GRID')
disp(' ')
disp(['number of points=   ',num2str(dime),',    box size (in A)=  ',num2str(glen), ',    Number of atoms=  ',num2str(atomN)])
disp(' ')
    disp('The available memory of your computer is insufficient to perform the required calculations')
    
argg2=sprintf('%s %g %s','The available RAM memory in your computer (',availmemm,' Mbs) seems to be insufficient to run MPBEC. Please consider to abort the calculations. You may need to close unused programs, free cached memory, and/or use the MPBEC advanced user level mode to manually set lower number of grid points, and try it again. If the problem persists, you would need to increase the RAM memory of your computer. Do you want to abort the calculations?? ');
button= questdlg(argg2, ...
        'Exit Dialog','abort','continue','abort');
         switch button
            case 'abort',
              disp('');
              disp('Exiting MPBEC calculation. This may take a few seconds...');
              exitmpbec=1 ;
    return
              stop
              pause(5);
            case 'continue',
                disp('');
               disp('resume calculations...');   
         end
end

% increasing allocation
fc1=zeros(dime(1),dime(2));
memory_in_use1 = monitor_MPBEC_used_memory;
memory_grid1=int16(memory_in_use1);
partmem=memory_grid1*16;
clear fc1

if (partmem>0.75*availmemm)
    disp('GRID')
disp(' ')
disp(['number of points=   ',num2str(dime),',    box size (in A)=  ',num2str(glen), ',    Number of atoms=  ',num2str(atomN)])
disp(' ')
    disp('The available memory of your computer is insufficient to perform the required calculations')
    
argg2=sprintf('%s %g %s','The available RAM memory in your computer (',availmemm,' Mbs) seems to be insufficient to run MPBEC. Please consider to abort the calculations. You may need to close unused programs, free cached memory, and/or use the MPBEC advanced user level mode to manually set lower number of grid points, and try it again. If the problem persists, you would need to increase the RAM memory of your computer. Do you want to abort the calculations?? ');
button= questdlg(argg2, ...
        'Exit Dialog','abort','continue','abort');
         switch button
            case 'abort',
              disp('');
              disp('Exiting MPBEC calculation. This may take a few seconds...');
              exitmpbec=1 ;
    return
              stop
              pause(5);
            case 'continue',
                disp('');
               disp('resume calculations...');   
         end
end

% increasing allocation
fc=zeros(dime(1)-1,dime(2)-1,dime(3)-1);
memory_in_use1 = monitor_MPBEC_used_memory;
memory_grid1=int16(memory_in_use1);
partmem=memory_grid1*16;
clear fc

if (partmem>0.75*availmemm)
    disp('GRID')
disp(' ')
disp(['number of points=   ',num2str(dime),',    box size (in A)=  ',num2str(glen), ',    Number of atoms=  ',num2str(atomN)])
disp(' ')
    disp('The available memory of your computer is insufficient to perform the required calculations')
    
argg2=sprintf('%s %g %s','The available RAM memory in your computer (',availmemm,' Mbs) seems to be insufficient to run MPBEC. Please consider to abort the calculations. You may need to close unused programs, free cached memory, and/or use the MPBEC advanced user level mode to manually set lower number of grid points, and try it again. If the problem persists, you would need to increase the RAM memory of your computer. Do you want to abort the calculations?? ');
button= questdlg(argg2, ...
        'Exit Dialog','abort','continue','abort');
         switch button
            case 'abort',
              disp('');
              disp('Exiting MPBEC calculation. This may take a few seconds...');
              exitmpbec=1 ;
    return
              stop
              pause(5);
            case 'continue',
                disp('');
               disp('resume calculations...');   
         end
end

% increasing allocation
AAMain=zeros(prod(dime-2),1);
memory_in_use2 = monitor_MPBEC_used_memory;
memory_grid2=int16(memory_in_use2);
memory_grid2=memory_grid2*6;
clear AAMain
partmem=partmem+memory_grid2;
if (partmem>0.75*availmemm)
    disp('GRID')
disp(' ')
disp(['number of points=   ',num2str(dime),',    box size (in A)=  ',num2str(glen), ',    Number of atoms=  ',num2str(atomN)])
disp(' ')
    disp('The available memory of your computer is insufficient to perform the required calculations')
    
argg2=sprintf('%s %g %s','The available RAM memory in your computer (',availmemm,' Mbs) seems to be insufficient to run MPBEC. Please consider to abort the calculations. You may need to close unused programs, free cached memory, and/or use the MPBEC advanced user level mode to manually set lower number of grid points, and try it again. If the problem persists, you would need to increase the RAM memory of your computer. Do you want to abort the calculations?? ');
button= questdlg(argg2, ...
        'Exit Dialog','abort','continue','abort');
         switch button
            case 'abort',
              disp('');
              disp('Exiting MPBEC calculation. This may take a few seconds...');
              exitmpbec=1 ;
    return
            case 'continue',
                disp('');
               disp('resume calculations...');   
         end

end

% increasing allocation
AAMain=zeros(prod(dime-2),1);
diagonals=horzcat(AAMain,AAMain,AAMain,AAMain);
memory_in_use3 = monitor_MPBEC_used_memory;
memory_grid3=int16(memory_in_use3);
partmem=partmem+memory_grid3;
if (partmem>0.75*availmemm)
    disp('GRID')
disp(' ')
disp(['number of points=   ',num2str(dime),',    box size (in A)=  ',num2str(glen), ',    Number of atoms=  ',num2str(atomN)])
disp(' ')
     disp('The available memory of your computer is insufficient to perform the required calculations')
    
argg2=sprintf('%s %g %s','The available RAM memory in your computer (',availmemm,' Mbs) seems to be insufficient to run MPBEC. Please consider to abort the calculations. You may need to close unused programs, free cached memory, and/or use the MPBEC advanced user level mode to manually set lower number of grid points, and try it again. If the problem persists, you would need to increase the RAM memory of your computer. Do you want to abort the calculations?? ');
button= questdlg(argg2, ...
        'Exit Dialog','abort','continue','abort');
         switch button
            case 'abort',
              disp('');
              disp('Exiting MPBEC calculation. This may take a few seconds...');
              exitmpbec=1 ;
return
            case 'continue',
                disp('');
               disp('resume calculations...');   
         end

end
clear AAMain diagonals

% increasing allocation
dcdiag=0;
dcx=-1;
dcy=-(dime(1)-2);
dcz=-(dime(1)-2)*(dime(2)-2);
pcx=-(dime(1)-3);
pcy=-((dime(1)-2)*(dime(2)-3));
pcz=-((dime(1)-2)*(dime(2)-2)*(dime(3)-3));
d=[dcdiag dcx dcy dcz];
AAMain=zeros(prod(dime-2),1);
diagonals=horzcat(AAMain,AAMain,AAMain,AAMain);
AA=spdiags(diagonals, d, prod(dime-2), prod(dime-2));
memory_in_use4 = monitor_MPBEC_used_memory;
memory_grid4=int16(memory_in_use4);
memory_grid4=memory_grid4*4;
clear AA AAMain diagonals
partmem=partmem+memory_grid4;
if (partmem>0.75*availmemm)
    disp('GRID')
disp(' ')
disp(['number of points=   ',num2str(dime),',    box size (in A)=  ',num2str(glen), ',    Number of atoms=  ',num2str(atomN)])
disp(' ')
        disp('The available memory of your computer is insufficient to perform the required calculations')
    
argg2=sprintf('%s %g %s','The available RAM memory in your computer (',availmemm,' Mbs) seems to be insufficient to run MPBEC. Please consider to abort the calculations. You may need to close unused programs, free cached memory, and/or use the MPBEC advanced user level mode to manually set lower number of grid points, and try it again. If the problem persists, you would need to increase the RAM memory of your computer. Do you want to abort the calculations?? ');
button= questdlg(argg2, ...
        'Exit Dialog','abort','continue','abort');
         switch button
            case 'abort',
              disp('');
              disp('Exiting MPBEC calculation. This may take a few seconds...');
              exitmpbec=1 ;
return
            case 'continue',
                disp('');
               disp('resume calculations...');   
         end

end

% increasing allocation

X=linspace(-glen(1)/2,glen(1)/2,dime(1));
Y=linspace(-glen(2)/2,glen(2)/2,dime(2));
Z=linspace(-glen(3)/2,glen(3)/2,dime(3));

area=zeros(atomN,1);
memory_in_use5 = monitor_MPBEC_used_memory;
memory_grid5=int16(memory_in_use5);
memory_grid5=memory_grid5*2;
clear area X Y Z
partmem=partmem+memory_grid5;
% total estimation on allocated memory
if (partmem>0.75*availmemm)
    disp('GRID')
disp(' ')
disp(['number of points=   ',num2str(dime),',    box size (in A)=  ',num2str(glen), ',    Number of atoms=  ',num2str(atomN)])
disp(' ')
    disp('The available memory of your computer is insufficient to perform the required calculations')
    
argg2=sprintf('%s %g %s','The available RAM memory in your computer (',availmemm,' Mbs) seems to be insufficient to run MPBEC. Please consider to abort the calculations. You may need to close unused programs, free cached memory, and/or use the MPBEC advanced user level mode to manually set lower number of grid points, and try it again. If the problem persists, you would need to increase the RAM memory of your computer. Do you want to abort the calculations?? ');
button= questdlg(argg2, ...
        'Exit Dialog','abort','continue','abort');
         switch button
            case 'abort',
              disp('');
              disp('Exiting MPBEC calculation. This may take a few seconds...');
              exitmpbec=1 ;
 return
            case 'continue',
                disp('');
               disp('resume calculations...');   
         end

end

% if check memory is OK, provide the following information to the user
argg2=sprintf('%s %g %s','The available RAM memory in your computer is (',availmemm,' Mbs) and the minimum memory required by MPBEC to perform');
argg3=sprintf('%s %g %s','the electrostatic calculations is estimated around (',partmem,' Mbs).');    
disp(argg2);
disp(argg3);
disp(' ');
disp('You should not have memory allocation issues. The computing time usually depends on your computer processor speed, ');
disp('the size of the biological system under consideration, the grid resolution and the required calculations');
disp(' ');  

 
