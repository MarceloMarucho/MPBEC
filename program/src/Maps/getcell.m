% lets built the cell list
%%%%%%%%%%%%%%%%%%%%
% %get molecular dimensions
lower_corner=zeros(VAPBS_DIM,1);
upper_corner=lower_corner;
for i=0:VAPBS_DIM-1
    lower_corner(i+1)=VLARGE;
    upper_corner(i+1)=-VLARGE;
end
r_max=-1.0;
% check each atom
for i=0:atomN-1
    for j=0:VAPBS_DIM-1
        post=atomP(i+1,j+1);
        if (post<lower_corner(j+1)) 
            lower_corner(j+1)=post;
        end
        if (post>upper_corner(j+1)) 
            upper_corner(j+1)=post;
        end
    end
        if (atomR(i+1)>r_max)
            r_max=atomR(i+1);
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert to grid based coordinates (pos=position to test, see Vacc_atomSurf)

% set up grid spacing
deltagrid=VCLIST_INFLATE*(r_max+maxionR+0.5); % check this
upper_cornerh=upper_corner+deltagrid;
lower_cornerh=lower_corner-deltagrid;
lengthh=upper_cornerh-lower_cornerh;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cell list constructor

Lc=zeros(VAPBS_DIM,1);
for i=1:VAPBS_DIM
%    Lc(i)=floor(2.*solutelen(i));
    Lc(i)=floor(lengthh(i)/r_max/2.);
if (Lc(i)<= 3)
    Lc(i)=3;
end
if (Lc(i)> MAXIHTABLEDIM)
    Lc(i)=MAXIHTABLEDIM;
end
end

disp(' ')
disp(['Constructing Cell list with ',num2str(Lc(1)),',',num2str(Lc(2)),',',num2str(Lc(3)), '  table'])
disp(' ')
disp(['Grid lengths =  ',num2str(lengthh(1)),',',num2str(lengthh(2)),',',num2str(lengthh(3))])
disp(' ')
disp(['Grid lower corner =  ',num2str(lower_cornerh(1)),',',num2str(lower_cornerh(2)),',',num2str(lower_cornerh(3))])
disp(' ')

% reset headers
head=empty*ones(Lc(1)+1,Lc(2)+1,Lc(3)+1);
ll=zeros(atomN,1);
poscor=zeros(atomN,3);
rn=zeros(3,1);
mc=zeros(atomN,3);
% scan atoms to contruct headers, list, and linker list
    for i=0:atomN-1
        % for each atom let obtain the cell index mc(1),mc(2),mc(3)
        for a=0:2
        poscor(i+1,a+1)=atomP(i+1,a+1)-lower_cornerh(a+1);
         rn(a+1)=lengthh(a+1)/double(Lc(a+1));
        mc(i+1,a+1)=floor(poscor(i+1,a+1)/rn(a+1));
        end
        % atom counter
        cx=mc(i+1,1)+1;
        cy=mc(i+1,2)+1;
        cz=mc(i+1,3)+1;
        %linked list previous ocuupant
        ll(i+1)=head(cx,cy,cz);
        % the last one goes to the header
        head(cx,cy,cz)=i+1;
    end
    clear mc poscor
   for cx=1:Lc(1)+1 
   for cy=1:Lc(2)+1
       for cz=1:Lc(3)+1
           if (cx==(Lc(1)+1) || cy==(Lc(2)+1) || cz==(Lc(3)+1))
    head(cx,cy,cz)=-1;
           end
       end
   end
   end
