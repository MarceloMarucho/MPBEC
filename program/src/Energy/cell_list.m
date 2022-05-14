% lets built the cell list
%%%%%%%%%%%%%%%%%%%%
% %get molecular dimensions
lowerCorner=zeros(VAPBS_DIM,1);
upperCorner=lowerCorner;
for i=0:VAPBS_DIM-1
    lowerCorner(i+1)=VLARGE;
    upperCorner(i+1)=-VLARGE;
end
r_max=-1.0;
% check each atom
for i=0:atomN-1
    for j=0:VAPBS_DIM-1
        post=atomP(i+1,j+1);
        if (post<lowerCorner(j+1)) 
            lowerCorner(j+1)=post;
        end
        if (post>upperCorner(j+1)) 
            upperCorner(j+1)=post;
        end
    end
        if (atomR(i+1)>r_max)
            r_max=atomR(i+1);
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert to grid based coordinates (pos=position to test, see Vacc_atomSurf)

% set up grid spacing
% deltagrid=VCLIST_INFLATE*(r_max+maxionR+0.5); % check this
% upperCorner=upperCorner-h/2.;
% lowerCorner=lowerCorner+h/2.;

 upperCorner=[xmax ymax zmax];
 lowerCorner=[xmin ymin zmin];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    xok = 0;
    yok = 0;
    zok = 0;
    sdime12=dime(1)*dime(2);
    VAPBS_LEFT=4;
    VAPBS_RIGHT=1;
    VAPBS_DOWN=6;
    VAPBS_UP=3;
    VAPBS_FRONT=2;
    VAPBS_BACK=5;

bflags=zeros(6,1);
partID=zeros(atomN,1);

    for i=1:atomN
   
        if ((atomP(i,1) < upperCorner(1)) && ...
            ( atomP(i,1)> lowerCorner(1)))
        xok = 1;
        else 
            if ((abs(atomP(i,1) - lowerCorner(1)) < VPMGSMALL) && ... 
                (bflags(VAPBS_LEFT) == 0)) 
            xok = 1;
            elseif ((abs(atomP(i,1) - lowerCorner(1)) < VPMGSMALL) && ... 
                (bflags(VAPBS_LEFT) == 1))  
            xok = 0.5;
            elseif ((abs(atomP(i,1) - lowerCorner(1)) < VPMGSMALL) && ... 
                (bflags(VAPBS_RIGHT) == 0)) 
            xok = 1;
            elseif ((abs(atomP(i,1) - lowerCorner(1)) < VPMGSMALL) && ... 
                (bflags(VAPBS_RIGHT) == 1)) 
            xok = 0.5;
            else
                xok = 0;
            end
        end
        if ((atomP(i,2) < upperCorner(2)) && ...
            ( atomP(i,2)> lowerCorner(2)))
        yok = 1;
        else 
            if ((abs(atomP(i,2) - lowerCorner(2)) < VPMGSMALL) && ... 
                (bflags(VAPBS_BACK) == 0)) 
            yok = 1;
            elseif ((abs(atomP(i,2) - lowerCorner(2)) < VPMGSMALL) && ... 
                (bflags(VAPBS_BACK) == 1))  
            yok = 0.5;
            elseif ((abs(atomP(i,2) - lowerCorner(2)) < VPMGSMALL) && ... 
                (bflags(VAPBS_FRONT) == 0)) 
            yok = 1;
            elseif ((abs(atomP(i,2) - lowerCorner(2)) < VPMGSMALL) && ... 
                (bflags(VAPBS_FRONT) == 1)) 
            yok = 0.5;
            else
                yok = 0;
            end
        end
        if ((atomP(i,3) < upperCorner(3)) && ...
            ( atomP(i,3)> lowerCorner(3)))
        zok = 1;
        else 
            if ((abs(atomP(i,3) - lowerCorner(3)) < VPMGSMALL) && ... 
                (bflags(VAPBS_DOWN) == 0)) 
            zok = 1;
            elseif ((abs(atomP(i,3) - lowerCorner(3)) < VPMGSMALL) && ... 
                (bflags(VAPBS_DOWN) == 1))  
            zok = 0.5;
            elseif ((abs(atomP(i,3) - lowerCorner(3)) < VPMGSMALL) && ... 
                (bflags(VAPBS_UP) == 0)) 
            zok = 1;
            elseif ((abs(atomP(i,3) - lowerCorner(3)) < VPMGSMALL) && ... 
                (bflags(VAPBS_UP) == 1)) 
            zok = 0.5;
            else
                zok = 0;
            end
        end

        partID(i) = xok*yok*zok; 

    end
   
   pvec=zeros(dime(1),dime(2),dime(3));
   
    for i=1:dime(1)
 %       xok = 0.0;
        x = (i-1)*h(1) + xmin;
        if ( (x < (upperCorner(1)-h(1)/2.)) && ...
             (x > (lowerCorner(1)+h(1)/2.))) 
           xok = 1.0;
        elseif ( (abs(x - lowerCorner(1)) < VPMGSMALL) && ...
                  (bflags(VAPBS_LEFT) == 0)) 
              xok = 1.0;
        elseif ((abs(x - lowerCorner(1)) < VPMGSMALL) && ...
                 (bflags(VAPBS_LEFT) == 1)) 
             xok = 0.5;
        elseif ((abs(x - upperCorner(1)) < VPMGSMALL) && ...
                 (bflags(VAPBS_RIGHT) == 0)) 
             xok = 1.0;
        elseif ((abs(x - upperCorner(1)) < VPMGSMALL) && ...
                 (bflags(VAPBS_RIGHT) == 1)) 
             xok = 0.5;
        elseif ((x > (upperCorner(1) + h(1)/2.)) || (x < (lowerCorner(1) - h(1)/2))) 
                xok = 0.0;
        elseif ((x < (upperCorner(1) + h(1)/2.)) || (x > (lowerCorner(1) - h(1)/2))) 
            x0 = max(x - h(1)/2., lowerCorner(1));
            x1 = min(x + h(1)/2., upperCorner(1));
            xok = abs(x1-x0)/h(1);

            if (xok < 0.0) 
                if (abs(xok) < VPMGSMALL) 
                    xok = 0.0;
                end
            end
            if (xok > 1.0) 
                if (abs(xok - 1.0) < VPMGSMALL) 
                    xok = 1.0;
                end 
            end
        %end
          else
                xok = 0.0;
        end
     
        for j=1: dime(2)
 %          yok = 0.0;
            y = (j-1)*h(2) + ymin;
        if ( (y < (upperCorner(2)-h(2)/2.)) && ...
             (y > (lowerCorner(2)+h(2)/2.))) 
           yok = 1.0;
        elseif ( (abs(y - lowerCorner(2)) < VPMGSMALL) && ...
                  (bflags(VAPBS_LEFT) == 0)) 
              yok = 1.0;
        elseif ((abs(y - lowerCorner(2)) < VPMGSMALL) && ...
                 (bflags(VAPBS_LEFT) == 1)) 
             yok = 0.5;
        elseif ((abs(y - upperCorner(2)) < VPMGSMALL) && ...
                 (bflags(VAPBS_RIGHT) == 0)) 
             yok = 1.0;
        elseif ((abs(y - upperCorner(2)) < VPMGSMALL) && ...
                 (bflags(VAPBS_RIGHT) == 1)) 
             yok = 0.5;
        elseif ((y > (upperCorner(2) + h(2)/2.)) || (y < (lowerCorner(2) - h(2)/2.))) 
                yok = 0.0;
        elseif ((y < (upperCorner(2) + h(2)/2.)) || (y > (lowerCorner(2) - h(2)/2.))) 
            y0 = max(y - h(2)/2., lowerCorner(2));
            y1 = min(y + h(2)/2., upperCorner(2));
            yok = abs(y1-y0)/h(2);

            if (yok < 0.0) 
                if (abs(yok) < VPMGSMALL) 
                    yok = 0.0;
                end
            end
            if (yok > 1.0) 
                if (abs(yok - 1.0) < VPMGSMALL) 
                    yok = 1.0;
                end 
            end
           % end
            else
                yok = 0.0;
        end
        
            for k=1:dime(3)
 %              zok = 0.0; 
                z = (k-1)*h(3) + zmin;
        if ( (z < (upperCorner(3)-h(3)/2.)) && ...
             (z > (lowerCorner(3)+h(3)/2.))) 
           zok = 1.0;
        elseif ( (abs(z - lowerCorner(3)) < VPMGSMALL) && ...
                  (bflags(VAPBS_LEFT) == 0)) 
             zok = 1.0;
        elseif ((abs(z - lowerCorner(3)) < VPMGSMALL) && ...
                 (bflags(VAPBS_LEFT) == 1)) 
             zok = 0.5;
        elseif ((abs(z - upperCorner(3)) < VPMGSMALL) && ...
                 (bflags(VAPBS_RIGHT) == 0)) 
             zok = 1.0;
        elseif ((abs(z - upperCorner(3)) < VPMGSMALL) && ...
                 (bflags(VAPBS_RIGHT) == 1)) 
             zok = 0.5;
        elseif ((z > (upperCorner(3) + h(3)/2.)) || (z < (lowerCorner(3) - h(3)/2.))) 
                zok = 0.0;
        elseif ((z < (upperCorner(3) + h(3)/2.)) || (z > (lowerCorner(3) - h(3)/2.))) 
            z0 = max(z - h(3)/2., lowerCorner(3));
            z1 = min(z + h(3)/2., upperCorner(3));
            zok = abs(z1-z0)/h(3);

            if (zok < 0.0) 
                if (abs(zok) < VPMGSMALL) 
                    zok = 0.0;
                end
            end
            if (zok > 1.0) 
                if (abs(zok - 1.0) < VPMGSMALL) 
                    zok = 1.0;
                end 
            end
           % end
            else
                zok = 0.0;
        end
        
                if (abs(xok*yok*zok) < VPMGSMALL) 
                   pvec(i,j,k) = 0.0;
                else
                   pvec(i,j,k) = xok*yok*zok;
                end
            end
        end
    end
