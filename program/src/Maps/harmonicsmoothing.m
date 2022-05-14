% fillcoCoefMolDielSmooth (smooth function based on Bruccoleri article)
if (strcmp(srfm,'SMOL')==1)
% define work arrays

a1cf=dielx;
a2cf=diely;
a3cf=dielz;


numptsin=3;
% smooth the dielectric values
for i=1:dime(1)
    for j=1:dime(2)
        for k=1:dime(3)
            % get the 8 points that are 1/sqrt(2) grid spacing away
%            numpts=0;
%            frac=0;
            fracin=1.0/a1cf(i,j,k)+1.0/a2cf(i,j,k)+1.0/a3cf(i,j,k);

            %Points for X-shifted array
            
            frac=fracin;
            numpts=numptsin;
            if (j>1)
                frac=frac+1.0/a2cf(i,j-1,k);
                numpts=numpts+1;
            end
            if (k>1) 
                frac=frac+1.0/a3cf(i,j,k-1);
                numpts=numpts+1;
            end
            if (i < (dime(1)-1))
                frac=frac+1.0/a2cf(i+1,j,k)+1.0/a3cf(i+1,j,k);
                numpts=numpts+2;
                if (j>1)
                    frac=frac+1.0/a2cf(i+1,j-1,k);
                    numpts=numpts+1;
                end
                if (k>1)
                    frac=frac+1.0/a3cf(i+1,j,k-1);
                    numpts=numpts+1;
                end
            end
            dielx(i,j,k)=numpts/frac;

            %Points for Y-shifted array
            
            frac=fracin;
            numpts=numptsin;
            
            if (i>1)   
                frac=frac+1.0/a1cf(i-1,j,k);
                numpts=numpts+1;
            end
            if (k>1) 
                frac=frac+1.0/a3cf(i,j,k-1);
                numpts=numpts+1;
            end
            if (j < (dime(2)-1))
                frac=frac+1.0/a1cf(i,j+1,k)+1.0/a3cf(i,j+1,k);
                numpts=numpts+2;
                if (i>1)
                    frac=frac+1.0/a1cf(i-1,j+1,k);
                    numpts=numpts+1;
                end
                if (k>1)
                    frac=frac+1.0/a3cf(i,j+1,k-1);
                    numpts=numpts+1;
                end
            end
            diely(i,j,k)=numpts/frac;
            
             %Points for Z-shifted array
            
            frac=fracin;
            numpts=numptsin;
            
            if (i>1)    
                frac=frac+1.0/a1cf(i-1,j,k);
                numpts=numpts+1;
            end
            if (j>1)    
                frac=frac+1.0/a2cf(i,j-1,k);
                numpts=numpts+1;
            end
            if (k < (dime(3)-1))
                pepeku=pepe+dime12;
                frac=frac+1.0/a1cf(i,j,k+1)+1.0/a2cf(i,j,k+1);
                numpts=numpts+2;
                if (i>1)
                    frac=frac+1.0/a1cf(i-1,j,k+1);
                    numpts=numpts+1;
                end
                if (j>1)
                    frac=frac+1.0/a2cf(i,j-1,k+1);
                    numpts=numpts+1;
                end
            end
            dielz(i,j,k)=numpts/frac;               
        end
    end
end

clear a1cf a2cf a3cf
end