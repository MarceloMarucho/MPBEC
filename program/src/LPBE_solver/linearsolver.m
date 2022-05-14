
if  (strcmp(selectsolver, 'biconjgrad')==1)

% biconjugated gradient first and then
% gmres solver
    disp('Solving the linearized PB equation using the')
    disp('Biconjugated gradient method stabilized by LU matrices')
    disp(' ')

disp('Performing the LU decomposition....')
 setup.type='nofill';
 setup.milu='row';
% 
 [L U]=ilu(A,setup);

disp('Done!....')
disp(' ')
%disp('Solving the linear matrix equation....')
%disp(' ')
 
     
 [pote,flag,relres,iter]=bicgstabmpbec(A,bb,accuracy, max_iteration,L,U);

 disp(' ')


if flag ~= 0
disp(' ')
disp('Sorry!!. The bicgstab function could not solve the linear equation due to....') 

disp(' ')
flag
iter
relres
 disp('flag = 1, bicgstab iterated maxit times but did not converge.')
 disp('flag = 2, Preconditioner M was ill-conditioned.')
 disp('flag = 3, bicgstab stagnated. Two consecutive iterates were the same.')
 disp('flag = 4, One of the scalar quantities calculated during bicgstab became too small or too large to continue computing.')
 disp(' ')
 disp('Use the help function browser to get more information on bicgstab')
disp(' ')
disp('Changing the linear solver.......')
disp('Solving the linearized PB equation using the')
disp('gmres method stabilized by LU matrices')
[pote,flag2,relres,iter]=gmresmpbec(A,bb,10,accuracy,max_iteration,L,U);
if flag2 ~= 0
disp(' ')
%disp('Sorry!!. The bicgstab function could not solve the linear equation due to....') 
disp('Sorry!!. The gmres function could not solve the linear equation due to....') 
disp(' ')
flag2
iter
relres
disp(' ')

disp('flag2 = 1, gmres iterated maxit times but did not converge.')
disp('flag2 = 2, Preconditioner M was ill-conditioned.')
disp('flag2 = 3, gmres stagnated. Two consecutive iterates were the same.')
disp('')
disp('Use the help function browser to get more information on gmres')

end
disp(' ')
disp('Sorry!. Neither gmres nor bicgstab were able to solve the linear system')
disp('Please check the inputfile.m and pqr files.')
disp('If nothing is wrong with those files, please use a smaller number for the value of the digits of precision')
disp('in the file inputfile.m and try it again')
disp(' ')
disp('If you are still having the same problem, please use a smaller number for the value of the number of grid points')
disp('in the file inputfile.m and try it again')
disp(' ')
disp('Thanks!!!')
fail=1;
return
end
%
fail=0;
errorh=relres;
iteration_number=iter;
% disp('Done!....')
disp(' ')

elseif  (strcmp(selectsolver, 'gmres')==1)
disp('Solving the linearized PB equation using the')
disp('gmres method stabilized by LU matrices')
disp(' ')

disp('Performing the LU decomposition....')
 setup.type='nofill';
 setup.milu='row';
% 
 [L U]=ilu(A,setup);

disp('Done!....')
%disp(' ')
%disp('Solving the linear matrix equation....')

 disp(' ')
%  dispstat(sprintf('Begining the process...'),'keepthis','timestamp');
%  for i = 97:100
%      dispstat(sprintf('Progress %d%%',i),'timestamp');
%      %doing some heavy stuff here
%     [pote,flag2,relres,iter]=gmres(A,bb,20,accuracy,max_iteration,L,U);
%  end
%  dispstat('Finished.','keepprev');
  [pote,flag2,relres,iter]=gmresmpbec(A,bb,20,accuracy,max_iteration,L,U);
 disp(' ')

if flag2 ~= 0
disp(' ')
%disp('Sorry!!. The bicgstab function could not solve the linear equation due to....') 
disp('Sorry!!. The gmres function could not solve the linear equation due to....') 
disp(' ')
flag2
iter
relres
disp(' ')

disp('flag2 = 1, gmres iterated maxit times but did not converge.')
disp('flag2 = 2, Preconditioner M was ill-conditioned.')
disp('flag2 = 3, gmres stagnated. Two consecutive iterates were the same.')
disp('')
disp('Use the help function browser to get more information on gmres')
fail=1
return
end
fail=0;
errorh=relres;
iteration_number=iter;
% disp('Done!....')
disp(' ')
else
disp('Solving the linearized PB equation using the')
    disp('minres method ')
    disp(' ')
%    disp('Solving the linear matrix equation....')
ene  = length(bb);
show   = false;
check  = true;
itnlim = ene*5;
[ pote, istop, itn, rnorm, Arnorm, Anorm, Acond, ynorm ] = minresmpbec( A, bb, [], 0.0, show, check, itnlim, accuracy );
%[pote,flag,relres,iter]=gmres(A,bb,10,accuracy,max_iteration);
relres=accuracy;
iter=itn;

% if flag ~= 0
% disp(' ')
% disp('Sorry!!. The bicgstab function could not solve the linear equation due to....') 
% 
% disp(' ')
% flag
% iter
% relres
%  disp('flag = 1, bicgstab iterated maxit times but did not converge.')
%  disp('flag = 2, Preconditioner M was ill-conditioned.')
%  disp('flag = 3, bicgstab stagnated. Two consecutive iterates were the same.')
%  disp('flag = 4, One of the scalar quantities calculated during bicgstab became too small or too large to continue computing.')
%  disp(' ')
%  disp('Use the help function browser to get more information on bicgstab')
% disp(' ')    
% fail=1
% return
% end
fail=0;
errorh=relres;
iteration_number=iter;
% disp('Done!....')
disp(' ')
end


    