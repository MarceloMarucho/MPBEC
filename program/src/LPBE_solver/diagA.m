function [A,b] = diagA( n );
  
%        [A,b] = diagA( n );
% makes a diagonal matrix A = diag(1:n) and rhs vector b.
% Use test1( A,b,shift ) with integer shift to make (A - shift*I) singular.

% 17 Oct 2003: First version for testing minres on singular A.

  
  d = (1:n)';
  A = spdiags(d,0,n,n);
  b = ones(n,1);
