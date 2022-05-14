function [x,flag,relres,iter,resvec] = bicgstabmpbec(A,b,tol,maxit,M1,M2,x0,varargin)
%BICGSTAB   BiConjugate Gradients Stabilized Method.
%   X = BICGSTAB(A,B) attempts to solve the system of linear equations
%   A*X=B for X. The N-by-N coefficient matrix A must be square and the
%   right hand side column vector B must have length N.
%
%   X = BICGSTAB(AFUN,B) accepts a function handle AFUN instead of the
%   matrix A. AFUN(X) accepts a vector input X and returns the
%   matrix-vector product A*X. In all of the following syntaxes, you can
%   replace A by AFUN.
%
%   X = BICGSTAB(A,B,TOL) specifies the tolerance of the method. If TOL is
%   [] then BICGSTAB uses the default, 1e-6.
%
%   X = BICGSTAB(A,B,TOL,MAXIT) specifies the maximum number of iterations.
%   If MAXIT is [] then BICGSTAB uses the default, min(N,20).
%
%   X = BICGSTAB(A,B,TOL,MAXIT,M) and X = BICGSTAB(A,B,TOL,MAXIT,M1,M2) use
%   preconditioner M or M=M1*M2 and effectively solve the system
%   A*inv(M)*X = B for X. If M is [] then a preconditioner is not
%   applied. M may be a function handle returning M\X.
%
%   X = BICGSTAB(A,B,TOL,MAXIT,M1,M2,X0) specifies the initial guess.  If
%   X0 is [] then BICGSTAB uses the default, an all zero vector.
%
%   [X,FLAG] = BICGSTAB(A,B,...) also returns a convergence FLAG:
%    0 BICGSTAB converged to the desired tolerance TOL within MAXIT iterations.
%    1 BICGSTAB iterated MAXIT times but did not converge.
%    2 preconditioner M was ill-conditioned.
%    3 BICGSTAB stagnated (two consecutive iterates were the same).
%    4 one of the scalar quantities calculated during BICGSTAB became
%      too small or too large to continue computing.
%
%   [X,FLAG,RELRES] = BICGSTAB(A,B,...) also returns the relative residual
%   NORM(B-A*X)/NORM(B). If FLAG is 0, then RELRES <= TOL.
%
%   [X,FLAG,RELRES,ITER] = BICGSTAB(A,B,...) also returns the iteration
%   number at which X was computed: 0 <= ITER <= MAXIT. ITER may be an
%   integer + 0.5, indicating convergence half way through an iteration.
%
%   [X,FLAG,RELRES,ITER,RESVEC] = BICGSTAB(A,B,...) also returns a vector
%   of the residual norms at each half iteration, including NORM(B-A*X0).
%
%   Example:
%      n = 21; A = gallery('wilk',n);  b = sum(A,2);
%      tol = 1e-12;  maxit = 15; M = diag([10:-1:1 1 1:10]);
%      x = bicgstab(A,b,tol,maxit,M);
%   Or, use this matrix-vector product function
%      %-----------------------------------------------------------------%
%      function y = afun(x,n)
%      y = [0; x(1:n-1)] + [((n-1)/2:-1:0)'; (1:(n-1)/2)'].*x+[x(2:n); 0];
%      %-----------------------------------------------------------------%
%   and this preconditioner backsolve function
%      %------------------------------------------%
%      function y = mfun(r,n)
%      y = r ./ [((n-1)/2:-1:1)'; 1; (1:(n-1)/2)'];
%      %------------------------------------------%
%   as inputs to BICGSTAB:
%      x1 = bicgstab(@(x)afun(x,n),b,tol,maxit,@(x)mfun(x,n));
%
%   Class support for inputs A,B,M1,M2,X0 and the output of AFUN:
%      float: double
%
%   See also BICG, BICGSTABL, CGS, GMRES, LSQR, MINRES, PCG, QMR, SYMMLQ,
%   TFQMR, ILU, FUNCTION_HANDLE.

%   Copyright 1984-2011 The MathWorks, Inc.
%   $Revision: 1.19.4.8 $ $Date: 2011/05/17 02:33:02 $

% Check for an acceptable number of input arguments
if nargin < 2
    error(message('MATLAB:bicgstab:NotEnoughInputs'));
end

% Determine whether A is a matrix or a function.
[atype,afun,afcnstr] = iterchk(A);
if strcmp(atype,'matrix')
    % Check matrix and right hand side vector inputs have appropriate sizes
    [m,n] = size(A);
    if (m ~= n)
        error(message('MATLAB:bicgstab:NonSquareMatrix'));
    end
    if ~isequal(size(b),[m,1])
        error(message('MATLAB:bicgstab:RSHsizeMatchCoeffMatrix', m));
    end
else
    m = size(b,1);
    n = m;
    if ~iscolumn(b)
        error(message('MATLAB:bicgstab:RSHnotColumn'));
    end
end

% Assign default values to unspecified parameters
if nargin < 3 || isempty(tol)
    tol = 1e-6;
end
warned = 0;
if tol < eps
    warning(message('MATLAB:bicgstab:tooSmallTolerance'));
    warned = 1;
    tol = eps;
elseif tol >= 1
    warning(message('MATLAB:bicgstab:tooBigTolerance'));
    warned = 1;
    tol = 1-eps;
end
if nargin < 4 || isempty(maxit)
    maxit = min(n,20);
end

% Check for all zero right hand side vector => all zero solution
n2b = norm(b);                     % Norm of rhs vector, b
if (n2b == 0)                      % if    rhs vector is all zeros
    x = zeros(n,1);                % then  solution is all zeros
    flag = 0;                      % a valid solution has been obtained
    relres = 0;                    % the relative residual is actually 0/0
    iter = 0;                      % no iterations need be performed
    resvec = 0;                    % resvec(1) = norm(b-A*x) = norm(0)
    if (nargout < 2)
        itermsg('bicgstab',tol,maxit,0,flag,iter,NaN);
    end
    return
end

if ((nargin >= 5) && ~isempty(M1))
    existM1 = 1;
    [m1type,m1fun,m1fcnstr] = iterchk(M1);
    if strcmp(m1type,'matrix')
        if ~isequal(size(M1),[m,m])
            error(message('MATLAB:bicgstab:WrongPrecondSize', m));
        end
    end
else
    existM1 = 0;
    m1type = 'matrix';
end

if ((nargin >= 6) && ~isempty(M2))
    existM2 = 1;
    [m2type,m2fun,m2fcnstr] = iterchk(M2);
    if strcmp(m2type,'matrix')
        if ~isequal(size(M2),[m,m])
            error(message('MATLAB:bicgstab:WrongPrecondSize', m));
        end
    end
else
    existM2 = 0;
    m2type = 'matrix';
end

if ((nargin >= 7) && ~isempty(x0))
    if ~isequal(size(x0),[n,1])
        error(message('MATLAB:bicgstab:WrongInitGuessSize', n));
    else
        x = x0;
    end
else
    x = zeros(n,1);
end

if ((nargin > 7) && strcmp(atype,'matrix') && ...
        strcmp(m1type,'matrix') && strcmp(m2type,'matrix'))
    error(message('MATLAB:bicgstab:TooManyInputs'));
end

% Set up for the method
flag = 1;
xmin = x;                          % Iterate which has minimal residual so far
imin = 0;                          % Iteration at which xmin was computed
tolb = tol * n2b;                  % Relative tolerance
r = b - iterapp('mtimes',afun,atype,afcnstr,x,varargin{:});
normr = norm(r);                   % Norm of residual
normr_act = normr;

if (normr <= tolb)                 % Initial guess is a good enough solution
    flag = 0;
    relres = normr / n2b;
    iter = 0;
    resvec = normr;
    if (nargout < 2)
        itermsg('bicgstab',tol,maxit,0,flag,iter,relres);
    end
    return
end

rt = r;                            % Shadow residual
resvec = zeros(2*maxit+1,1);       % Preallocate vector for norm of residuals
resvec(1) = normr;                 % resvec(1) = norm(b-A*x0)
normrmin = normr;                  % Norm of residual from xmin
rho = 1;
omega = 1;
stag = 0;                          % stagnation of the method
alpha = [];                        % overshadow any functions named alpha
moresteps = 0;
maxmsteps = min([floor(n/50),10,n-maxit]);
maxstagsteps = 3;
% loop over maxit iterations (unless convergence or failure)
fprintf('Solving the linear matrix equation')
for ii = 1 : maxit
        if mod (ii,10)==0 
        fprintf('.')
        end
    rho1 = rho;
    rho = rt' * r;
    if (rho == 0.0) || isinf(rho)
        flag = 4;
        resvec = resvec(1:2*ii-1);
        break
    end
    if ii == 1
        p = r;
    else
        beta = (rho/rho1)*(alpha/omega);
        if (beta == 0) || ~isfinite(beta)
            flag = 4;
            break
        end
        p = r + beta * (p - omega * v);
    end
    if existM1
        ph = iterapp('mldivide',m1fun,m1type,m1fcnstr,p,varargin{:});
        if ~all(isfinite(ph))
            flag = 2;
            resvec = resvec(1:2*ii-1);
            break
        end
    else
        ph = p;
    end
    if existM2
        ph = iterapp('mldivide',m2fun,m2type,m2fcnstr,ph,varargin{:});
        if ~all(isfinite(ph))
            flag = 2;
            resvec = resvec(1:2*ii-1);
            break
        end
    end
    v = iterapp('mtimes',afun,atype,afcnstr,ph,varargin{:});
    rtv = rt' * v;
    if (rtv == 0) || isinf(rtv)
        flag = 4;
        resvec = resvec(1:2*ii-1);
        break
    end
    alpha = rho / rtv;
    if isinf(alpha)
        flag = 4;
        resvec = resvec(1:2*ii-1);
        break
    end
    
    if abs(alpha)*norm(ph) < eps*norm(x)
        stag = stag + 1;
    else
        stag = 0;
    end
    
    xhalf = x + alpha * ph;        % form the "half" iterate
    s = r - alpha * v;             % residual associated with xhalf
    normr = norm(s);
    normr_act = normr;
    resvec(2*ii) = normr;
    
    % check for convergence
    if (normr <= tolb || stag >= maxstagsteps || moresteps)
        s = b - iterapp('mtimes',afun,atype,afcnstr,xhalf,varargin{:});
        normr_act = norm(s);
        resvec(2*ii) = normr_act;
        if normr_act <= tolb
            x = xhalf;
            flag = 0;
            iter = ii - 0.5;
            resvec = resvec(1:2*ii);
            break
        else
            if stag >= maxstagsteps && moresteps == 0
                stag = 0;
            end
            moresteps = moresteps + 1;
            if moresteps >= maxmsteps
                if ~warned
                    warning(message('MATLAB:bicgstab:tooSmallTolerance'));
                end
                flag = 3;
                x = xhalf;
                resvec = resvec(1:2*ii);
                break;
            end
        end
    end
    
    if stag >= maxstagsteps
        flag = 3;
        resvec = resvec(1:2*ii);
        break
    end
    
    if normr_act < normrmin        % update minimal norm quantities
        normrmin = normr_act;
        xmin = xhalf;
        imin = ii - 0.5;
    end
    
    if existM1
        sh = iterapp('mldivide',m1fun,m1type,m1fcnstr,s,varargin{:});
        if ~all(isfinite(sh))
            flag = 2;
            resvec = resvec(1:2*ii);
            break
        end
    else
        sh = s;
    end
    if existM2
        sh = iterapp('mldivide',m2fun,m2type,m2fcnstr,sh,varargin{:});
        if ~all(isfinite(sh))
            flag = 2;
            resvec = resvec(1:2*ii);
            break
        end
    end
    t = iterapp('mtimes',afun,atype,afcnstr,sh,varargin{:});
    tt = t' * t;
    if (tt == 0) || isinf(tt)
        flag = 4;
        resvec = resvec(1:2*ii);
        break
    end
    omega = (t' * s) / tt;
    if isinf(omega)
        flag = 4;
        resvec = resvec(1:2*ii);
        break
    end
    
    if abs(omega)*norm(sh) < eps*norm(xhalf)
        stag = stag + 1;
    else
        stag = 0;
    end
    
    x = xhalf + omega * sh;        % x = (x + alpha * ph) + omega * sh
    r = s - omega * t;
    normr = norm(r);
    normr_act = normr;
    resvec(2*ii+1) = normr;
    
    % check for convergence        
    if (normr <= tolb || stag >= maxstagsteps || moresteps)
        r = b - iterapp('mtimes',afun,atype,afcnstr,x,varargin{:});
        normr_act = norm(r);
        resvec(2*ii+1) = normr_act;
        if normr_act <= tolb
            flag = 0;
            iter = ii;
            resvec = resvec(1:2*ii+1);
            break
        else
            if stag >= maxstagsteps && moresteps == 0
                stag = 0;
            end
            moresteps = moresteps + 1;
            if moresteps >= maxmsteps
                if ~warned
                    warning(message('MATLAB:bicgstab:tooSmallTolerance'));
                end
                flag = 3;
                resvec = resvec(1:2*ii+1);
                break;
            end
        end        
    end
    
    if normr_act < normrmin        % update minimal norm quantities
        normrmin = normr_act;
        xmin = x;
        imin = ii;
    end
    
    if stag >= maxstagsteps
        flag = 3;
        resvec = resvec(1:2*ii+1);
        break
    end
    
end                                % for ii = 1 : maxit
fprintf('\nDone!\n')
% returned solution is first with minimal residual
if flag == 0
    relres = normr_act / n2b;
else
    r = b - iterapp('mtimes',afun,atype,afcnstr,xmin,varargin{:});
    if norm(r) <= normr_act
        x = xmin;
        iter = imin;
        relres = norm(r) / n2b;
    else
        iter = ii;
        relres = normr_act / n2b;
    end
end

% only display a message if the output flag is not used
if nargout < 2
    itermsg('bicgstab',tol,maxit,ii,flag,iter,relres);
end
%
function [atype,afun,afcnstr] = iterchk(A)
%ITERCHK  Checks arguments to iterative methods.
%   [ATYPE,AFUN,AFCNSTR] = ITERCHK(A) returns the following:
%   ATYPE is either 'matrix', 'function', 'expression' or 'inline object'.
%   AFUN is the function name or inline object.
%   AFUN is '' if ATYPE is 'matrix'.
%   AFCNSTR is the function name if ATYPE is 'function'.
%   AFCNSTR is the formula of the function if ATYPE is 'expression' or
%   'inline object'.  AFCNSTR is '' if ATYPE is 'matrix'.
%
%   See also BICG, BICGSTAB, CGS, GMRES, LSQR, MINRES, PCG, QMR, SYMMLQ.

%   Copyright 1984-2004 The MathWorks, Inc. 
%   $Revision: 1.8.4.3 $ $Date: 2010/08/23 23:12:57 $


[afun,afunmsg] = fcnchk(A);
if isempty(afunmsg)
   if isa(afun,'inline')      
      if isa(A,'inline')
         atype = 'inline object';
      else
         atype = 'expression';
      end
      afcnstr = formula(afun);
   else % both function_handles @fun and function names 'fun'
      atype = 'function';
      if isa(A,'function_handle')
          afcnstr = func2str(A);
      else
          afcnstr = A;
      end
   end
elseif isa(A,'float')
   afun = A;
   atype = 'matrix';
   afcnstr = '';
else
   error(message('MATLAB:iterchk:InvalidInput'));
end
%
function y = iterapp(op,afun,atype,afcnstr,x,varargin)
%ITERAPP   Apply matrix operator to vector and error gracefully.
%   ITERAPP(OP,AFUN,ATYPE,AFCNSTR,X) applies matrix operator AFUN to vector
%   X. If ATYPE is 'matrix, then AFUN is a matrix and the OP is applied
%   directly. OP is either 'mtimes' or 'mldivide'.
%   ATYPE and AFCNSTR are used in case of error.
%   ITERAPP(OP,AFUN,ATYPE,AFCNSTR,X,P1,P2,...) allows extra arguments to
%   AFUN(X,P1,P2,...) although this usage is now discouraged in favor of
%   using anonymous functions.
%   AFUN(X,P1,P2,...,PN,TFLAG) should accept a TFLAG as its final input if
%   the calling function is BICG, LSQR or QMR. TFLAG is either 'transp' or
%   'notransp' depending on whether A' OP X or A OP X is required.
%   ITERAPP is designed for use by iterative methods like PCG which
%   require matrix operators AFUN representing matrices A to operate on
%   vectors X and return A*X and may also have operators MFUN representing
%   preconditioning matrices M operate on vectors X and return M\X.
%
%   See also BICG, BICGSTAB, CGS, GMRES, LSQR, MINRES, PCG, QMR, SYMMLQ.

%   Copyright 1984-2011 The MathWorks, Inc.
%   $Revision: 1.7.4.8 $ $Date: 2011/05/17 02:33:16 $

if strcmp(atype,'matrix')
    switch lower(op)
        case 'mtimes'
            if (nargin >= 6) && isequal(varargin{end},'transp')
                y = afun' * x;
            else
                y = afun * x;
            end
        case 'mldivide'
            if (nargin >= 6) && isequal(varargin{end},'transp')
                y = afun' \ x;
            else
                y = afun \ x;
            end
        otherwise
            error(message('MATLAB:iterapp:InvalidOp'))
    end
else
    try
        if (nargin >= 6) && isequal(varargin{end},'notransp')
            % A request for A*x coming from BICG, LSQR and QMR
            try
                % New syntax: we now request afun(x,P1,P2,...,PN,'notransp')
                y = afun(x,varargin{:});
            catch
                % Old syntax: we used to accept afun(x,P1,P2,...,PN)
                y = afun(x,varargin{1:end-1});
            end
        else
            % A request for A*x
            % coming from BICGSTAB, CGS, GMRES, MINRES, PCG or SYMMLQ
            % with the call always afun(P1,P2,...,PN)
            % or a request for A'*x coming from
            % BICG, LSQR and QMR in the afun(x,P1,P2,...,PN,'transp') case
            y = afun(x,varargin{:});
        end
    catch ME
        error(message('MATLAB:iterapp:InvalidInput', atype,afcnstr, ME.message));
    end

    if ~iscolumn(y)
        error(message('MATLAB:iterapp:MustReturnColumn', atype, afcnstr));
    end
end