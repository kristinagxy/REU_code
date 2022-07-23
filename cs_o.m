function [X,info,lambda] = cs_o(A,b,tol,varargin)
% IRcgls_stripped removes several things from IRcgls:
%                 Reorthogonalization, Enrichment, Preconditioning, Restart
%
% options  = IRcgls_stripped('defaults')
% [X,info] = IRcgls_stripped(A,b)
% [X,info] = IRcgls_stripped(A,b,K)
% [X,info] = IRcgls_stripped(A,b,options)
% [X,info] = IRcgls_stripped(A,b,K,options)
%
% This function applies the CG algorithm implicitly to the normal equations
% for the least squares problem.  We obtain a regularized solution by
% terminating the iterations.  Alternatively this function can be used to
% computed a Tikhonov solution.
%
% With 'defaults' as input returns the default options.  Otherwise outputs
% the iterates specified in K, using max(K) as MaxIter, and using all other
% default options.  With options as input: uses the user-specified options
% and all the other default options.
%
% Inputs:
%  A : either (a) a full or sparse matrix
%             (b) a matrix object that performs the matrix*vector operation
%             (c) user-defined function handle
%  b : right-hand side vector
%  s_l : lower bound of non-zero singular values of A
%  s_u : upper bound of non-zero singular values of A
%  tol : tolerance
%  K : (optional) integer vector that specifies which iterates are returned
%      in X; the maximum number of iterations is assumed to be max(K)
%      [ positive integer | vector of positive components ]
%  options : structure with the following fields (optional)
%      x0         - initial guess for the iterations; default = zero vector
%                   [ array | {'none'} ]
%      MaxIter    - maximum allowed number of iterations
%                   [ positive integer | {100} ]
%                   NOTE: K overrules MaxIter if both are assigned
%      x_true     - true solution; allows us to returns error norms with
%                   respect to x_true at each iteration
%                   [ array | {'none'} ]
%      NoiseLevel - norm of noise in rhs divided by norm of rhs 
%                   [ {'none'} | nonnegative scalar]
%      eta        - safety factor for the discrepancy principle
%                   [ {1.01} | scalar greater than (and close to) 1 ]
%      NE_Rtol    - relative tolerance on the normal equation residual norm
%                   [ {1e-12} | positive integer ]
%      RegParam   - regularization parameter lambda, to be employed if
%                   CGLS is used to solve the regularized problem
%                   (A'*A + lambda^2*L'*L)*x = A'*b;
%                   [ {0} | nonnegative scalar ]
%      RegMatrix  - regularization matrix L, used either as a
%                   regularization matrix to solve the regularized
%                   problem (A'*A + lambda^2*L'*L)*x = A'*b, or as a
%                   priorconditioner for the problem (in this case
%                   one should set RegParam = 0)
%                   [ {'Identity'} | 'Laplacian1D' | 'Laplacian2D' |
%                     matrix | function handle ]
%      IterBar    - shows the progress of the iterations
%                   [ {'on'} | 'off' ]
%      NoStop     - specifies whether the iterations should proceed after
%                   a stopping criterion has been satisfied
%                   [ 'on' | {'off'} ]
% Note: the options structure can be created using the function IRset.

% Outputs:
%   X : computed solutions, stored column-wise (at the iterations listed in K)
%   info: structure with the following fields:
%      its      - number of the last computed iteration
%      saved_iterations - iteration numbers of iterates stored in X 
%      StopFlag - a string that describes the stopping condition:
%                   * Reached maximum number of iterations
%                   * Residual tolerance satisfied (discrepancy principle) 
%                   * Normal equation residual tolerance satisfied
%      StopReg  - struct containing information about the solution that
%                 satisfies the stopping criterion, with the fields:
%                   It   : iteration where the stopping criterion is satisfied
%                   X    : the solution satisfying the stopping criterion
%                   Enrm : the corresponding relative error (requires x_true)
%      Rnrm     - relative residual norms at each iteration
%      NE_Rnrm  - normal eqs relative residual norms
%      Xnrm     - solution norms at each iteration
%      Enrm     - relative error norms (requires x_true) at each iteration
%      BestReg  - struct containing information about the solution that
%                 minimizes Enrm (requires x_true), with the fields:
%                   It   : iteration where the minimum is attained
%                   X    : best solution
%                   Enrm : best relative error
%
% See also: IRcgls, IRhybrid_lsqr, IRnnfcgls, IRrrgmres, IRget, IRset

% Emore REU Program, Summer 2022

% This file is part of the IR Tools package and is distributed under the 
% 3-Clause BSD License. A separate license file should be provided as part 
% of the package.
%​
% Set default values for options.
defaultopt = struct('x0','none', 'MaxIter',500, 'x_true','none',...
    'NoiseLevel','none', 'eta',1.01, 'NE_Rtol',1e-12, 'IterBar','on', ...
    'stdCGLS_out','off', 'NoStop','off', ...
    'RegParam',0, 'RegMatrix','Identity');
  
% If input is 'defaults,' return the default options in X.
if nargin==1 && nargout <= 1 && isequal(A,'defaults')
    X = defaultopt;
    return;
end

defaultopt.restart    = 'off';
defaultopt.verbosity  = 'on';
defaultopt.enrichment = 'none';
% Check for acceptable number of optional input arguments.
switch length(varargin)
    case 0
        K = []; options = [];
    case 1
        if isa(varargin{1}, 'double')
            K = varargin{1}; options = [];
        else
            K = []; options = varargin{1};
        end
    case 2
        if isa(varargin{1}, 'double')
            K = varargin{1}; options = varargin{2};
        else
            K = varargin{2}; options = varargin{1};
        end
        if isfield(options, 'MaxIter') && ~isempty(options.MaxIter) && ... 
                (~isempty(K) && options.MaxIter ~= max(K))
            warning('The value of MaxIter is discarded; the maximum value in K is taken as MaxIter')
        end 
    otherwise
        error('Too many input parameters')
end

if isempty(options)
    options = defaultopt;
end

options = IRset(defaultopt, options);

MaxIter    = IRget(options, 'MaxIter',    [], 'fast');
x_true     = IRget(options, 'x_true',     [], 'fast');
NoiseLevel = IRget(options, 'NoiseLevel', [], 'fast');
eta        = IRget(options, 'eta',        [], 'fast');
NE_Rtol    = IRget(options, 'NE_Rtol',    [], 'fast');
IterBar    = IRget(options, 'IterBar',    [], 'fast');
L          = IRget(options, 'RegMatrix',  [], 'fast');
TikParam   = IRget(options, 'RegParam',   [], 'fast');
NoStop     = IRget(options, 'NoStop',     [], 'fast');
verbose    = IRget(options, 'verbosity',  [], 'fast');

verbose = strcmp(verbose, 'on');

if ischar(TikParam), TikParam = 0; end

if isempty(K)
    K = MaxIter;
end

% Sorting the iteration numbers (in case they are shuffled in input).
K = K(:); K = sort(K,'ascend'); K = unique(K);
if ~((isreal(K) && (all(K > 0)) && all(K == floor(K))))
    error('K must be a vector of positive real integers')
end
if K(end) ~= MaxIter
    MaxIter = K(end); 
end

StopIt = MaxIter;

%  We need to find the number of columns in matrix A, but if A is not given 
%  as a matrix, and no initial guess is given, then we can find it by 
%  computing A'*b.  We need this anyway, so doesn't cost additional work.

[~, info_lsqr] = IRhybrid_lsqr(A, b);
lambda = info_lsqr.StopReg.RegP;

s_u = svds(A,1,'largest')*1.1;

s_l = lambda;

d = (s_u^2+s_l^2)/2;
c = (s_u^2-s_l^2)/2;

[m,n]=size(A);

b = [b;zeros(n,1)];

nrmb = norm(b(:));
nrmAtb = norm(A'*b(1:m));

v = zeros(n,1);
r = b;

if isempty(NoiseLevel) || strcmp(NoiseLevel,'none')
    Rtol = 0;
else
    Rtol = eta*NoiseLevel;
end

% See if an initial guess is given. If not, then use 0 as the initial guess.  
x = IRget(options, 'x0', [], 'fast');

if strcmp(x,'none')
    r = b;
    x = zeros(n,1);
else
    try
        Ax = A_times_vec(A, x);
        r = b(:) - Ax;
        d = d - Atransp_times_vec(A, Ax);
    catch
        error('Check the length of x')
    end
end

% Tikhonov regularization?
%
% If TikParam is set to 'off', then any regularization operators are 
% implemented through preconditioning.  In this case, a regularization
% parameter is not specified, and regularization is enforced through 
% termination of the iteration (i.e., we use iterative regularization).
%
% If TikParam is 'on', then the least squares system is augmented
% with the regularization operator, weighted by the given regularization
% parameter, and the resulting over-determined least squares problem is
% solved.  This is the usual implementation of Tikhonov regularization.
if TikParam == 0
    tik = false;
else
    tik = true;
    if strcmpi(L,'identity')
        L = speye(n);
    elseif strcmpi(L, 'Laplacian1D')
        L = LaplacianMatrix1D(n);
    elseif strcmpi(L, 'Laplacian2D')
        L = LaplacianMatrix2D(n);
    %else
    %   Assume a user has supplied a matrix of function handle for L.
    end
    Lx = A_times_vec(L, TikParam*x);
    r = [ r ; -Lx ];
    d = d - Atransp_times_vec(L, TikParam*Lx);
end

if (tik) && Rtol ~= 0
    warning('With these input options IRcs solves the normal equations associated to the Tikhonov regularized problem. The solution to this problem should be computed with high accuracy, so ''NoiseLevel'' should ideally be ''none'' or 0.')
end

% Declare matrices.
X = zeros(n,length(K));
Xnrm    = zeros(MaxIter,1);
Rnrm    = zeros(MaxIter,1);
NE_Rnrm = zeros(MaxIter,1);

if strcmp(x_true,'none')
    errornorms = false;
else
    errornorms = true;
    Enrm = zeros(MaxIter,1);
    nrmtrue = norm(x_true(:));
    BestReg.It = [];
    BestReg.X = [];
    BestReg.Xnrm = [];
    BestReg.Rnrm = [];
    BestReg.NE_Rnrm = [];
    BestReg.Enrm = [];
    BestEnrm = 1e10;
end

NoStop = strcmp(NoStop,'on');

saved_iterations = zeros(1, length(K));


% Prepare for iterations.
% normr2 = norm(d(:))^2;

% Iterate.
noIterBar = strcmp(IterBar,{'off'});
if ~noIterBar
    h_wait = waitbar(0, 'Running iterations, please wait ...');
end
j = 0;

num_iter = ceil((log(tol)-log(2))/log((s_u-s_l)/(s_u+s_l)))+1;
for k=1:min(num_iter,MaxIter)
    
    if ~noIterBar
        waitbar(k/MaxIter, h_wait)
    end
   


    % update alpha and beta
    if k == 1
        beta = 0;
        alpha = 1/d;
    elseif k == 2
        beta = 1/2*(c/d)*(c/d);
        alpha = 1/(d-c*c/2/d);
    else
        beta = (alpha*c/2)^2;
        alpha = 1/(d-alpha*c*c/4);
    end
    

    % update v, x, and r
    % A'*r(1:ncols) + lambda*r(ncols+1:end)
    s = A'*r(1:m)+lambda*r(m+1:end);
    v = beta * v + s;
    % CHECK THIS FOR TIK:
    x = x + alpha*v;
    r = r - alpha * [A*v; lambda*v];

   
    % Compute norms.
    Xnrm(k)    = norm(x(:));
    Rnrm(k)    = norm(r(:))/nrmb;
    NE_Rnrm(k) = norm(s(:))/nrmAtb;
    if errornorms
        Enrm(k) = norm(x_true(:)-x(:))/nrmtrue;
        if Enrm(k)<BestEnrm
            BestReg.It = k;
            BestReg.X = x;
            BestReg.Xnrm = Xnrm(k);
            BestReg.Rnrm = Rnrm(k);
            BestReg.NE_Rnrm = NE_Rnrm(k);
            BestEnrm = Enrm(k);
            BestReg.Enrm = BestEnrm;
        end
    end
    AlreadySaved = 0;
    if any(k==K)
        j = j+1;
        X(:,j) = x;
        saved_iterations(j) = k;
        AlreadySaved = 1; 
    end




    if (Rnrm(k) <= Rtol) && (StopIt == MaxIter)
        % Stop because residual satisfies ||b-A*x||/||b|| <= Rtol.
        if verbose
            disp('Residual tolerance satisfied')
        end
        StopFlag = 'Residual tolerance satisfied';
        if ~AlreadySaved && ~NoStop
            j = j+1;
            X(:,j) = x;
            saved_iterations(j) = k;
            AlreadySaved = 1;
        end
        StopIt = k;
        StopReg.It = k;
        StopReg.X = x;
        StopReg.Xnrm = Xnrm(k);
        StopReg.Rnrm = Rnrm(k);
        StopReg.NE_Rnrm = NE_Rnrm(k);
        if errornorms
            StopReg.Enrm = Enrm(k);
        end
        if ~ NoStop
            Xnrm    = Xnrm(1:k);
            Rnrm    = Rnrm(1:k);
            NE_Rnrm = NE_Rnrm(1:k);
            if errornorms
                Enrm = Enrm(1:k);
            end
            X = X(:,1:j);
            saved_iterations = saved_iterations(1:j);
            break
        end
    end

    if (NE_Rnrm(k) <= NE_Rtol) && (StopIt == MaxIter)
        if verbose
            disp('Normal equations residual tolerance satisfied')
        end
        StopFlag = 'Normal equations residual tolerance satisfied';
        if ~AlreadySaved && ~NoStop
            j = j+1;
            X(:,j) = x;
            saved_iterations(j) = k;
            AlreadySaved = 1;
        end
        StopIt = k;
        StopReg.It = k;
        StopReg.X = x;
        StopReg.Xnrm = Xnrm(k);
        StopReg.Rnrm = Rnrm(k);
        StopReg.NE_Rnrm = NE_Rnrm(k);
        if errornorms
            StopReg.Enrm = Enrm(k);
        end
        if ~ NoStop
            Xnrm    = Xnrm(1:k);
            Rnrm    = Rnrm(1:k);
            NE_Rnrm = NE_Rnrm(1:k);
            if errornorms
                Enrm = Enrm(1:k);
            end
            X = X(:,1:j);
            saved_iterations = saved_iterations(1:j);
            return
        end
    end
end


if k == MaxIter
    if StopIt == MaxIter
        % Stop because max number of iterations reached.
        if verbose
            disp('Reached maximum number of iterations')
        end
        StopFlag = 'Reached maximum number of iterations';
        if ~AlreadySaved
            j = j+1;
            X(:,j) = x;
            saved_iterations(j) = k;
        end
        StopReg.It = k;
        StopReg.X = x;
        StopReg.Xnrm = Xnrm(k);
        StopReg.Rnrm = Rnrm(k);
        StopReg.NE_Rnrm = NE_Rnrm(k);
        if errornorms
            StopReg.Enrm = Enrm(k);
        end
        Xnrm    = Xnrm(1:k);
        Rnrm    = Rnrm(1:k);
        NE_Rnrm = NE_Rnrm(1:k);
        if errornorms
            Enrm = Enrm(1:k);
        end
        X = X(:,1:j);
        saved_iterations = saved_iterations(1:j);
    end 
end

if k == num_iter
    if StopIt == MaxIter
        % Stop as the algorithm has reached the end
        if verbose
            disp('the algorithm has reached the end')
        end
        StopFlag = 'the algorithm has reached the end';
        if ~AlreadySaved
            j = j+1;
            X(:,j) = x;
            saved_iterations(j) = k;
        end
        StopReg.It = k;
        StopReg.X = x;
        StopReg.Xnrm = Xnrm(k);
        StopReg.Rnrm = Rnrm(k);
        StopReg.NE_Rnrm = NE_Rnrm(k);
        if errornorms
            StopReg.Enrm = Enrm(k);
        end
        Xnrm    = Xnrm(1:k);
        Rnrm    = Rnrm(1:k);
        NE_Rnrm = NE_Rnrm(1:k);
        if errornorms
            Enrm = Enrm(1:k);
        end
        X = X(:,1:j);
        saved_iterations = saved_iterations(1:j);
    end 
end

if ~noIterBar, close(h_wait), end
if nargout==3
  info.its = k;
  info.saved_iterations = saved_iterations(1:j);
  info.StopFlag = StopFlag;
  info.StopReg = StopReg;
  info.Rnrm = Rnrm(1:k);
  info.NE_Rnrm = NE_Rnrm(1:k);
  info.Xnrm = Xnrm(1:k);
  if errornorms
    info.Enrm = Enrm(1:k);
    info.BestReg = BestReg;
  end
end

end