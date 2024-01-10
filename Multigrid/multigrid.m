function [x,res,ii,grids] = multigrid(A,b,pre,post,cycle,smooth,grids,maxit,tol,x0,dat)
% Recursive Multigrid Algorithm, Solves Ax = b for Symmetric Diagonally
% Dominant A Matrix
%
% Author: Ben Beaudry
% 
% INPUTS:
%    A = A matrix (n x n)
%    b = Right hand side vector (n x 1)
%    pre = Number of presmoothing iterations
%    post = Number of postsmoothing iterations
%    cycle = Type of multigrid cycle (1=V-cycle, 2=W-cycle, 3=F-cycle)
%    smooth = Smoother type (1=Jacobi, 2=Gauss-Seidel)
%    grids = Max grids in cycle, used for grids level during recursion
%    maxit = Max iterations of solver
%    tol = Tolerance of solver
% OPTIONAL:
%    x0 = Solution Guess
% NOT USER INPUT:
%    dat = grids data for recursion
%
% OUTPUTS:
%    x = Solution vector
%    res = Residual vector
%    ii = Number of top level iterations
%    grids = Max grids in cycle, used for grids level during recursion

if grids==0

    % solve exactly at coarsest grids
    x = A\b;
    res = 0;
    ii = 0;

else

    % inital guess
    if nargin<10
        x = b*0;
    else
        x = x0;
    end

    % create grids data
    if nargin<11
        dat = CreateGridHierarchy(A,grids,smooth);
    end

    % pre allocate
    res = zeros(maxit,1);

    % number of grids
    n = length(dat.RAI);

    for ii = 1:maxit

        % presmoothing
        x = smoothing(A,b,x,dat,smooth,grids,n,pre);

        % compute residual error and retrict to coarse grids
        rh = restrict(A,b,x,dat,grids,n);

        % inital coarse grids guess
        if ii == 1
            eh = rh*0;
        end

        % coarse grids recursion
        [eh,~,~,grids] = multigrid(dat.RAI{n-grids+1},rh,pre,post,...
            cycle,smooth,grids-1,1,tol,eh,dat);

        if cycle == 2 && grids ~= 1 && grids ~= n % W cycle

            % prolongation / interpolate error and correct solution
            x = prolong(x,dat,grids,n,eh);

            % resmoothing
            x = smoothing(A,b,x,dat,smooth,grids,n,1);

            % compute residual error and restrict to coarse grids
            rh = restrict(A,b,x,dat,grids,n);

            % coarse grids recursion
            [eh,~,~,grids] = multigrid(dat.RAI{n-grids+1},rh,pre,post,...
                cycle,smooth,grids-1,1,tol,eh,dat);

        end

        if cycle == 3 && grids ~= 1 % F cycle

            % prolongation / interpolate error and correct solution
            x = prolong(x,dat,grids,n,eh);

            % resmoothing
            x = smoothing(A,b,x,dat,smooth,grids,n,1);

            % compute residual error and restrict to coarse grids
            rh = restrict(A,b,x,dat,grids,n);

            % coarse grids recursion V cycle
            [eh,~,~,grids] = multigrid(dat.RAI{n-grids+1},rh,pre,post,...
                1,smooth,grids-1,1,tol,eh,dat);

        end

        % prolongation / interpolate error and correct solution
        x = prolong(x,dat,grids,n,eh);

        % postsmoothing
        x = smoothing(A,b,x,dat,smooth,grids,n,post);

        % compute residual
        res(ii) = norm(b-A*x);

        % check for convergance
        if res(ii)<tol
            break;
        end

    end
end
grids = grids+1;
end


% Function calls

function dat = CreateGridHierarchy(A,grids,smooth)

for i = 1:grids
    if i == 1
        if smooth == 1

            % create Jacobi
            dat.d{i} = diag(A);
            dat.dinv{i} = 1./dat.d{i};

        elseif smooth == 2

            % create Gauss-Seidel
            dat.L{i} = tril(A);
            dat.U{i} = triu(A,1);

        end

        [m,~] = size(A);

    else
        if smooth == 1

            % create Jacobi
            dat.d{i} = diag(dat.RAI{i-1});
            dat.dinv{i} = 1./dat.d{i};

        elseif smooth == 2

            % create Gauss-Seidel
            dat.L{i} = tril(dat.RAI{i-1});
            dat.U{i} = triu(dat.RAI{i-1},1);

        end

        [m,~] = size(dat.RAI{i-1});

    end

    % create interpolation matrices
    dat.I{i} = spdiagsSpecial([1 2 1],-2:0,m);
    dat.I{i} = dat.I{i}/2;

    % create restriction matrices
    dat.R{i} = dat.I{i}'/2;

    % coarse grid matrix
    if i == 1
        dat.RAI{i} = dat.R{i}*A*dat.I{i};
    else
        dat.RAI{i} = dat.R{i}*dat.RAI{i-1}*dat.I{i};
    end
end
end

function x = smoothing(A,b,x,dat,smooth,grids,n,iters)
for i = 1:iters
    if smooth == 1
        % Jacobi
        x = x + dat.dinv{n-grids+1}.*(b-A*x);
    elseif smooth == 2
        % Guass-Siedel
        x = dat.L{n-grids+1}\(b-dat.U{n-grids+1}*x);
    end
end
end

function rh = restrict(A,b,x,dat,grids,n)
r = b-A*x; % residual
rh = dat.R{n-grids+1}*r;
end

function x = prolong(x,dat,grids,n,eh)
x = x+dat.I{n-grids+1}*eh;
end

function [I] = spdiagsSpecial(v,d,n)
% modified version of spdiags to cut down on matrix size and improve
% performance

if rem(n,2)
    m = (n-1)+1;
    nn = (n-3)/2+1;
    B = ones((n-3)/2+1,1)*v;
else
    m = (n);
    nn = (n-2)/2;
    B = ones((n-2)/2,1)*v;
end
d = d(:);
p = length(d);

% Compute lengths of diagonals:
len = max(0, min(m, nn-d) - max(1, 1-d) + 1);
len = [0; cumsum(len)];

a = zeros(len(p+1), 3);
for k = 1:p
    i = (max(1,1-d(k)):min(m,nn-d(k)))';
    a((len(k)+1):len(k+1),:) = [i i+d(k) B(i+(m>=nn)*d(k),k)];
end

a(:,1) = a(:,1)+(a(:,2)-1);

I = sparse(a(:,1),a(:,2),a(:,3),m,nn);
end