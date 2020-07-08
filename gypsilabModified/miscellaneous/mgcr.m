function [xk,flag,err,iter] = mgcr(varargin)
% Copyright (c) 20015, Matthieu Aussal, Ecole Polytechnique       
% GNU General Public License v3.0. 
% Generalized minimum residual method for multiple right hand side
% "Soudais, P. (1994). Iterative solution of a 3-D scattering problem from 
% arbitrary shaped multidielectric and multiconducting bodies. 
% IEEE transactions on antennas and propagation, 42(7), 954-959."

% Infos
disp('Start MGCR : Multiple Generalized Conjugate Residual (no restart)')

% Default parameters
A     = varargin{1};
B     = varargin{2};
tol   = 1e-6;
maxit = 10;
Mm1   = @(V) V;
X0    = zeros(size(B)); 

% Input analysis
if (nargin >= 4) && ~isempty(varargin{4})
    tol = varargin{4};
end
if (nargin >= 5) && ~isempty(varargin{5})
    maxit = varargin{5};
end
if (nargin >= 6) && ~isempty(varargin{6})
    M = varargin{6};
    if isnumeric(M)
        Mm1 = @(V) varargin{6} \ V;
    else
        Mm1 = M;
    end    
end
if (nargin >= 7) && ~isempty(varargin{7})
    L = varargin{6};
    U = varargin{7};    
    Mm1 = @(V) U \ (L \ V); 
end
if (nargin == 8) && ~isempty(varargin{8})
    X0 = varargin{8};
end

% Matrix-Vector product is handle function
if isnumeric(A)
    A = @(V) A * V;
end

% Initialisation
xk   = X0;
rk   = Mm1(A(xk) - B);
nrk1 = sqrt(sum(conj(rk).*rk));
p    = -rk(:,1);
Ap   = Mm1(A(p));
n2Ap = Ap'*Ap;
err  = max(sqrt(sum(conj(rk).*rk))./nrk1);
iter = 1;

% Disk save
save('mgcr_tmp.mat','xk','err','iter')

% Iterative resolution
while (err > tol) && (iter <= maxit)
    % Clock
    t = tps;
    
    % Update all right hand side
    alphak = - (Ap(:,iter)'*rk) ./ n2Ap(iter);
    xk     = xk + p(:,iter)*alphak;
    rk     = rk + Ap(:,iter)*alphak;
    
    % Error
    [err,kmax]  = max(sqrt(sum(conj(rk).*rk))./nrk1);
    
    % Weights
    Ar    = Mm1(A(rk(:,kmax)));
    sigma = (Ap'*Ar) ./ n2Ap.';
    
    % Incrementation
    p(:,iter+1)  = - rk(:,kmax) + p*sigma;
    Ap(:,iter+1) = - Ar + Ap*sigma;
    n2Ap(iter+1) = Ap(:,iter+1)'*Ap(:,iter+1);
    iter = iter+1;
    
    % Disk save
    save('mgcr_tmp.mat','xk','err','iter')
    
    % Infos
    disp([' + Iteration ',num2str(iter-1),' in ', num2str(tps()-t,'%6.2f'), ...
        ' seconds with relative residual ',num2str(err,'%6.2e'),'.'])
end

% Disk clean
delete('mgcr_tmp.mat')

% Infos
if (iter <= maxit)
    flag = 0;
    disp(['mgcr converged at iteration ',num2str(iter-1), ...
        ' to a solution with relative residual ',num2str(err,'%6.3e'),'.']);
else
    flag = 1;
    disp(['mgcr stopped at iteration ',num2str(iter-1),...
        ' without converging to the desired tolerance ',num2str(tol)])
    disp('because the maximum number of iterations was reached.');
    disp(['The iterate returned (number ',num2str(iter-1),...
        ') has relative residual ',num2str(err,'%6.3e').'.'])
end
end


function t = tps()
t = clock;
t = t(4)*3600 + t(5)*60 + t(6);
end
