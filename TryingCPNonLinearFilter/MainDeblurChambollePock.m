% Generate a synthetic 1D piecewise constant signal
n = 1500; % Length of the signal
f = zeros(n, 1);

segLength = round(n / 4); % Divide signal into 4 equal parts
x = linspace(0,1,n);
f(1:segLength) = 0; %First semgent
f(segLength+1:2*segLength) = 1;        % Second segment
f(2*segLength+1:3*segLength) = -1;     % Third segment
f(3*segLength+1:end) = 0.5;                 % Fourth segment

% Add Gaussian noise to the signal
noise_std = 0.2;      % Standard deviation of the noise
f_noisy = f + noise_std * randn(size(f));

% Plot the original clean signal and the noisy signal
figure;
plot(x,f, 'k-', 'LineWidth', 2, 'DisplayName', 'Original Signal');
hold on;
plot(x,f_noisy, 'r--', 'DisplayName', 'Noisy Signal');
title('Original Signal and Noisy Signal');
legend show;

D = derivative(x);% Forward finite difference operator
L = normest(D);

h = (x(end)-x(1)) /(N-1); % uniform grid spacing
alpha = (10*h)^2;
% Set Chambolle-Pock algorithm parameters
e = 4*h;
lambda = 1/e^2;%1000000000;%%%%%% To think in terms of f    /(1*h)^2;%0.2;         % Regularization parameter
tau = 0.5*1/L;           % Step size for the primal variable
sigma = 0.5*1/L;         % Step size for the dual variable
theta = 1;            % Over-relaxation parameter
maxIter = 150000;       % Maximum number of iterations
tol = 1e-15;           % Convergence tolerance

% Run Chambolle-Pock 1D ROF algorithm

alpha = 1;
beta  = 0.1;
k = 1;
%proxf = @(p) proxLsquared(p,sigma, 0);
%proxf = @(p) proxIndicatorLinf(p,1);
%proxf = @(p) proxIndicatorL2(p, 1);
proxf = @(p) proxDualSquareSuportSegment(p,sigma,k,alpha,beta);
proxg = @(u) proxLsquared(u,tau*lambda, f); %%%General prox


u0 = f_noisy;

[u, iter] = chambollePock(u0, tau, sigma, theta,maxIter,tol,D,proxf,proxg);

% Plot the denoised signal
figure;
plot(x,f, 'k-', 'LineWidth', 2, 'DisplayName', 'Original Signal');
hold on;
plot(x,f_noisy, 'r--', 'DisplayName', 'Noisy Signal');
plot(x,u, 'b', 'LineWidth', 2, 'DisplayName', 'Denoised Signal (ROF)');
title('Original, Noisy, and Denoised Signals');
legend show;

% Display the number of iterations
fprintf('Chambolle-Pock algorithm converged in %d iterations.\n', iter);



function [u, iter] = chambollePock(u0, tau, sigma, theta, maxIter, tol,D,proxSigmaF,proxTauG)
u = u0;
p = zeros(size(u0));
uNew = u;
for iter = 1:maxIter
    p = proxSigmaF(p + sigma*D*uNew);
    uOld = u;
    u = proxTauG(u - tau * D'*p);
    uNew = u + theta * (u - uOld);
    if norm(u - uOld) / norm(u) < tol
        break;
    end
end
end

function p = proxDualSquareSuportSegment(p,sigma,k,alpha,beta)
pdotk = p*k;
s = max(0,pdotk)/(1+sigma/alpha^2) + min(0,pdotk)/(1+sigma/beta^2);
p = s*k;
end

function p = proxIndicatorLinf(p, sigma)
    p = (sigma*p)./ max(sigma,abs(p));
end

function p = proxIndicatorL2(p, sigma)
    p = p * min(1, sigma / norm(p, 2));
end

function u = proxLsquared(u, tau, f)
    % Proximal operator for the l2-norm squared term
    u = (u + tau * f) / (1 + tau);
end

function D = derivative(x)
%D = derivativeCentral(x);
D = derivativeBackward(x);
end

function divP = divergence(p,D)
divP = D'*p;
end


function [D] = derivativeBackward(x)
N = length(x);
h = (x(end)-x(1)) /(N-1); 
xb = x(2:end)-x(1:end-1);
B = zeros(N,2);
B(2:N,1)=-1./xb; B(2:N,2)=1./xb;
dx = x(2)-x(1); B(1,1)=-1/dx;B(1,2)=1/dx;
D = spdiags(B,[-1 0],N,N);
%D = D*h;
D = D;
end

function [D] = derivativeCentral(x)
N = length(x);
h = (x(end)-x(1)) /(N-1); 
dxm = x(2:end-1)-x(1:end-2);
dxp = x(3:end)-x(2:end-1);
a = -dxp./(dxm.*(dxm+dxp));
b = (dxp - dxm)./(dxm.*dxp);
c = dxm./(dxp.*(dxm+dxp));
A = zeros(N,3); A(2:N-1,1)=a; A(2:N-1,2)=b; A(2:N-1,3)=c;
dx = x(2)-x(1); A(1,1)=-1/dx;A(1,2)=0;A(1,3)=1/dx;
dx = x(end)-x(end-1); A(end,1)=-1/dx;A(end,2)=0;A(end,3)=1/dx;
D = spdiags(A,[-1 0 1],N,N);
%D = D*h;
D = D;
end
