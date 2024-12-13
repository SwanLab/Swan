function MainDeblurChambollePock
close all
N = 15000;
xmax = 0;xmin=1;
normD = xmax-xmin;
x = linspace(xmax,xmin,N);
[f,f_noisy] = createNoisyStepFunction(N);

e = normD/20;%4000*h;
lambda = 1/e^2;
D = derivative(x);
L = normest(D);
%h = (x(end)-x(1))/(N-1);
%tau = 1/L;
%sigma = 1/L;
tau(:,1)   = 1./sum(abs(D),1);
sigma(:,1) = 1./sum(abs(D),2);

gamma = 100;
maxIter = 150000;
tol = 1e-6;

%proxf = @(p) proxIndicatorLinf(p,1);
%proxf = @(p) proxIndicatorL2(p, 1);
%proxf = @(p,sigma) proxLsquared(p,sigma, 0);
proxf = @(p,sigma) proxDualSquareSuportSegment(p,sigma);
proxg = @(u,tau) proxLsquared(u,tau*lambda, f);

u0 = f_noisy;
fig_handle = figure;
set(fig_handle, 'WindowState', 'maximized');
plotUF = @(u) plotU(x,u,f,f_noisy,fig_handle);
[u, iter] = chambollePock(u0, tau, sigma, gamma,maxIter,tol,D,proxf,proxg,plotUF);
fprintf('Chambolle-Pock algorithm converged in %d iterations.\n', iter);
plotUF(u)
end

function [f,f_noisy] = createNoisyStepFunction(N)
f = zeros(N, 1);
segLength = round(N / 4);
f(1:segLength) = 0;
f(segLength+1:2*segLength) = 1;
f(2*segLength+1:3*segLength) = -1;
f(3*segLength+1:end) = 0.5;
noise_std = 0.2;
f_noisy = f + noise_std * randn(size(f));
end

function plotU(x, u, f, f_noisy, fig_handle)
if ~isvalid(fig_handle)
fig_handle = figure;
end
figure(fig_handle);
cla;
plot(x, f, 'k-', 'LineWidth', 2, 'DisplayName', 'Original Signal');
hold on;
plot(x, f_noisy, 'r--', 'DisplayName', 'Noisy Signal');
plot(x, u, 'b', 'LineWidth', 2, 'DisplayName', 'Denoised Signal (ROF)');
title('Original, Noisy, and Denoised Signals');
legend show;
drawnow;
end

function [u, iter] = chambollePock(u0, tau, sigma, gamma, maxIter, tol,D,proxSigmaF,proxTauG,plotU)
u = u0;
p = zeros(size(u0));
uNew = u;
for iter = 1:maxIter
    p = proxSigmaF(p + sigma.*(D*uNew),sigma);
    uOld = u;
    u = proxTauG(u - tau.*(D'*p),tau);
    theta = 1./sqrt(1 + 2*gamma*tau);
    %theta = gamma;
    tau = theta.*tau;
    sigma = sigma./theta;    
    uNew = u + theta.*(u - uOld);
    if norm(u - uOld) / norm(u) < tol
        break;
    end
    if  mod(iter,1000) == 0
        plotU(u)
    end
end
end

function p = proxDualSquareSuportSegment(p,sigma)
alpha = 1;
beta  = 0.1;
k = 1;
pdotk = p*k;
s = max(0,pdotk)./(1+sigma/alpha^2) + min(0,pdotk)./(1+sigma/beta^2);
p = s*k;
end

function p = proxIndicatorLinf(p, sigma)
p = (sigma*p)./ max(sigma,abs(p));
end

function p = proxIndicatorL2(p, sigma)
p = p * min(1, sigma./norm(p, 2));
end

function u = proxLsquared(u, tau, f)
u = (u + tau.*f)./(1 + tau);
end

function D = derivative(x)
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
D = D;
end
