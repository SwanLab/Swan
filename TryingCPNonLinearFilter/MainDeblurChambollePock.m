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

h = (x(end)-x(1)) /(N-1); % uniform grid spacing
alpha = (10*h)^2;
% Set Chambolle-Pock algorithm parameters
lambda = 0.25;%/(1*h)^2;%0.2;         % Regularization parameter
tau = 0.0025;           % Step size for the primal variable
sigma = 0.0025;         % Step size for the dual variable
theta = 1;            % Over-relaxation parameter
maxIter = 50000;       % Maximum number of iterations
tol = 1e-6;           % Convergence tolerance

% Run Chambolle-Pock 1D ROF algorithm
[u, iter] = chambolle_pock_1d_ROF(f_noisy, lambda, tau, sigma, theta, maxIter, tol,x,h);

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



function [u, iter] = chambolle_pock_1d_ROF(f, lambda, tau, sigma, theta, maxIter, tol,x,h)

% Initialize variables
n = length(f);                 % Length of the signal
u = f;                         % Primal variable (initial guess is the input signal)
p = zeros(n, 1);             % Dual variable (n-1 because of finite differences)
ubar = u;                      % Extrapolated variable for u

D = derivativeBackward(x);% Forward finite difference operator
%h = 1;

for iter = 1:maxIter
    
    % 1. Update dual variable (p) using gradient of ubar
    grad_ubar = D*ubar;    % Forward difference for 1D gradient
    p = p + sigma * grad_ubar; 
    p = p ./ max(1, abs(p));   % Proximal operator for dual variable (projection onto l-infinity norm ball)
    
    % 2. Update primal variable (u) using divergence of p
    div_p = D'*p; % Divergence is the negative adjoint of the gradient
    u_old = u;
    u = (u - tau * div_p + tau * lambda * f) / (1 + tau * lambda);
    
    % 3. Extrapolation step for ubar
    ubar = u + theta * (u - u_old);
    
    % Check convergence (using relative change in u)
    if norm(u - u_old) / norm(u) < tol
        break;
    end
end

end

function D = derivative(x)
D = derivativeCentral(x);
%D = derivativeBackward(x);
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
D = D*h;
%D = D;
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
D = D*h;
%D = D;
end
