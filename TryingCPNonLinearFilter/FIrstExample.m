function MainExample
L = 1;
N = 2000;
x = linspace(0,L,N)';


solveLaplacianNonUniform(x);
projectAndSmoothStepInverse(x);
runningChambollePock(x)


end

function  solveLaplacianNonUniform(x)
% Example PDE: -u''=f with f=sin(pi*x), exact u=sin(pi*x)/pi^2
f = sin(pi*x);
u_true = sin(pi*x)/(pi^2);

D_cen = derivative(x);
D2 = D_cen' * D_cen; % approximates -u''

% Impose Dirichlet BC by removing first and last rows/cols
D2_interior = D2(2:end-1,2:end-1);
f_interior = f(2:end-1);

u_interior = D2_interior \ f_interior;

u = zeros(length(x),1);
u(2:end-1) = u_interior;

rel_err = max(abs(u - u_true)) / max(abs(u_true));
fprintf('Relative error PDE solve: %e\n', rel_err);

figure;
plot(x, u, 'b-', 'LineWidth', 1.5); hold on;
plot(x, u_true, 'r--', 'LineWidth', 2);
xlabel('x'); ylabel('u(x)');
legend('Computed','True');
title('Solution of -u''''=sin(pi*x) with Non-uniform Grid');
grid on;

end

function runningChambollePock(x)

% Example usage of Chambolle's projection algorithm for TV denoising

% Create a noisy signal
N = length(x);
f_clean = sin(2*pi*x) + 0.1*sin(10*pi*x);
noise = 0.02*randn(N,1);
f_noisy = f_clean + noise;

% Set parameters for Chambolle's projection
lambda = 0.1;     % Regularization parameter (adjust as needed)
maxIter = 1000;     % Number of iterations

% Apply Chambolle's TV denoising
u_denoised = chambolle_projection(f_noisy, lambda, maxIter);

% Plot results
figure;
plot(x, f_clean, 'g-', 'LineWidth', 2); hold on;
plot(x, f_noisy, 'r-', 'LineWidth', 1);
plot(x, u_denoised, 'b-', 'LineWidth', 2);
legend('Original Clean Signal', 'Noisy Signal', 'Denoised Signal');
xlabel('x'); ylabel('u(x)');
title('Chambolle TV Denoising Example');
grid on;

end

function u = chambolle_projection(f, lambda, maxIter)
    % Chambolle's projection for TV denoising in 1D
    % Solves: min_u (1/2)*||u - f||^2 + lambda * TV(u)
    
    N = length(f);
    p = zeros(N-1,1);  % Dual variable has length N-1
    tau = 0.15;         % Step size (should satisfy tau <= 0.25 for stability)
    
    for k = 1:maxIter
        % Compute divergence of p
        divp = [p(1); p(2:end) - p(1:end-1); -p(end)];
        
        % Update primal variable u
        u = f - lambda * divp;
        
        % Compute gradient of u
        grad_u = u(2:end) - u(1:end-1);
        
        % Update dual variable p without scaling by lambda
        p_new = p + tau * grad_u;
        
        % Project p to ensure |p| <= 1
        p = p_new ./ max(1, abs(p_new));
    end
    
    % Final update to u after iterations
    divp = [p(1); p(2:end) - p(1:end-1); -p(end)];
    u = f - lambda * divp;
end



function [u_original, u_smoothed] = projectAndSmoothStepInverse(x)
N = length(x);
h = (x(end)-x(1)) /(N-1); % uniform grid spacing
alpha = (10*h)^2;
D_cen = derivative(x);
D2 = D_cen' * D_cen; % approximates -u''

% Original step function: u=0 for x<0.5, u=1 for x>=0.5, with BC u(0)=0,u(1)=0
u_original = zeros(N,1);
u_original(x>=0.25 & x<=0.75) = 1;
u_original(1)=0; u_original(end)=0;

% Backward Euler-like step: (I + alpha D2)*u_new = u_old
I = speye(N);

% Impose Dirichlet BC:
I(1,:) = 0; I(1,1)=1;
I(end,:) = 0; I(end,end)=1;

D2(1,:) = 0; D2(1,1)=1;
D2(end,:) = 0; D2(end,end)=1;

A = I + alpha*D2;
u_smoothed = A \ u_original;


figure;
plot(x, u_original, 'r-', 'LineWidth', 1.5); hold on;
plot(x, u_smoothed, 'b-', 'LineWidth', 1.5);
xlabel('x'); ylabel('u(x)');
legend('Original Step','Smoothed Step');
title('Step Function Before and After Inverse-Based Laplacian Smoothing with alpha scaled by h^2');
grid on;
end




function D_cen = derivative(x)
D_cen = derivativeCentral(x);
end


function [D_cen] = derivativeCentral(x)
N = length(x);
dxm = x(2:end-1)-x(1:end-2);
dxp = x(3:end)-x(2:end-1);
a = -dxp./(dxm.*(dxm+dxp));
b = (dxp - dxm)./(dxm.*dxp);
c = dxm./(dxp.*(dxm+dxp));
A = zeros(N,3); A(2:N-1,1)=a; A(2:N-1,2)=b; A(2:N-1,3)=c;
dx = x(2)-x(1); A(1,1)=-1/dx;A(1,2)=0;A(1,3)=1/dx;
dx = x(end)-x(end-1); A(end,1)=-1/dx;A(end,2)=0;A(end,3)=1/dx;
D_cen = spdiags(A,[-1 0 1],N,N);
end

function [D_bwd] = derivativeBackward(x)
N = length(x);
dxm = x(2:end-1)-x(1:end-2);
dxp = x(3:end)-x(2:end-1);
a = -dxp./(dxm.*(dxm+dxp));
b = (dxp - dxm)./(dxm.*dxp);
c = dxm./(dxp.*(dxm+dxp));
A = zeros(N,3); A(2:N-1,1)=a; A(2:N-1,2)=b; A(2:N-1,3)=c;
dx = x(2)-x(1); A(1,1)=-1/dx;A(1,2)=0;A(1,3)=1/dx;
dx = x(end)-x(end-1); A(end,1)=-1/dx;A(end,2)=0;A(end,3)=1/dx;
xb = x(2:end)-x(1:end-1);
B = zeros(N,2);
B(2:N,1)=-1./xb; B(2:N,2)=1./xb;
dx = x(2)-x(1); B(1,1)=-1/dx;B(1,2)=1/dx;
D_bwd = spdiags(B,[-1 0],N,N);
end