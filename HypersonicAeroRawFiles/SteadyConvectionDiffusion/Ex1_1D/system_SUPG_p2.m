function [K,f] = system_SUPG_p2(tau,tau_c,a,nu,xnode,source)
% [K,f] = system_SUPG_p2(a,nu,xnode)
% System obtained by discretizing the weak form associated to
% the convection-diffusion equation
%               a ux - nu uxx = f
% with the stabilized SUPG method using parabolic interpolation.
% Boundary conditions are not considered
%
% Input:
%   tau, tau_c: stabilization parameters
%   a, nu:      equation parameters
%   xnode :     nodes coordinates
%


% Gauss points on the reference element [-1,1]
xipg = [-sqrt(15)/5 0 sqrt(15)/5]'; 
wpg = [5/9 8/9 5/9]';

% Shape functions and its derivatives on the reference element
N_mef   =  [(xipg-1).*xipg/2  1-xipg.^2  (xipg+1).*xipg/2];  
Nxi_mef =  [ xipg-1/2 -2*xipg xipg+1/2];
Nxxi_mef=  [ 1 -2 1; 1 -2 1; 1 -2 1];

% Number of nodes and elements
numnp = size(xnode,2); 
numel = (numnp-1)/2; 

% Number of Gauss points
ngaus = size(wpg,1);

% Allocation of storage for the matrix and vector
K = zeros(numnp,numnp);
f = zeros(numnp,1);

% Matrix of stabilization parameters
% tau_c for the end nodes and tau for the middle one
Tau = diag([tau_c,tau,tau_c]);

% MATRIX AND VECTOR'S CALCULATION
% Loop on the elements
for i = 1:numel
    h = xnode(2*i+1)-xnode(2*i);
    weigth = wpg*h;
    isp = [2*i-1 2*i 2*i+1]; % Global number of the nodes
    % Loop on the Gauss points
    for ig=1:ngaus
        N = N_mef(ig,:);
        Nx = Nxi_mef(ig,:)/h;
        Nxx= Nxxi_mef(ig,:)/(h^2);
        w_ig = weigth(ig);
        x = xnode(2*i) + h*xipg(ig); % x-coordinate of the Gaussian point
        s = source.compute(x);
        % Assembly
        K(isp,isp) = K(isp,isp) + w_ig*(N'*a*Nx+Nx'*nu*Nx) ...
                                + w_ig*Tau*(a*Nx)'*(a*Nx-nu*Nxx);
        f(isp) = f(isp) + w_ig*(N'+Tau*(a*Nx)')*s;
    end
end

