function [K,f] = system_SU_p2(tau,tau_c,a,nu,xnode)
% [K,f] = system_SU_p2(a,nu,xnode)
% System obtained by discretizing the weak form associated to
% the convection-diffusion equation
%               a ux - nu uxx = f
% with the stabilized SU method usong parabolic interpolation.
% Boundary conditions are not considered.
%
% Input:
%   tau, tau_c: stabilization parameters
%   a, nu:      equation parameters
%   xnode :     nodes coordinates
%


% Gauss points on the reference element [-1,1]
xipg = [-sqrt(15)/5 0 sqrt(15)/5]'; 
wpg = [5/9 8/9 5/9]';

% definición de las funciones de forma y derivadas en el elemento de referencia
N_mef   =  [(xipg-1).*xipg/2  1-xipg.^2  (xipg+1).*xipg/2];  
Nxi_mef =  [ xipg-1/2 -2*xipg xipg+1/2];

% Number of elements and nodes
numnp = size(xnode,2); 
numel = (numnp-1)/2; 

% Number of Gauss points
ngaus = size(wpg,1);

% Allocation of storage for the matrix and the vector
K = zeros(numnp,numnp);
f = zeros(numnp,1);

% Matrix of stabilization parameters
% tau_c for the end nodes and tau for the middle one
Tau = diag([tau_c,tau,tau_c]);


% MATRIX ANVECTOR CALCULATION
% Loop on the elements
for i=1:numel
    h = xnode(2*i+1)-xnode(2*i);
    weigth = wpg*h;
    isp = [2*i-1 2*i 2*i+1]; % Global number of the element's nodes
    % Loop on the Gauss points
    for ig=1:ngaus
        N = N_mef(ig,:);
        Nx = Nxi_mef(ig,:)/h;
        w_ig = weigth(ig);
        x = xnode(2*i) + h*xipg(ig); % x-coordinate of the Gaussina point
        % Assembly
        K(isp,isp) = K(isp,isp) + w_ig*(N'*a*Nx+Nx'*nu*Nx) ...
                                + w_ig*Tau*(a*Nx)'*(a*Nx);
        f(isp) = f(isp) + w_ig*(N')*SourceTerm(x);
   end
 end
