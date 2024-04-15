function [K,f] = system_SGS_p2(tau,tau_c,a,nu,xnode,problem)
% [K,f] = system_SGS_p2(a,nu,xnode)
% System obtained by discretizing the weak form associated 
% to the convection-diffusion equation
%               a ux - nu uxx = f
% with gthe stabilized SGS method using parabolic interpolation.
% Boundary conditions are not considered.
% 
% Input:
%   tau, tau_c: stabilization parameters
%   a, nu:      equation parameters
%   xnode :     nodes coordinates
%



% Gauss points and weigths on the reference element [-1,1]
xipg = [-sqrt(15)/5 0 sqrt(15)/5]'; wpg = [5/9 8/9 5/9]';

% Shape funcitons on the reference element
N_mef   =  [(xipg-1).*xipg/2  1-xipg.^2  (xipg+1).*xipg/2];  
Nxi_mef =  [ xipg-1/2 -2*xipg xipg+1/2];
Nxxi_mef=  [ 1 -2 1; 1 -2 1; 1 -2 1];

% Numer of nodes and elements
numnp = size(xnode,2); 
numel = (numnp-1)/2; 

% Number of Gauss points in each element
ngaus = size(wpg,1);

% Allocation of storage for the matrix and vector
K = zeros(numnp,numnp);
f = zeros(numnp,1);

% Matrix of stabilization parameters
% tau_c for the corner nodes and tau for the middle one
Tau = diag([tau_c,tau,tau_c]);

% MATRIX AND VECTOR CALCULATIONS
% Loop on the elements
for i=1:numel
    h = xnode(2*i+1)-xnode(2*i);
    weigth = wpg*h;
    isp = [2*i-1 2*i 2*i+1]; % Global number of the element's nodes
    % Loop on the Gauss points of the element
    for ig=1:ngaus
        N = N_mef(ig,:);
        Nx = Nxi_mef(ig,:)/h;
        Nxx= Nxxi_mef(ig,:)/(h^2);
        w_ig = weigth(ig);
        x = xnode(2*i) + h*xipg(ig); % x-coordinate of the gaussian point
        % Assembly
        K(isp,isp) = K(isp,isp) + w_ig*(N'*a*Nx+Nx'*nu*Nx) ...
                               + w_ig*Tau*(a*Nx+nu*Nxx)'*(a*Nx-nu*Nxx);
        f(isp) = f(isp) + w_ig*(N'+Tau*(a*Nx+nu*Nxx)')*SourceTerm(x,problem);
    end
 end

