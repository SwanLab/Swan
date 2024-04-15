function [K,f] = system_p2(a,nu,xnode,source)
% [K,f] = system_p2(a,nu,xnode)
% System obtained by discretizing the weak form associated to
% the convection-diffusion equation
%               a ux - nu uxx = f
% using a Galerkin formulation and parabolic interpolation.
% Boundary conditions are not considered
% 
% Input:
%   a, nu:  equation parameters
%   xnode:  nodes coordinates
%



% Gauss points on the reference element [-1,1]
 xipg = [-sqrt(15)/5 0 sqrt(15)/5]'; wpg = [5/9 8/9 5/9]';

% Shape functions 
N_mef   =  [(xipg-1).*xipg/2  1-xipg.^2  (xipg+1).*xipg/2];  
Nxi_mef =  [ xipg-1/2 -2*xipg xipg+1/2];

% Number of nodes and elements
numnp = size(xnode,2); numel = (numnp-1)/2; 

% Number of Gauss points on each element
 ngaus = size(wpg,1);

% Allocation of storage 
K = zeros(numnp,numnp);
f = zeros(numnp,1);


% MATRIX AND VECTOR CALCULATION
% Loop on elements
for i=1:numel
    h = xnode(2*i+1)-xnode(2*i);
    weigth = wpg*h;
    isp = [2*i-1 2*i 2*i+1]; % Global number of the element's nodes
    % Loop on Gauss points
    for ig=1:ngaus
        N = N_mef(ig,:);
        Nx = Nxi_mef(ig,:)/h;
        w_ig = weigth(ig);
        x = xnode(2*i) + h*xipg(ig); %x-coordinate of the gaussian point
        s = source.compute(x);
        % Assembly
        K(isp,isp) = K(isp,isp) + w_ig*(N'*a*Nx+Nx'*nu*Nx);
        f(isp) = f(isp) + w_ig*(N')*s;
    end
end
