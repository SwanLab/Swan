function [K,f] = system_SU_p1(tau,a,nu,xnode,source)
% [K,f] = system_SU_p1(tau,a,nu,xnode)
% System obtained by discretizing the weak form associated to
% the convection-diffusion equation
%               a ux - nu uxx = f
% with gthe stabilized SU method using linear interpolation.
% Boundary conditions are not considered.
%
% Input:
%   tau, tau_c: stabilization parameters
%   a, nu:      equation parameters
%   xnode :     nodes coordinates
%



% Gauss points and weights on the reference element [-1,1]
xipg = [-1/sqrt(3) 1/sqrt(3)]'; 
wpg  = [1 1]';

% Shape functions on the reference element
N_mef   =  [(1-xipg)/2 (1+xipg)/2];  
Nxi_mef =  [-1/2 1/2; -1/2 1/2];

% Number of nodes and elements
npoin = size(xnode,2); 
nelem = npoin-1; 

% Number of Gauss points on each element
ngaus = size(wpg,1);

% Allocation of storage for the matrix and vector
K = zeros(npoin,npoin);
f = zeros(npoin,1);

% MATRIX AND VECTOR CALCULATION
% Loop on the elements
for i=1:nelem
    h = xnode(i+1)-xnode(i);
    xm = (xnode(i+1)+xnode(i))/2;
    weigth = wpg*h/2;
    isp = [i i+1]; % Global numer of the element's nodes
    % Loop on the gaussina points
    for ig=1:ngaus
        N = N_mef(ig,:);
        Nx = Nxi_mef(ig,:)*2/h;
        w_ig = weigth(ig);
        x = xm + h/2*xipg(ig); % x-coordinate of the Gauus point
        s = source.compute(x);
        % Assembly
        K(isp,isp) = K(isp,isp) + w_ig*(N'*a*Nx+Nx'*nu*Nx) ...
                                + w_ig*(tau*a*Nx)'*(a*Nx);
        f(isp) = f(isp) + w_ig*(N)'*s;
    end
end
