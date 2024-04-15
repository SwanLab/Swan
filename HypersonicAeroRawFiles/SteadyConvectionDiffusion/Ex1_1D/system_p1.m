function [K,f] = system_p1(a,nu,xnode,problem)
% [K,f] = system_p1(a,nu,xnode)
% System obtained by discretizing the weak form associated to
% the convection-diffusion equation
%               a ux - nu uxx = f
% using the Galerkin formulation and linear interpolation.
% Boundary conditions are not considered
% 
% Input:
%   a, nu:      equation parameters
%   xnode :     nodes coordinates
%



% Gauss points on the reference element [-1,1]
xipg = [-1/sqrt(3) 1/sqrt(3)]'; 
wpg = [1 1]';

% Shape functions on the reference element
N_mef   =  [(1-xipg)/2 (1+xipg)/2];  
Nxi_mef =  [-1/2 1/2; -1/2 1/2];

% Number of nodes and elements
numnp = size(xnode,2); numel = numnp-1; 

% Number of Gauss points on each element
ngaus = size(wpg,1);

% Allocation
K = zeros(numnp,numnp);
f = zeros(numnp,1);


% MATRIX AND VECTOR CALCULATION

% Loop on elements
for i=1:numel
    h = xnode(i+1)-xnode(i);
    xm= (xnode(i)+xnode(i+1))/2;
    weigth = wpg*h/2;
    isp = [i i+1]; % Global number of the nodes of the current element
    % Lopo on Gauss points
    for ig=1:ngaus
        N = N_mef(ig,:);
        Nx = Nxi_mef(ig,:)*2/h;
        w_ig = weigth(ig);
        x = xm + h/2*xipg(ig); % x-coordinate of the gauss point
        % Assmebly
        K(isp,isp) = K(isp,isp) + w_ig*(N'*(a*Nx)+Nx'*(nu*Nx));
        f(isp) = f(isp) + w_ig*(N')*SourceTerm(x,problem);
    end
end

