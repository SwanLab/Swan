function [weig,posgp,shapef,dershapef] = Quadrilateral4NInPoints(TypeIntegrand) ;
%%
% =========================================================================
% Quadrilateral4NInPoints — Gauss rule & shape functions for 4-node (bilinear) quadrilateral element
% =========================================================================
% PURPOSE
%   Provide Gauss weights/points, bilinear shape functions, and their derivatives
%   with respect to parent coordinates (ξ, η) for a 4-node quadrilateral element.
%   By default, a 2×2 Gauss–Legendre rule is used, suitable for both stiffness and
%   mass/load integrations.
%
% SIGNATURE
%   [weig, posgp, shapef, dershapef] = Quadrilateral4NInPoints(TypeIntegrand)
%
% INPUT
%   TypeIntegrand :
%     • char  : 'K' or 'RHS' (both use 2×2 Gauss–Legendre rule)
%     • cell  : {posgp, weig} to override default integration scheme
%               - posgp : 2×ngaus parent coordinates [ξ; η]
%               - weig  : 1×ngaus Gauss weights
%
% OUTPUT
%   weig        : 1×ngaus Gauss weights
%   posgp       : 2×ngaus Gauss points (ξ, η in parent coordinates)
%   shapef      : ngaus×4 shape functions N(ξ,η)
%   dershapef   : 2×4×ngaus ∂N/∂ξ, ∂N/∂η at each Gauss point
%
% DEFAULT QUADRATURE
%   • 2×2 Gauss–Legendre rule (order 2 exactness for bilinear forms)
%       posgp = (1/√3) × [−1  1  1 −1;
%                         −1 −1  1  1]
%       weig  = [1 1 1 1]
%
% IMPLEMENTATION NOTES
%   • Shape functions and derivatives are computed using Quadrilateral4N(ξV),
%     returning:
%         Ne(ξ,η)  : [N1 N2 N3 N4]
%         BeXi(ξ,η): [∂N1/∂ξ ∂N2/∂ξ ∂N3/∂ξ ∂N4/∂ξ;
%                     ∂N1/∂η ∂N2/∂η ∂N3/∂η ∂N4/∂η]
%   • dershapef is expressed in parent coordinates; mapping to physical (x,y)
%     coordinates requires the Jacobian matrix of the isoparametric mapping.
%
% USAGE
%   [w, xg, N, dN] = Quadrilateral4NInPoints('K');    % stiffness/stress assembly
%   [w, xg, N, dN] = Quadrilateral4NInPoints('RHS');  % mass/loads assembly
%   [w, xg, N, dN] = Quadrilateral4NInPoints({x_custom, w_custom}); % custom rule
%
% .mlx references: (none)
%
% Written by Joaquín A. Hernández (JAHO), UPC/CIMNE
% Contact: jhortega@cimne.upc.edu
% Automatically commented by ChatGPT — 07-Nov-2025
% =========================================================================


%switch TypeIntegrand
%    case {'K','RHS'}
% Four integration points

if ~iscell(TypeIntegrand)
    weig  = [1 1 1 1] ;
    posgp = 1/sqrt(3)*[-1 1 1 -1
        -1 -1 1 1 ];
else
    weig = TypeIntegrand{2}' ;
    posgp = TypeIntegrand{1} ;
    
end


ndim = 2; nnodeE = 4 ;
ngaus = length(weig) ;
shapef = zeros(ngaus,nnodeE) ;
dershapef = zeros(ndim,nnodeE,ngaus) ;
for g=1:length(weig) ;
    xiV = posgp(:,g) ;
    [Ne BeXi] = Quadrilateral4N(xiV) ;
    shapef(g,:) = Ne ;
    dershapef(:,:,g) = BeXi ;
end
%end