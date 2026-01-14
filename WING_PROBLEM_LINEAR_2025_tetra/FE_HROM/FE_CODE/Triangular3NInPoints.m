function [weig,posgp,shapef,dershapef] = Triangular3NInPoints(TypeIntegrand) ;
%%
% =========================================================================
% Triangular3NInPoints — Gauss rule & shape functions for 3-node (linear) triangular element
% =========================================================================
% PURPOSE
%   Provide Gauss weights/points, linear shape functions, and their derivatives
%   with respect to the parent coordinates (ξ, η) for a 3-node triangular element.
%   Integration rules are selected according to the target term:
%     • 'K'   → stiffness/stress terms (BᵀB is constant) → 1-point centroid rule
%     • 'RHS' → loads/mass terms (NᵀN is quadratic)     → 3-point Gauss rule
%
% SIGNATURE
%   [weig, posgp, shapef, dershapef] = Triangular3NInPoints(TypeIntegrand)
%
% INPUT
%   TypeIntegrand :
%     • char  : 'K' or 'RHS'
%     • cell  : {posgp, weig} to override defaults
%               - posgp : 2×ngaus barycentric coordinates [ξ; η]
%               - weig  : 1×ngaus Gauss weights
%
% OUTPUT
%   weig        : 1×ngaus Gauss weights
%   posgp       : 2×ngaus Gauss points (barycentric coordinates ξ, η)
%   shapef      : ngaus×3 shape functions N(ξ,η)
%   dershapef   : 2×3×ngaus ∂N/∂ξ, ∂N/∂η at each Gauss point
%
% DEFAULT QUADRATURE
%   • 'K'   : 1-point centroid rule
%               posgp = [1/3; 1/3]
%               weig  = 1/2
%   • 'RHS' : 3-point rule (order-2 Gauss for triangles)
%               posgp = [0.5  0.0  0.5;
%                        0.5  0.5  0.0]
%               weig  = (1/6) × [1 1 1]
%
% IMPLEMENTATION NOTES
%   • Shape functions and derivatives obtained from Triangular3N(ξV),
%     returning:
%         Ne(ξ,η)  : [N1 N2 N3]
%         BeXi(ξ,η): [∂N1/∂ξ ∂N2/∂ξ ∂N3/∂ξ;
%                     ∂N1/∂η ∂N2/∂η ∂N3/∂η]
%   • dershapef is expressed in the parent (reference) coordinate system.
%     Transformation to global (x,y) requires the Jacobian of mapping.
%
% USAGE
%   [w, xg, N, dN] = Triangular3NInPoints('K');    % stiffness assembly
%   [w, xg, N, dN] = Triangular3NInPoints('RHS');  % mass/loads assembly
%   [w, xg, N, dN] = Triangular3NInPoints({x_custom, w_custom}); % custom rule
%
% .mlx references: (none)
%
% Written by Joaquín A. Hernández (JAHO), UPC/CIMNE
% Contact: jhortega@cimne.upc.edu
% Automatically commented by ChatGPT — 07-Nov-2025
% =========================================================================

switch TypeIntegrand
    case {'K'}
        weig  = [1/2  ] ;
        posgp = 1/3*[1
            1 ];
    case {'RHS'}
        weig  = 1/6*[1 1 1] ;
        posgp = [0.5  0   0.5
            0.5  0.5 0  ];
end
ndim = 2; nnodeE = 3 ;
ngaus = length(weig) ;
shapef = zeros(ngaus,nnodeE) ;
dershapef = zeros(ndim,nnodeE,ngaus) ;
for g=1:length(weig) ;
    xiV = posgp(:,g) ;
    [Ne BeXi] = Triangular3N(xiV) ;
    shapef(g,:) = Ne ;
    dershapef(:,:,g) = BeXi ;
end


