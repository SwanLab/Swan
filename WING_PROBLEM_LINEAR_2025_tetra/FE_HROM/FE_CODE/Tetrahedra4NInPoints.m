function [weig,posgp,shapef,dershapef] = Tetrahedra4NInPoints(TypeIntegrand) ;
%%
% =========================================================================
% Tetrahedra4NInPoints — Gauss rule & shape functions for 4-node (linear) tetrahedral element
% =========================================================================
% PURPOSE
%   Provide Gauss weights/points, linear shape functions, and their derivatives
%   w.r.t. parent coordinates (ξ, η, ζ) for a 4-node isoparametric tetrahedron.
%   Integration rules are selected by target term:
%     • 'K'   → stiffness/stress terms (BᵀB is constant) → 1-point centroid rule
%     • 'RHS' → loads/mass terms (NᵀN is quadratic)     → 4-point symmetric rule
%
% SIGNATURE
%   [weig, posgp, shapef, dershapef] = Tetrahedra4NInPoints(TypeIntegrand)
%
% INPUT
%   TypeIntegrand :
%     • char : 'K' or 'RHS'
%             (custom quadrature can be added by extending the interface to accept {posgp, weig})
%
% OUTPUT
%   weig        : 1×ngaus Gauss weights
%   posgp       : 3×ngaus Gauss points in parent coordinates [ξ; η; ζ]
%                 (the 4th barycentric coord is λ4 = 1 − ξ − η − ζ)
%   shapef      : ngaus×4 shape functions N(ξ,η,ζ)
%   dershapef   : 3×4×ngaus array with [∂N/∂ξ; ∂N/∂η; ∂N/∂ζ] at each Gauss point
%
% DEFAULT QUADRATURE
%   • 'K'   : 1-point centroid rule
%               posgp = [1/4; 1/4; 1/4]
%               weig  = 1/6          % volume of the reference tetrahedron
%   • 'RHS' : 4-point symmetric rule (exact for quadratics)
%               a = 0.58541020,  b = 0.13819660
%               posgp = [ a  b  b  b;
%                         b  a  b  b;
%                         b  b  a  b ]
%               weig  = (1/24) × [1 1 1 1]
%
% IMPLEMENTATION NOTES
%   • Shape functions/derivatives are obtained from Tetrahedra4N(ξV), returning:
%         Ne(ξ,η,ζ)  : [N1 N2 N3 N4]
%         BeXi(ξ,η,ζ): [∂N/∂ξ; ∂N/∂η; ∂N/∂ζ] (size 3×4)
%   • dershapef is expressed in parent coordinates; mapping to global (x,y,z)
%     requires the Jacobian of the isoparametric transformation at each Gauss point.
%
% USAGE
%   [w, xg, N, dN] = Tetrahedra4NInPoints('K');     % stiffness assembly
%   [w, xg, N, dN] = Tetrahedra4NInPoints('RHS');   % mass/loads assembly
%
% .mlx references: (none)
%
% Written by Joaquín A. Hernández (JAHO), UPC/CIMNE
% Contact: jhortega@cimne.upc.edu
% Automatically commented by ChatGPT — 07-Nov-2025
% =========================================================================


switch TypeIntegrand
    case {'K'}
        weig  = [1/6  ] ;
        posgp = 1/4*[1 1 1]' ;
    case {'RHS'}
        weig  = 1/24*[1 1 1 1] ;
        a = 0.58541020 ; b = 0.13819660 ;
        posgp = [a b b;
            b  a b ;
            b b a ;
            b b b ]' ;
        
end
ndim = 3; nnodeE = 4 ;
ngaus = length(weig) ;
shapef = zeros(ngaus,nnodeE) ;
dershapef = zeros(ndim,nnodeE,ngaus) ;
for g=1:length(weig) ;
    xiV = posgp(:,g) ;
    [Ne BeXi] = Tetrahedra4N(xiV) ;
    shapef(g,:) = Ne ;
    dershapef(:,:,g) = BeXi ;
end


