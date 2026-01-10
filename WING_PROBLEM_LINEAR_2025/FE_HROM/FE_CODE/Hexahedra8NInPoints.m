function [weig,posgp,shapef,dershapef] = Hexahedra8NInPoints(TypeIntegrand) ;
%%
% =========================================================================
% Hexahedra8NInPoints — Gauss rule & shape functions for 8-node (trilinear) hexahedral element
% =========================================================================
% PURPOSE
%   Provide Gauss weights/points, trilinear (Q1) shape functions, and their
%   derivatives w.r.t. parent coordinates (ξ, η, ζ) for an 8-node isoparametric
%   hexahedron. By default, a 2×2×2 Gauss–Legendre rule is used, suitable for
%   stiffness and mass/load integrations.
%
% SIGNATURE
%   [weig, posgp, shapef, dershapef] = Hexahedra8NInPoints(TypeIntegrand)
%
% INPUT
%   TypeIntegrand :
%     • char  : 'K' or 'RHS'  (both use the default 2×2×2 Gauss–Legendre rule)
%     • cell  : {posgp, weig} to override the default quadrature
%               - posgp : 3×ngaus parent coordinates [ξ; η; ζ]
%               - weig  : 1×ngaus Gauss weights
%
% OUTPUT
%   weig        : 1×ngaus Gauss weights
%   posgp       : 3×ngaus Gauss points in the parent domain
%   shapef      : ngaus×8 matrix of trilinear shape functions N(ξ,η,ζ)
%   dershapef   : 3×8×ngaus array with [∂N/∂ξ; ∂N/∂η; ∂N/∂ζ] at each Gauss point
%
% DEFAULT QUADRATURE (2×2×2 Gauss–Legendre)
%   • Points  : ξ,η,ζ ∈ {−1/√3, +1/√3} on a tensor grid (8 total)
%   • Weights : all equal to 1
%     The points are ordered in the code as columns of:
%         posgp = (1/√3) * [−1  1  1 −1  −1  1  1 −1;
%                            −1 −1  1  1  −1 −1  1  1;
%                            −1 −1 −1 −1   1  1  1  1]
%
% IMPLEMENTATION NOTES
%   • Shape functions and derivatives are computed via Hexahedra8N(ξV), returning:
%         Ne(ξ,η,ζ)  : [N1 … N8]
%         BeXi(ξ,η,ζ): [∂N/∂ξ; ∂N/∂η; ∂N/∂ζ]  (size 3×8)
%   • dershapef is in parent coordinates; mapping to global (x,y,z) requires
%     the Jacobian of the isoparametric transformation.
%   • Standard node order (counter-clockwise, bottom to top) is assumed, with
%     corner 1 at (−1,−1,−1) in the parent domain.
%
% USAGE
%   [w, xg, N, dN] = Hexahedra8NInPoints('K');         % stiffness/stress assembly
%   [w, xg, N, dN] = Hexahedra8NInPoints('RHS');       % mass/loads assembly
%   [w, xg, N, dN] = Hexahedra8NInPoints({Xc, Wc});    % custom quadrature
%
% .mlx references: (none)
%
% Written by Joaquín A. Hernández (JAHO), UPC/CIMNE
% Contact: jhortega@cimne.upc.edu
% Automatically commented by ChatGPT — 07-Nov-2025
% =========================================================================

if ~iscell(TypeIntegrand)
    
    weig  = [1 1 1 1 1 1 1 1] ;
    posgp = 1/sqrt(3)*[-1  1 1 -1   -1  1 1 -1
        -1 -1 1  1   -1 -1 1  1
        -1 -1 -1 -1   1 1  1  1];
else
    weig = TypeIntegrand{2}' ;
    posgp = TypeIntegrand{1} ;
    
end
ndim = 3; nnodeE = 8 ;
ngaus = length(weig) ;
shapef = zeros(ngaus,nnodeE) ;
dershapef = zeros(ndim,nnodeE,ngaus) ;
for g=1:length(weig) ;
    xiV = posgp(:,g) ;
    [Ne BeXi] = Hexahedra8N(xiV) ;
    shapef(g,:) = Ne ;
    dershapef(:,:,g) = BeXi ;
end
%end