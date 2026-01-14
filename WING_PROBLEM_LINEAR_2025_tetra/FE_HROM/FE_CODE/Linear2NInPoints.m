function [weig,posgp,shapef,dershapef] = Linear2NInPoints(TypeIntegrand) ;
%%
% =========================================================================
% Linear2NInPoints — Gauss rule & shape functions for 2-node 1D element
% =========================================================================
% PURPOSE
%   Provide Gauss weights/points, shape functions, and derivatives for a
%   2-node linear 1D finite element, with quadrature tailored to the target
%   integrand:
%     • 'K'   → stiffness/stress terms (1-point Gauss)
%     • 'RHS' → loads/mass terms     (2-point Gauss)
%   Custom quadrature can be injected by passing a cell {posgp, weig}.
%
% SIGNATURE
%   [weig, posgp, shapef, dershapef] = Linear2NInPoints(TypeIntegrand)
%
% INPUT
%   TypeIntegrand :
%     • char  : 'K' or 'RHS'
%     • cell  : {posgp, weig} to override default quadrature
%               - posgp : 1×ngaus parent coords
%               - weig  : 1×ngaus weights
%
% OUTPUT
%   weig        : 1×ngaus Gauss weights
%   posgp       : 1×ngaus Gauss points in parent coord (ξ)
%   shapef      : ngaus×2   shape functions at each Gauss point
%   dershapef   : 1×2×ngaus derivatives ∂N/∂ξ at each Gauss point
%
% DETAILS
%   Defaults:
%     • TypeIntegrand='K'   → posgp = [0],        weig = [2]
%     • TypeIntegrand='RHS' → posgp = ±1/√3,      weig = [1 1]
%   Shape functions and derivatives are obtained from Linear2N(xi):
%     Ne(ξ)  = [ (1-ξ)/2, (1+ξ)/2 ]
%     BeXi   = [ -1/2, 1/2 ]
%
% PRACTICAL NOTES
%   • dershapef is w.r.t. parent coordinate ξ; mapping to physical space
%     requires element Jacobian assembled elsewhere.
%   • Using the denser RHS quadrature for mass/loads improves consistency
%     with standard 2-point Gauss integration in 1D.
%
% .mlx references: (none)
%
% Written by Joaquín A. Hernández (JAHO), UPC/CIMNE
% Contact: jhortega@cimne.upc.edu
% Automatically commented by ChatGPT — 07-Nov-2025
% =========================================================================



if ~iscell(TypeIntegrand)
    
    switch TypeIntegrand
        case {'K'}
            % 1 integration point
            weig  = [2] ;
            posgp = 0;
        case {'RHS'}
            % Two integration points
            weig  = [1 1] ;
            posgp = 1/sqrt(3)*[-1 1 ] ;
    end
    
else
    weig = TypeIntegrand{2}' ;
    posgp = TypeIntegrand{1} ;
end

%
% if ~iscell(TypeIntegrand)
%   weig  = [1 1 1 1] ;
%         posgp = 1/sqrt(3)*[-1 1 1 -1
%             -1 -1 1 1 ];
% else
%     weig = TypeIntegrand{2}' ;
%     posgp = TypeIntegrand{1} ;
%
% end
%



ndim = 1; nnodeE = 2 ;
ngaus = length(weig) ;
shapef = zeros(ngaus,nnodeE) ;
dershapef = zeros(ndim,nnodeE,ngaus) ;
for g=1:length(weig) ;
    xi = posgp(g) ;
    [Ne BeXi] = Linear2N(xi) ;
    shapef(g,:) = Ne ;
    dershapef(:,:,g) = BeXi ;
end