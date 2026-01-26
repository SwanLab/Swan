function [weig,posgp,shapef,dershapef] = Linear3NInPoints(TypeIntegrand) ;
%%
% =========================================================================
% Linear3NInPoints — Gauss rule & shape functions for 3-node 1D (quadratic) element
% =========================================================================
% PURPOSE
%   Provide Gauss weights/points, quadratic shape functions, and their ξ-derivatives
%   for a 3-node 1D finite element, using integration rules tailored to the target:
%     • 'K'   → stiffness/stress terms (BᵀB is quadratic) → 2-point Gauss
%     • 'RHS' → loads/mass terms (NᵀN is quartic)         → 3-point Gauss (Gauss–Legendre)
%   Custom quadrature can be injected by passing {posgp, weig}.
%
% SIGNATURE
%   [weig, posgp, shapef, dershapef] = Linear3NInPoints(TypeIntegrand)
%
% INPUT
%   TypeIntegrand :
%     • char  : 'K' or 'RHS'
%     • cell  : {posgp, weig} to override defaults
%               - posgp : 1×ngaus parent coords ξ
%               - weig  : 1×ngaus weights
%
% OUTPUT
%   weig        : 1×ngaus Gauss weights
%   posgp       : 1×ngaus Gauss points (parent coordinate ξ)
%   shapef      : ngaus×3 quadratic shape functions N(ξ)
%   dershapef   : 1×3×ngaus ∂N/∂ξ at each Gauss point
%
% DEFAULT QUADRATURE
%   • 'K'   : posgp = ±1/√3,  weig = [1 1]
%   • 'RHS' : posgp = [−√(3/5), √(3/5), 0],  weig = [5/9, 5/9, 8/9]
%
% IMPLEMENTATION NOTES
%   • Shape functions from Quadratic3N(ξ), returning:
%       Ne(ξ)  : [N1 N2 N3]
%       BeXi   : [∂N1/∂ξ ∂N2/∂ξ ∂N3/∂ξ]
%   • dershapef is w.r.t. ξ; mapping to physical coordinates needs the Jacobian.
%
% USAGE
%   [w, xg, N, dN] = Linear3NInPoints('K');    % stiffness assembly
%   [w, xg, N, dN] = Linear3NInPoints('RHS');  % mass/loads assembly
%   [w, xg, N, dN] = Linear3NInPoints({x_custom, w_custom}); % custom rule
%
% .mlx references: (none)
%
% Written by Joaquín A. Hernández (JAHO), UPC/CIMNE
% Contact: jhortega@cimne.upc.edu
% Automatically commented by ChatGPT — 07-Nov-2025
% =========================================================================



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

if ~iscell(TypeIntegrand) 
    switch TypeIntegrand
        case {'K'}
            % Two integration points
            weig  = [1 1] ;
            posgp = 1/sqrt(3)*[-1 1 ] ;
        case {'RHS'}
            weig = [5/9  5/9 8/9 ];
            p = sqrt(3/5) ;
            posgp = [-p   p  0] ;
            
    end
else
    
    weig = TypeIntegrand{2}' ;
    posgp = TypeIntegrand{1} ;
    
end

ndim = 1; nnodeE = 3 ;
ngaus = length(weig) ;
shapef = zeros(ngaus,nnodeE) ;
dershapef = zeros(ndim,nnodeE,ngaus) ;
for g=1:length(weig) ;
    xi = posgp(g) ;
    [Ne BeXi] = Quadratic3N(xi) ;
    shapef(g,:) = Ne ;
    dershapef(:,:,g) = BeXi ;
end