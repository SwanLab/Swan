function  [Ne BeXi] = Quadrilateral4N(xiV) ; 
%%
% =========================================================================
% Quadrilateral4N — Shape functions & derivatives for 4-node (bilinear) quadrilateral element
% =========================================================================
% PURPOSE
%   Compute bilinear shape functions N(ξ,η) and their derivatives with respect
%   to the parent coordinates (ξ, η) for a 4-node isoparametric quadrilateral.
%
% SIGNATURE
%   [Ne, BeXi] = Quadrilateral4N(xiV)
%
% INPUT
%   xiV : 2×1 vector of parent coordinates [ξ; η]
%
% OUTPUT
%   Ne   : 1×4 vector of shape functions [N₁ N₂ N₃ N₄]
%   BeXi : 2×4 matrix of derivatives of shape functions with respect to (ξ, η)
%           BeXi = [∂N₁/∂ξ  ∂N₂/∂ξ  ∂N₃/∂ξ  ∂N₄/∂ξ;
%                    ∂N₁/∂η  ∂N₂/∂η  ∂N₃/∂η  ∂N₄/∂η]
%
% FORMULATION
%   Shape functions (bilinear Lagrange interpolation in parent domain):
%       N₁(ξ,η) = ¼(1−ξ)(1−η)
%       N₂(ξ,η) = ¼(1+ξ)(1−η)
%       N₃(ξ,η) = ¼(1+ξ)(1+η)
%       N₄(ξ,η) = ¼(1−ξ)(1+η)
%
%   Their derivatives:
%       ∂N₁/∂ξ = −¼(1−η),   ∂N₁/∂η = −¼(1−ξ)
%       ∂N₂/∂ξ =  +¼(1−η),  ∂N₂/∂η = −¼(1+ξ)
%       ∂N₃/∂ξ =  +¼(1+η),  ∂N₃/∂η =  +¼(1+ξ)
%       ∂N₄/∂ξ =  −¼(1+η),  ∂N₄/∂η =  +¼(1−ξ)
%
% NOTES
%   • The mapping to global (x,y) coordinates requires the Jacobian matrix.
%   • Node numbering follows the standard counterclockwise order:
%         (ξ,η) = (−1,−1) → N₁
%                  (1,−1) → N₂
%                  (1, 1) → N₃
%                 (−1, 1) → N₄
%
% USAGE
%   [N, dN] = Quadrilateral4N([ξ; η]);
%
% .mlx references: (none)
%
% Written by Joaquín A. Hernández (JAHO), UPC/CIMNE
% Contact: jhortega@cimne.upc.edu
% Automatically commented by ChatGPT — 07-Nov-2025
% =========================================================================

xi = xiV(1) ; eta = xiV(2) ; 
% Matrix of shape functions
Ne =0.25*[(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta) ]; 
% Matrix of the gradient of shape functions 
BeXi = 0.25*[ -(1-eta),  (1-eta),  (1+eta), -(1+eta) ; 
             -(1-xi) , -(1+xi) , (1+xi), (1-xi)   ] ; 
