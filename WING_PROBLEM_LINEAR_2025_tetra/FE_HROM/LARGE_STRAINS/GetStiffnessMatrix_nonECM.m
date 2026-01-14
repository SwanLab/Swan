function K = GetStiffnessMatrix_nonECM(OPERFE,DATA,VAR,FgradST,celastST,wMSTeq,ResREDeqMST,tauNONallDER_q,tauNONder2_q)
%--------------------------------------------------------------------------
% function K = GetStiffnessMatrix_nonECM(OPERFE, DATA, VAR, FgradST, celastST, ...
%                                        wMSTeq, ResREDeqMST, tauNONallDER_q, tauNONder2_q)
%
% PURPOSE:
%   Computes the **reduced stiffness matrix K** for a hyperreduced-order model (HROM)
%   based on a **nonlinear solution manifold** τ(q). This function builds the tangent
%   stiffness matrix projected onto the tangent space of τ(q), incorporating geometric
%   and material contributions. It replaces standard FE stiffness assembly by:
%
%     1. Projecting the full-order B-operator onto the tangent space using ∂τ/∂q,
%     2. Replacing standard Gauss weights with ECM-based weights (`wMSTeq`),
%     3. Optionally including geometric stiffness from manifold curvature (∂²τ/∂q²).
%
% INPUTS:
% --------
%   OPERFE : struct
%       Finite element operators and reduction basis, including:
%       - Bst          : full-order stress projection matrix (modified in-place)
%       - wSTs         : integration weights (replaced by wMSTeq)
%       - DOFl         : indices of reduced (free) DOFs
%
%   DATA : struct
%       Global analysis parameters:
%       - SMALL_STRAIN_KINEMATICS : toggle for large/small strain formulation
%       - MESH.ndim               : spatial dimension
%
%   VAR : struct
%       Current state variables, especially:
%       - PK2STRESS : second Piola–Kirchhoff stress tensor (for geometric terms)
%
%   FgradST : [ngaus × ndim^2]
%       Deformation gradient at integration points (used in stiffness assembly).
%
%   celastST : [ngaus × nVoigt^2]
%       Elastic constitutive tensor in Voigt form, used in tangent stiffness computation.
%
%   wMSTeq : [ngaus × 1]
%       Equivalent integration weights at master ECM points (see Eq. 188 in *MLEARNstruct_1.pdf*).
%
%   ResREDeqMST : [nDOF × 1]
%       Residual at master ECM points used for geometric stiffness contribution (Eq. 191).
%
%   tauNONallDER_q : [nFullDOF × nRedDOF]
%       Jacobian of the nonlinear manifold mapping τ(q).
%
%   tauNONder2_q : cell array or struct
%       Second derivatives of τ(q) with respect to q, used to assemble geometric correction terms.
%
% OUTPUT:
% -------
%   K : [nRedDOF × nRedDOF]
%       Reduced tangent stiffness matrix projected onto the manifold.
%
% FEATURES:
% ---------
%   - Projects the B-operator to the reduced basis via multiplication by ∂τ/∂q.
%   - Uses equivalent integration weights from ECM for accurate quadrature.
%   - Supports geometric stiffness from the curvature of τ(q) via ∂²τ/∂q².
%   - Compatible with both small and large strain models (although only large strain is active).
%
% REFERENCES:
%   - *MLEARNstruct_1.pdf*, Section 9.2:
%     Eq. (184): projected tangent matrix
%     Eq. (188): equivalent integration weights
%     Eq. (191–193): geometric correction from τ″(q)
%
% DEPENDENCIES:
%   - KstiffLargeStrains
%   - KstiffSmallStrains (not yet tested)
%
% REMARKS:
% ---------
%   - If `OPERFE.KinternalFORCES_given` is set, it may override automatic stiffness assembly.
%   - Only one master ECM point is currently supported (see Newton caller).
%   - For future extension: allow accumulation over multiple ECM points and test small strain case.
%
% AUTHOR:
%   Joaquín A. Hernández Ortega, UPC – CIMNE
%   Last updated: 22 July 2025, Pedralbes, Barcelona
%--------------------------------------------------------------------------

if nargin == 0
    load('tmp3.mat')
end



OPERFE.Bst = OPERFE.Bst*tauNONallDER_q ;
OPERFE.wSTs = wMSTeq ;

%Here $\KmatREDeqMSTl$
%is the standard  large strain stiffness matrix (more specifically, the contribution of the master point assuming
%that the integration weight is $\wMSTeq$


if DATA.SMALL_STRAIN_KINEMATICS == 0
    
    K = KstiffLargeStrains(OPERFE,VAR.PK2STRESS,FgradST,DATA.MESH.ndim,celastST,DATA) ;
else
    error('Option not tested so far, 22-July-2025')
    K = KstiffSmallStrains(OPERFE,FgradST,DATA.MESH.ndim,celastST,DATA) ;
end


%--------------------------------------------------------------------------
% 3. Optional: Add second-order geometric correction to the stiffness matrix
%    If τ″(q) is available, use Eq. (185) from the reference to include
%    geometric stiffness arising from the curvature of the manifold.
%    This term is: ∑_i [ τ″_i(q) * r_i(u) ]  (scalar-weighted Hessian contributions)
%  tauNONder2_q = [] ;
% disp('borrar')
% NO GEO CONTRIBUTION SO FAR (21-July-2025)
if ~isempty(tauNONder2_q)
    K_geo = 0 ;
    
    % Loop over generalized DOFs (excluding prescribed DOFs)
    DOFl_exted = length(ResREDeqMST)-length(OPERFE.DOFl) ; 
    for imodesEXT = 1:DOFl_exted
        % Geometric stiffness contribution: second derivative * residual component
        K_geo = K_geo + tauNONder2_q(imodesEXT) * ResREDeqMST(imodesEXT);
    end
    
    % Add geometric stiffness contribution to the reduced tangent matrix
    % Only affecting the reduced DOFs (DOFl)
    K(OPERFE.DOFl, OPERFE.DOFl) = K(OPERFE.DOFl, OPERFE.DOFl) + K_geo ;
end



%end