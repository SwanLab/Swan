function [Kll,K_geo] = ...
    GetJacobianMatrix_MHROM(OPERFE, DATA, VAR, FgradST, celastST,tauNONallDER_q,tauNONder2_q,DOFl,Kred_w_LL,...
    USE_GEOMETRIC_term_K,USE_WEIGHT_TERM_k) 
%--------------------------------------------------------------------------
% GetJacobianMatrix_MHROM
% -----------------------
% PURPOSE
%   Assemble the reduced Newton tangent for a manifold HROM at the current
%   iterate, starting from the **full-space** tangent and then adding
%   optional contributions:
%     1) Projection of the full tangent to the manifold’s tangent space:
%            K_ROM = τ′(q)^T · K_full · τ′(q)
%     2) Geometric (manifold-curvature) correction via τ″(q)·r(u) on DOFl
%     3) MAW-ECM weight-sensitivity term Kred_w_LL on DOFl
%
% WHAT THE ROUTINE DOES
%   • Full tangent in embedding coordinates:
%        Kextend = GetStiffnessMatrix(OPERFE, DATA, VAR, FgradST, celastST)
%     (contains material + geometric contributions in the **extended** space).
%   • Tangent-space projection to reduced coordinates using τ′:
%        K = τ′(q)^T · Kextend · τ′(q)      where τ′ ≡ tauNONallDER_q
%   • Optional geometric term (USE_GEOMETRIC_term_K == 1):
%        K_geo = Σ_i [ τ″_i(q) * r_i(u) ]   (assembled in embedding space)
%        K(DOFl,DOFl) ← K(DOFl,DOFl) + K_geo
%   • Optional weight term (USE_WEIGHT_TERM_k == 1):
%        K(DOFl,DOFl) ← K(DOFl,DOFl) + Kred_w_LL
%     where Kred_w_LL has been prebuilt from InternalForcesW as
%     (τ′)^T[∂Fint/∂w]_{DOFl}·(∂w/∂q), consistent with current τ′ and ECM weights.
%   • Optional dynamic correction (if DATA.ISDYNAMIC == 1) via
%        K = KstifflDynamicPart(...)
%   • Optional hydrostatic follower-load correction (if OPERFE.HYDRO present).
%   • If affine/periodic constraints exist (DOFm not empty), returns the
%     condensed block Kll = Aᵀ K A; otherwise Kll = K(DOFl,DOFl).
%
% INPUTS
%   OPERFE : FE/HROM operators and partitions
%            • DOFr, DOFm, A, G, uBAR (optional affine/periodic maps)
%            • HYDRO (optional)      : data for hydrostatic contribution
%   DATA   : Solver/model controls
%            • ISDYNAMIC (0/1), time-integration parameters (if dynamic)
%   VAR    : State at current iterate (includes RESID_extend used in K_geo)
%   FgradST, celastST : Kinematic/constitutive tensors to build Kextend
%   tauNONallDER_q    : τ′(q) (block Jacobian in the extended space)
%   tauNONder2_q      : τ″(q) (third-order tensor for geometric term)
%   DOFl              : Indices of free latent coordinates
%   Kred_w_LL         : MAW-ECM weight-sensitivity contribution on DOFl
%   USE_GEOMETRIC_term_K : logical flag → include τ″·r term when true
%   USE_WEIGHT_TERM_k    : logical flag → include Kred_w_LL when true
%
% OUTPUTS
%   Kll     : Reduced tangent on DOFl (or condensed Aᵀ K A, if DOFm present)
%   K_geo   : Geometric correction assembled for DOFl (0 if disabled)
%
% NUMERICAL NOTES
%   • Geometric correction loops over “generalized” components of the extended
%     residual VAR.RESID_extend. Ensure τ″ has compatible dimensions:
%       size(τ″) → [nEXT × nRED × nRED], so squeeze(τ″(i,:,:)) yields an
%       nRED×nRED slice to be scaled by the residual component r_i.
%   • Kred_w_LL must be computed with the **same** τ′ and ECM weights used
%     in the residual assembly (typically via InternalForcesW).
%   • Dynamic and hydro paths are pass-through augmentations to K (if enabled).
%
% SEE ALSO
%   GetStiffnessMatrix, InternalForcesW, KstifflDynamicPart
%
% HISTORY
%   • Default options: if last two flags are omitted, both geometric and
%     MAW-ECM weight terms are enabled by default.
%
% AUTHOR
%   Joaquín A. Hernández Ortega (JAHO) — jhortega@cimne.upc.edu
%
% COMMENTS UPDATED BY: ChatGPT-5
%   Date of update: 07-NOV-2025, Barcelona, Spain
%--------------------------------------------------------------------------


if nargin == 9
    USE_GEOMETRIC_term_K = 1; 
    USE_WEIGHT_TERM_k = 1; 
end


    %--------------------------------------------------------------------------
            % 1. Compute full-space tangent stiffness matrix based on current deformation state
            %    This includes material and geometric contributions in the full displacement space.
            Kextend = GetStiffnessMatrix(OPERFE, DATA, VAR, FgradST, celastST) ;
            
            %--------------------------------------------------------------------------
            % 2. Project the full stiffness matrix onto the tangent space of the nonlinear manifold
            %    Following Eq. (184):   K_ROM = τ′(q)^T · K · τ′(q)
            %    This yields the reduced tangent matrix governing the update of q in the ROM.
            K = tauNONallDER_q' * Kextend * tauNONallDER_q ;
            
            
            %--------------------------------------------------------------------------
            % 3. Optional: Add second-order geometric correction to the stiffness matrix
            %    If τ″(q) is available, use Eq. (185) from the reference to include
            %    geometric stiffness arising from the curvature of the manifold.
            %    This term is: ∑_i [ τ″_i(q) * r_i(u) ]  (scalar-weighted Hessian contributions)
            %                tauNONder2_q = [] ;
            %                disp('borrar esta parte, en la que pongo vacio tauNONder2_q')
            %   tauNONder2_q = [] ;
              K_geo =0 ;
            if USE_GEOMETRIC_term_K == 1
              
                
                % Loop over generalized DOFs (excluding prescribed DOFs)
                for imodesEXT = 1:(length(VAR.RESID_extend) - length(OPERFE.DOFr))
                    % Geometric stiffness contribution: second derivative * residual component
                    K_geo = K_geo + squeeze(tauNONder2_q(imodesEXT,:,:) * VAR.RESID_extend(imodesEXT));
                end
                
                % Add geometric stiffness contribution to the reduced tangent matrix
                % Only affecting the reduced DOFs (DOFl)
                K(DOFl, DOFl) = K(DOFl, DOFl) + K_geo ;
                
            end
            
     
            
            % CONTRIBUTION OF WEIGHTS
            %  USE_symmetric_form = 1;
            if USE_WEIGHT_TERM_k == 1
                %     if  USE_symmetric_form == 1
                %         K(DOFl, DOFl) =  K(DOFl, DOFl) + 0.5*(Kred_w_LL+Kred_w_LL') ;
                %   else                
                K(DOFl, DOFl) =  K(DOFl, DOFl) + Kred_w_LL ;
                %  end
            end
            
        
            
            if DATA.ISDYNAMIC == 1
                % Dynamic part of the Jacobian
                K  =  KstifflDynamicPart(VAR,OPERFE,DATA.TIME_INTEGRATION_PARAMETERS,DATA.DeltaT,K) ;
            end
            
            if ~isempty(OPERFE.HYDRO)
                % Contribution of hydrostatic forces
                % \FextPR{}{}   = - \NbST{T}  \wSTft{}
                % diag(kPRESSst)*Lbool
                % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/DEVELOPMENTS/DynamicFloatingBodies.tex
                % ------------------------
                K = K - OPERFE.HYDRO.NbST_w'*(ConvertBlockDiag_general(kPRESSst,DATA.MESH.ndim,OPERFE.HYDRO.irowsNDIM,OPERFE.HYDRO.icolsNDIM)*OPERFE.HYDRO.Lbool) ;
            end
            
            
             if  ~isempty(OPERFE.DOFm)
                Kll = OPERFE.A'*(K*OPERFE.A) ;                 
            else
                Kll = K(DOFl,DOFl) ; 
            end