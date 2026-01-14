function  [Fint,Kred_w_LL] = InternalForcesW(OPERFE,PoneST,PK2STRESS,DATA,VAR,tauNONallDER_q,DOFl) ;
%--------------------------------------------------------------------------
% InternalForcesW
% ----------------
% PURPOSE
%   Assemble the **internal force vector** Fint in a manifold HROM using
%   ECM-selected Gauss points with **latent-dependent weights** (MAW-ECM),
%   and return the reduced, ECM-weight-sensitivity block needed for the
%   Newton tangent:
%       • Fint        = ∑_g  (B_g^T P_g) w_g(q)
%       • Kred_w_LL   = (τ′)^T [∂Fint/∂w]_{DOFl} · (∂w/∂q)   (weight term)
%
% WHAT THIS IMPLEMENTATION DOES
%   1) Builds per-ECM-point contributions fINTredNON_ECM (columns = ECM points)
%        fINTredNON_ECM(:,g) =  ∑_{tensor comp i}  B_{g,i}^T P_{g,i}
%      where P_{g,i} is either the i-th stacked component of:
%        • PoneST (1st PK stress)   → large-strain path, or
%        • PK2STRESS (2nd PK stress)→ small-strain path,
%      and B_{g,i} is the corresponding row-block of OPERFE.Bst.
%   2) Accumulates the internal force with current MAW-ECM weights:
%        Fint = fINTredNON_ECM * OPERFE.wSTs
%   3) Projects the per-ECM matrix to the **tangent manifold** using τ′:
%        fINTredNON_ECM ← (τNONallDER_q)^T · fINTredNON_ECM
%   4) Forms the **weight-sensitivity block** on the free latents (DOFl):
%        Kred_w_LL = fINTredNON_ECM(DOFl,:) * OPERFE.dw_q
%      which is used by the Newton driver to augment the reduced Jacobian
%      when latent-dependent weight terms are enabled.
%
% INPUTS
%   OPERFE            : HROM operators with ECM selections
%                       • Bst      – stacked stress operator on ECM points
%                       • wSTs     – ECM weights at current latent state
%                       • dw_q     – ∂w/∂q (columns align with latent DOFs)
%   PoneST            : Stacked 1st PK stress (use for large-strain path).
%   PK2STRESS         : Stacked 2nd PK stress (use when PoneST is empty).
%   DATA              : Includes mesh dimensions:
%                       • MESH.ndim   – spatial dimension (sets nF = ndim^2)
%                       • MESH.nstrain– strain components for PK2 path
%   VAR               : (Not modified here; included for interface symmetry)
%   tauNONallDER_q    : Block Jacobian of extended embedding [τ′ 0; 0 I_r]
%   DOFl              : Indices of free latent coordinates (rows to extract)
%
% OUTPUTS
%   Fint              : Internal force vector in **embedding coordinates**
%                       (i.e., before the Newton driver projects/blocks).
%   Kred_w_LL         : Weight-sensitivity contribution for the reduced
%                       Jacobian on DOFl (consistent with current τ′ and w).
%
% NOTES / CONVENTIONS
%   • Path selection:
%       - If PoneST is nonempty → large-strain (PK1) path is used.
%       - Else PK2STRESS path is used (nF ← DATA.MESH.nstrain).
%   • Dimensions:
%       - mECM = length(OPERFE.wSTs) columns in fINTredNON_ECM,
%       - nREDall = size(OPERFE.Bst,2) rows (modes of the embedding basis).
%   • Consistency:
%       - Kred_w_LL must be assembled with the **same** τ′ and w used in the
%         residual; this function ensures that by reusing tauNONallDER_q and
%         OPERFE.(wSTs,dw_q) passed in.
%   • The geometric (τ″·r) and other tangent terms are not handled here; the
%     Newton driver/Jacobian builder adds them as needed.
%
% RELATED EQUATIONS (informal)
%   Fint(q) = ∑_g f_g(q) w_g(q),
%   with f_g(q) = ∑_i B_{g,i}^T P_{g,i}(q), and
%   ∂Fint/∂q |_w + ∂Fint/∂w · ∂w/∂q  → the latter is what Kred_w_LL captures
%   after projection to DOFl with τ′.
%
% DEPENDENCIES
%   • Requires OPERFE.Bst, OPERFE.wSTs, OPERFE.dw_q to be consistent with the
%     current latent state (typically set by RetrieveWeightsMAWECM_* helpers).
%
% AUTHOR
%   Joaquín A. Hernández Ortega (JAHO) — jhortega@cimne.upc.edu
%
% COMMENTS UPDATED BY: ChatGPT-5
%   Date of update: 07-NOV-2025, Barcelona, Spain
%--------------------------------------------------------------------------

if nargin == 0
    load('tmp2.mat')
end


if ~isempty(PoneST)
    
    % In this implementation (manifold HROM) we need the separate contribution of Gauss
    % point
    %     \begin{equation}
    % \begin{split}
    %   \cblue{\K^*_{w}} & \defeq \Par{\derpar{\FintREDnon}{\qALL{}{}} }_{\fINTredNON = fixed}  =    \sum_{g=1 }^{m} \overbrace{\fINTredNON_g (\qALL{}{})}^{\nREDall \times 1} \overbrace{\derpar{ w_g(\q)}{\qALL{}{}}}^{1 \times \nREDall}
    %   \end{split}
    % \end{equation}
    %
    
    % Let us compute a ma
    
    %Bst = OPERFE.Bst*tauNONallDER_q ;
    nF = DATA.MESH.ndim^2 ;
    
    
    nREDall = size( OPERFE.Bst,2) ;
    mECM =length(OPERFE.wSTs);  % Number of master integration points
    fINTredNON_ECM = zeros(nREDall,mECM);
    for icomp = 1:nF
        icol = icomp:nF:size(PoneST,1) ;
        
        BstLOC = OPERFE.Bst(icol,:) ;
        % The size of the above matrix is
        %  mECM x nmodesUall
        PoneSTloc = PoneST(icol) ;
        % The size of the above matrix mECMx1
        % Thus, we have to compute
        fINTloc = bsxfun(@times,BstLOC,PoneSTloc);
        % fINTloc is now  a mECM x nmodesUall matrix
        
        % Now we just add the contribution
        fINTredNON_ECM = fINTredNON_ECM + fINTloc';
    end
    
    Fint  = fINTredNON_ECM*OPERFE.wSTs ;
    
    
    
    fINTredNON_ECM = tauNONallDER_q'*fINTredNON_ECM ;
    Kred_w_LL = fINTredNON_ECM(DOFl,:)*OPERFE.dw_q ;
    
    %   Here fINTredNON_ECM is a matrix with as many columns as ECM points.
    % We are only interested in the unconstrained entries
    
    
    %
    %     %
    %         %  if  DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES == 0
    %         for icomp = 1:nF
    %             icol = icomp:nF:size(PoneST,1) ;
    %             PoneST(icol,:) = PoneST(icol,:).*OPERFE.wSTs;
    %         end
    %
    %         Fint_OLD = OPERFE.Bst'*PoneST ;
    %     Fint_OLD = tauNONallDER_q'*Fint_OLD ;
    %
    
    
    
    %     else
    %
    %         % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/108_EIFEM_metamat/04_EIFEM_NONL.mlx
    %          nF = DATA.MESH.ndim^2 ;
    %         for icomp = 1:nF
    %             icol = icomp:nF:length(PoneST) ;
    %             PoneST(icol,:) = VAR.PK1STRESS_incre(icol,:).*OPERFE.wSTs;
    %         end
    %         Fint = OPERFE.KstiffLINEAR*VAR.DISP + OPERFE.Bst'*PoneST ;
    %
    %     end
    
    
    
    
else
    
    
    nF = DATA.MESH.nstrain ;
    
    
    nREDall = size( OPERFE.Bst,2) ;
    mECM =length(OPERFE.wSTs);  % Number of master integration points
    fINTredNON_ECM = zeros(nREDall,mECM);
    for icomp = 1:nF
        icol = icomp:nF:size(PK2STRESS,1) ;
        
        BstLOC = OPERFE.Bst(icol,:) ;
        % The size of the above matrix is
        %  mECM x nmodesUall
        PoneSTloc = PK2STRESS(icol) ;
        % The size of the above matrix mECMx1
        % Thus, we have to compute
        fINTloc = bsxfun(@times,BstLOC,PoneSTloc);
        % fINTloc is now  a mECM x nmodesUall matrix
        
        % Now we just add the contribution
        fINTredNON_ECM = fINTredNON_ECM + fINTloc';
    end
    
    Fint  = fINTredNON_ECM*OPERFE.wSTs ;
    
    
    
    fINTredNON_ECM = tauNONallDER_q'*fINTredNON_ECM ;
    Kred_w_LL = fINTredNON_ECM(DOFl,:)*OPERFE.dw_q ;
    
    %     %   if  DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES == 0
    %     nF = DATA.MESH.nstrain ;
    %     for icomp = 1:nF
    %         icol = icomp:nF:length(PK2STRESS) ;
    %         PK2STRESS(icol,:) = PK2STRESS(icol,:).*OPERFE.wSTs;
    %     end
    %     Fint = OPERFE.Bst'*PK2STRESS ;
    %     else
    %         % Special implementation for EIFEM
    %         % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx
    %         nF = DATA.MESH.nstrain ;
    %         for icomp = 1:nF
    %             icol = icomp:nF:length(PK2STRESS) ;
    %             PK2STRESS(icol,:) = VAR.PoneST(icol,:).*OPERFE.wSTs;
    %         end
    %         Fint = OPERFE.KstiffLINEAR*VAR.DISP + OPERFE.Bst'*PK2STRESS ;
    %     end
end

