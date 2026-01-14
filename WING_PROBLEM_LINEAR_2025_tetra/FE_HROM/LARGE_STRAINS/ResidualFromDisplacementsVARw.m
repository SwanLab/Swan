function  [VAR,celastST,FgradST,detFgrad,Kred_w_LL] =...
    ResidualFromDisplacementsVARw(OPERFE,VAR,MATPRO,DATA,VARint_n,tauNONallDER_q,DOFl) 
%--------------------------------------------------------------------------
% ResidualFromDisplacementsVARw
% -----------------------------
% VARiation of ResidualFromDisplacementsVAR — **VAriable ECM weights**
% JAHO, 22-Sept-2025, HGs Pedralbes, BArcelona
%
% PURPOSE
%   Core residual/internal-force assembler for manifold HROM with **variable
%   ECM weights** (MAW-ECM). Compared to ResidualFromDisplacementsVAR, this
%   variant:
%     • accepts the block Jacobian of the extended embedding τ′ (tauNONallDER_q),
%     • computes internal forces using MAW-ECM-aware operator (InternalForcesW),
%     • returns the ECM-weighted reduced tangent block on DOFl (Kred_w_LL).
%
% INPUTS
%   OPERFE      : FE/HROM operators and options
%                 - (optional) KinternalFORCES_given : linear shortcut (matrix
%                   or piecewise-constant struct with MATRICES/INTERVALS)
%                 - (other standard fields used by StressesFromDisplacementsVAR
%                   and InternalForcesW)
%   VAR         : State container at the current iterate
%                 - DISP          : current reduced coordinates (extended state
%                                   is handled upstream in GetResidual_MHROM)
%                 - FEXT          : external forces (full/extended)
%                 - PK1/PK2STRESS : (filled here for nonlinear path)
%   MATPRO      : Material parameters / constitutive model descriptor
%   DATA        : Run controls and flags
%                 - INTERNAL_FORCES_USING_precomputed_Kstiff (0/1)
%                 - (other standard solver flags)
%   VARint_n    : Internal variables at previous time step (for history models)
%   tauNONallDER_q : Block-Jacobian of the extended embedding
%                    [ τ′(qL)   0;  0   I_r ]  (provided by caller)
%   DOFl        : Indices of free latent coordinates
%
% OUTPUTS
%   VAR         : Updated with
%                 - PK1/PK2STRESS   : stress measures at Gauss points
%                 - FINT            : internal force vector (full/extended)
%                 - RESID           : residual vector FINT - FEXT
%   celastST    : (elastic) tangent operator at Gauss points (Voigt)
%   FgradST     : deformation gradient(s) at Gauss points
%   detFgrad    : determinant(s) of F at Gauss points
%   Kred_w_LL   : ECM-weighted reduced tangent block on DOFl (returned by
%                 InternalForcesW for use in reduced Jacobian assembly)
%
% BEHAVIOR
%   • Nonlinear (default): calls StressesFromDisplacementsVAR to compute
%     stresses/constitutive/tangent quantities, then assembles FINT using
%     InternalForcesW(…, tauNONallDER_q, DOFl) so the result is consistent
%     with the current embedding and MAW-ECM weights.
%   • Linear shortcut: if OPERFE.KinternalFORCES_given is provided, bypass
%     constitutive evaluation and set FINT = K * DISP (piecewise option
%     supported via MATRICES/INTERVALS).
%   • Residual: VAR.RESID = VAR.FINT − VAR.FEXT (full/extended space).
%   • Negative-Jacobian safeguard: if stress path returns empty (e.g., detF ≤ 0)
%     and no linear shortcut is enabled, upper layers should abort the solve.
%
% NOTES
%   • This routine does not perform tangent projection to reduced space;
%     projection and DOFl extraction are handled by the Newton driver.
%   • Kred_w_LL is provided for efficiency so the Newton assembly can reuse
%     the ECM-weighted contributions consistent with the residual evaluation.
%   • Compatible with small- or large-strain paths selected in DATA.
%
% DEPENDENCIES
%   StressesFromDisplacementsVAR, InternalForcesW, DefaultField
%
% AUTHOR
%   Joaquín A. Hernández Ortega (JAHO) — jhortega@cimne.upc.edu
%
% COMMENTS UPDATED BY: ChatGPT-5
%   Date of update: 07-NOV-2025, Barcelona, Spain
%--------------------------------------------------------------------------


if nargin == 0
    load('tmp1.mat')
end



% Stresses from displacements
%[GLSTRAINS,PK2STRESS,celastST,FgradST,PoneST,detFgrad] = StressesFromDisplacements(OPERFE,d,MATPRO,DATA) ;
if DATA.INTERNAL_FORCES_USING_precomputed_Kstiff ==0
    [VAR,celastST,FgradST,detFgrad] = StressesFromDisplacementsVAR(OPERFE,VAR,MATPRO,DATA,VARint_n) ;
else
    celastST = [] ; FgradST = [] ; detFgrad = [] ;
end


if isempty(VAR.PK2STRESS) && DATA.INTERNAL_FORCES_USING_precomputed_Kstiff ==0
    VAR.PoneST = [] ; Fint = [] ; VAR.RESID = [] ;
else
    OPERFE = DefaultField(OPERFE,'KinternalFORCES_given',[]) ;
    % 6.1. Internal forces
    if isempty(OPERFE.KinternalFORCES_given)
     %   if ~isfield(OPERFE,'BstW')    
      
        [VAR.FINT,Kred_w_LL] = InternalForcesW(OPERFE,VAR.PK1STRESS,VAR.PK2STRESS,DATA,VAR,tauNONallDER_q,DOFl) ;
    
      %  else 
            % 3-Jan-2024, see /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/106_ECMtailoredTOmodes/01_VIABILITY.mlx
            % This option proved to be unreliable 
       %  VAR.FINT = InternalForces_TAILOREDWEIGHTS(OPERFE,VAR.PK1STRESS,VAR.PK2STRESS,DATA) ;
       % end
    else
        % FOR the dynamic mode decomposition, see
        %/home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/DynamicModeDecomposition/01_ForcedVibration_2D/RunDMD_aux.mlx
        
        % Internal forces are assumed to be linear with displacements (and given)
        if ~isstruct(OPERFE.KinternalFORCES_given)
            VAR.FINT = OPERFE.KinternalFORCES_given*VAR.DISP ;
        else
            iMATRIX = OPERFE.KinternalFORCES_given.INTERVALS(DATA.istep) ;
            Kloc = OPERFE.KinternalFORCES_given.MATRICES{iMATRIX} ;
            VAR.FINT = Kloc*VAR.DISP ;
        end
    end
    
    % 6.2. Residual
    VAR.RESID  = VAR.FINT- VAR.FEXT;
end


 