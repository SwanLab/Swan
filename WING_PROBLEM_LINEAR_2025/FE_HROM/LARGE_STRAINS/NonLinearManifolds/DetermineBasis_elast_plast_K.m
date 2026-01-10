function [PhiMaster_lin,PhiMaster_nonl,PhiSlave_nonl,qMASTER_nonl,...
    qSLAVE_nonl,OTHER_output,...
    DATA_evaluateTAU_and_DER, nREDcoor,Kll,DATA_interp] = DetermineBasis_elast_plast_K(SNAPdisp,DOFl,ind_elastic,ind_plastic,...
    DATAoffline,DATA_interp,OTHER_output,OPERFE,INTV_GLO)
% =========================================================================
% DETERMINEBASIS_ELAST_PLAST_K — Elastoplastic Basis with Kll Inner Product
% =========================================================================
% PURPOSE
%   Build an elastoplastic reduced basis from displacement snapshots using a
%   stiffness-weighted inner product (Kll-norm). The routine:
%     1) Extracts the ELASTIC basis from elastic snapshots via weighted SVD
%        (WSVDT) under the Kll metric.
%     2) Projects PLASTIC snapshots onto the Kll-orthogonal complement of the
%        elastic subspace.
%     3) Extracts a PLASTIC basis via randomized (weighted) SVD (SRSVD),
%        truncated by DATAoffline.nmodes_PLASTIC (and the chosen norm).
%     4) Splits PLASTIC basis into:
%           – PhiMaster_nonl : first plastic SVD mode = MASTER plastic variable,
%           – PhiSlave_nonl  : remaining SLAVE plastic modes.
%     5) Computes plastic coordinates (qMASTER_nonl, qSLAVE_nonl) and builds
%        spline mappings qSLAVE_nonl = f(qMASTER_nonl) with τ, τ′, τ″.
%
% CONTEXT / RATIONALE
%   • Kll = K(DOFl,DOFl) provides a strain-energy inner product that enforces
%     physically meaningful orthogonality for the elastic/plastic split.
%   • The MASTER plastic latent variable MUST be the first plastic SVD mode
%     to guarantee an invertible encoder and consistent stress recovery /
%     hyperreduction downstream.
%
% INPUTS
%   SNAPdisp     : [nDOF × nsnap] displacement snapshots (decompressed U*S*V^T).
%   DOFl         : Indices of free (unconstrained) DOFs.
%   ind_elastic  : Column indices of elastic snapshots.
%   ind_plastic  : Column indices of plastic snapshots.
%   DATAoffline  : Options / tolerances:
%                    • errorDISP  – WSVDT/SVD truncation tolerance.
%                    • nmodes_PLASTIC – cap on # plastic modes.
%                    • NormToUseForSVD_truncation_Displacements ∈
%                          {'StiffnessMatrix','GeometricMassMatrix'}
%                      selects the weighting used to truncate plastic modes.
%   DATA_interp  : Interpolation settings for slave mapping (Bsplines):
%                    • METHOD_SELECT_REFERENCE_MODE ∈ {'FIRST_MODE','FIRST_SVD_MODE','STANDARD_ROM'}
%                    • order_Bsplines, NSAMPLES, ratio_NSAMPLES_knots, …
%   OTHER_output : FE/meta data (K, BC mapping A, plotting buffers, etc.).
%   OPERFE       : FE operators; may include geometric mass matrix OPERFE.Mgeo.
%   INTV_GLO     : Global internal-variable snapshots (kept for alternative
%                  latent definitions; defaults here to FIRST_SVD_MODE path).
%
% OUTPUTS
%   PhiMaster_lin   : Elastic basis (Kll-orthonormal on DOFl).
%   PhiMaster_nonl  : Master plastic mode (first plastic SVD mode).
%   PhiSlave_nonl   : Plastic slave modes.
%   qMASTER_nonl    : Coordinates along PhiMaster_nonl.
%   qSLAVE_nonl     : Coordinates along PhiSlave_nonl.
%   OTHER_output    : Updated (adds Phi_To_Plot for GiD visualization, etc.).
%   DATA_evaluateTAU_and_DER : Handles for τ(q), τ′(q), τ″(q).
%   nREDcoor        : # generalized coordinates used by the (ELAST ⊕ master) model.
%   Kll             : Constrained stiffness (or periodic) metric on DOFl.
%   DATA_interp     : Potentially enriched with IndexLatentPlasticVariable, nREDcoor.
%
% METHOD (high level)
%   1) Build Kll:
%        – No affine map (A empty): Kll = K(DOFl,DOFl).
%        – With affine map (d = A*dL): Kll = Aᵀ K A  (periodic/affine BCs).
%   2) Elastic basis:
%        – WSVDT on SNAPdisp(DOFl,ind_elastic) with Kll → PhiMaster_lin.
%   3) Plastic residual:
%        – Project plastic snapshots Kll-orthogonal to PhiMaster_lin.
%   4) Plastic basis:
%        – If DATAoffline.NormToUseForSVD_truncation_Displacements == 'StiffnessMatrix'
%             → PlasticModesNMODESgiven(... Kll_chol, Kll)
%          Else if 'GeometricMassMatrix' and OPERFE.Mgeo available
%             → PlasticModesNMODESgivenM(... Kll_chol, Kll, MgeoLL)
%        – Truncate to nmodes_PLASTIC.
%   5) Split into master/slave plastic modes, compute qMASTER/qSLAVE via Kll projections.
%   6) If METHOD_SELECT_REFERENCE_MODE ∈ {'FIRST_MODE','FIRST_SVD_MODE'}
%        → Fit B-splines qSLAVE = f(qMASTER) and build τ, τ′, τ″ (Bsplines_EP_1stmode).
%      If 'STANDARD_ROM'
%        → pass-through (identity τ), i.e., no manifold coupling.
%
% DIAGNOSTICS / SIDE EFFECTS
%   • Prints counts of elastic/plastic modes kept.
%   • OTHER_output.Phi_To_Plot stores [Φ_ELAST, Φ_PLAST] expanded to full DOFs
%     (via A if affine BCs) for GiD visualization.
%
% NOTES / GOTCHAS
%   • Master = first plastic SVD mode is REQUIRED for a viable encoder; “linear
%     with force”, “arc-length”, etc. are present as experiments but not viable.
%   • DATAoffline.nmodes_PLASTIC limits plastic basis size regardless of tol.
%   • Choosing 'GeometricMassMatrix' requires OPERFE.Mgeo; otherwise fallback
%     to 'StiffnessMatrix'.
%
% VERSION / AUTHORSHIP
%   17-AUG-2025 — J.A. Hernández, Molinos Marfagones (Cartagena).
%   07-NOV-2025 — Comments refreshed/clarified; K/M-weighted truncation noted; 
%                  encoder viability emphasized. Barcelona.
%   Author: Joaquín Alberto Hernández Ortega (JAHO) — jhortega@cimne.upc.edu
%   (This header was generated automatically by ChatGPT on 07-NOV-2025.)
% =========================================================================


%--------------------------------------------------------------------------

if nargin == 0
    load('tmp1.mat')
    % DATA_interp.NSAMPLES = 100;
    %s close all
    
end
OTHER_output = DefaultField(OTHER_output,'DISP_CONDITIONS') ;
OTHER_output.DISP_CONDITIONS = DefaultField( OTHER_output.DISP_CONDITIONS,'A',[]) ;
if isempty(OTHER_output.DISP_CONDITIONS.A)
    Kll = OTHER_output.K(DOFl,DOFl) ;
else
    Kll = OTHER_output.DISP_CONDITIONS.A'*(OTHER_output.K*OTHER_output.DISP_CONDITIONS.A) ;
    % https://chatgpt.com/share/68e88e70-8f74-8013-afc3-e184680459bb
    % A corresponds to a mapping matrix such that d = A*dL; here dL are the
    % independent DOFs. Kll is the "periodic" stiffness.
end
OPERFE = DefaultField(OPERFE,'Mgeo',[]) ; 
if ~isempty(OPERFE.Mgeo)
    if isempty(OTHER_output.DISP_CONDITIONS.A)
        MgeoLL = OPERFE.Mgeo(DOFl,DOFl) ;
    else
        MgeoLL = OTHER_output.DISP_CONDITIONS.A'*(OPERFE.Mgeo*OTHER_output.DISP_CONDITIONS.A) ;
        
    end
    
else
    MgeoLL = [] ;
end



% LET US BEGIN BY DETERMINING the elastic (linear) mode
% For consistency, we shall use the norm induced by the constrained
% stiffness matrix
DATAs.TOL = 1e-10 ;
[PhiMaster_lin,SS,VV,Kll_chol] = WSVDT(SNAPdisp(DOFl,ind_elastic),Kll,DATAs);
% if length(SS) >1
%     error('This function only supports 1-parameter problems')
% end
% To be plotted in GID
%DOFr = OTHER_output.DISP_CONDITIONS.DOFr;
if ~isempty(OTHER_output.DISP_CONDITIONS.A)
    %  DOFm = OTHER_output.DISP_CONDITIONS.DOFm;
    PhiMaster_lin_plot = OTHER_output.DISP_CONDITIONS.A*PhiMaster_lin ;
else
    PhiMaster_lin_plot = zeros(size(SNAPdisp,1),size(PhiMaster_lin,2)) ;
    
    PhiMaster_lin_plot(DOFl,:) = PhiMaster_lin ;
    
end
% ---------------------------------------------------------------------------------

% % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/15_FORCE_plast.mlx
 DATAoffline = DefaultField(DATAoffline,'NormToUseForSVD_truncation_Displacements','StiffnessMatrix') ; %GeometricMassMatrix
% 
switch DATAoffline.NormToUseForSVD_truncation_Displacements
    case 'GeometricMassMatrix'
   [PhiALL,Phi_non,SNAPdisp_plast_L] =PlasticModesNMODESgivenM(PhiMaster_lin,SNAPdisp,DOFl,ind_plastic,DATAoffline,Kll_chol,Kll,MgeoLL) ;
    case 'StiffnessMatrix'
        % Before 18-Oct-2025
        [PhiALL,Phi_non,SNAPdisp_plast_L] =PlasticModesNMODESgiven(PhiMaster_lin,SNAPdisp,DOFl,ind_plastic,DATAoffline,Kll_chol,Kll) ;
    otherwise
        error('Option not implemented')
end




if ~isempty(OTHER_output.DISP_CONDITIONS.A)
    %  DOFm = OTHER_output.DISP_CONDITIONS.DOFm;
    Phi_non_plot = OTHER_output.DISP_CONDITIONS.A*Phi_non ;
else
    Phi_non_plot = zeros(size(SNAPdisp,1),size(Phi_non,2)) ;
    
    Phi_non_plot(DOFl,:) = Phi_non ;
    
end

OTHER_output.Phi_To_Plot = [PhiMaster_lin_plot,Phi_non_plot] ;



disp('***********************************************************')
disp(['Number of elastic displacement modes =',num2str(size(PhiMaster_lin,2)), ' '])
disp('***********************************************************')

disp('***********************************************************')
disp(['Number of inelastic displacement modes =',num2str(size(Phi_non,2)) ')'])
disp('***********************************************************')

DATA_interp = DefaultField(DATA_interp,'METHOD_SELECT_REFERENCE_MODE','FIRST_MODE') ;
DATA_interp.ChooseTypeLatentVariable_plasticity  = DATA_interp.METHOD_SELECT_REFERENCE_MODE ;
switch DATA_interp.ChooseTypeLatentVariable_plasticity
    case 'STANDARD_ROM'
        PhiMaster_nonl = Phi_non ;
        PhiSlave_nonl = [] ;
        qMASTER_nonl = [] ;
        qSLAVE_nonl = [] ;
        nREDcoor = size(PhiALL,2);
        DATA_evaluateTAU_and_DER.nameFunctionEvaluate = 'tauFUN_identity';
    case {'FIRST_MODE','FIRST_SVD_MODE'}
        
                nREDcoor = size(PhiMaster_lin,2) + 1;  

        DATA_interp.IndexLatentPlasticVariable = size(PhiMaster_lin,2) + 1; 
        DATA_interp.nREDcoor = nREDcoor ;  
        [PhiMaster_nonl,PhiSlave_nonl,qMASTER_nonl,...
            qSLAVE_nonl,...
            DATA_evaluateTAU_and_DER] = Bsplines_EP_1stmode(Phi_non,Kll,DATA_interp,SNAPdisp_plast_L) ;
        
        
        %     case 'LINEAR_WITH_FORCE'
        %         error('Option not viable')
        %
        %
        %     case 'ArcLength'
        %       error('Option not viable')
        % %         % OPTION 1
        %         % LATENT VARIABLE FOR THE PLASTIC COMPONET: VOLUME-AVERAGE OF THE INTERNAL VARIABLE
        %         % sEE /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/FE_CODE/NONLINEAR_Dyn/SmallStrainJ2Plasticity_LOCAL.m
        %         % alphaN   - [ngaus × 1] accumulated plastic strain at step n
        %         %qMASTER_nonl =   sum(INTV_GLO(:,ind_plastic).*OPERFE.wSTs,1)/sum(OPERFE.wSTs) ;
        %
        %         %error('Option not implemented')
        %         PhiMaster_nonl = [] ; % Phi_non(:,1) ;  % Master mode, plastic
        %         PhiSlave_nonl =  Phi_non ; % Slave modes, plastic
        %
        %         qSLAVE_nonl = PhiSlave_nonl'*(Kll*SNAPdisp_plast_L) ;
        %
        %         % OPTION 2
        %         % https://chatgpt.com/share/689f2694-af88-8013-a554-0d33dd96253a
        %         % Arclength
        %         d = vecnorm(diff(qSLAVE_nonl,1,2),2,1);      % Euclidean norm per segment
        %         s = [0, cumsum(d)];                  % cumulative
        %         qMASTER_nonl = s / s(end);                   % normalize to [0,1]
        %
        %
        %         figure(35)
        %         hold on
        %         xlabel('qMASTER')
        %         ylabel('q')
        %
        %         for iii  =1:size(qSLAVE_nonl,1)
        %             plot(qMASTER_nonl, qSLAVE_nonl(iii,:),'DisplayName',['q',num2str(iii)])
        %         end
        %         legend show
        %
        %                 q_firstMODE = qSLAVE_nonl(1,:) ;
        %
        %
        %         [UU,SS,VV] = SVDT(qSLAVE_nonl) ;
        %
        %         % figure(345)
        %         % hold on
        %         % title('Right singular vectors qSLAVE nonl as a function of qMASTER nonl')
        %         % xlabel('qMASTER nonl')
        %         % ylabel('qSLAVE nonl')
        %         %
        %         %
        %         % for iii = 1:size(VV,2)
        %         %     plot(qMASTER_nonl,VV(:,iii),'DisplayName',['\lambda',num2str(iii),'=',num2str(SS(iii)/SS(1))]) ;
        %         % end
        %
        %
        %         [DATA_evaluateTAU_and_DER, nREDcoor] = BsplinesLeastSquares_PLAST_al(DATA_interp, qMASTER_nonl, VV', UU, SS,q_firstMODE) ;
        %
        %
        %
        %
    otherwise
        error(['Option not implemented...'])
        
end




% Since we intend to construct qSLAVE_nonl = f(qMASTER_nonl) using

