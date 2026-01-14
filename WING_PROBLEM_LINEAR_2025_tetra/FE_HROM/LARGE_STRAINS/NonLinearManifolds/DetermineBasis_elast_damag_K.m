function [PhiMaster_lin,PhiMaster_nonl,PhiSlave_nonl,qMASTER_nonl,...
    qSLAVE_nonl,OTHER_output,...
    DATA_evaluateTAU_and_DER, nREDcoor,DATA_interp,GmatN] = DetermineBasis_elast_damag_K(SNAPdisp,DOFl,ind_linear,ind_nonlinear,...
    DATAoffline,DATA_interp,OTHER_output,OPERFE,INTV_GLO)
% DetermineBasis_elast_damag_K
% -------------------------------------------------------------------------
% PURPOSE
%   Build the *damage-aware* reduced basis and the associated 2-coordinate
%   manifold ingredients using displacement snapshots and a weighted inner
%   product G (typically the constrained stiffness K_ll). The routine:
%     (i)   extracts a 1-mode *elastic/linear* basis Φ_lin,
%     (ii)  extracts a *nonlinear/damage* basis Φ_non in the G-orthogonal
%           complement of Φ_lin,
%     (iii) designates a *nonlinear master* mode Φ_non,mst and its *slave*
%           complement Φ_non,slv,
%     (iv)  computes G-weighted coordinates (qMASTER_nonl, qSLAVE_nonl),
%     (v)   fits B-splines so that slave amplitudes become smooth functions
%           of the master coordinate (for τ(q), τ′(q), τ″(q) assembly).
%
% THEORETICAL CONTEXT (univariate damage manifold; doc pp. 138–146)
%   The decoder uses two latent coordinates q = (qLIN, qNON) and factors the
%   response as
%       τ(q) = [ qLIN ; qLIN*qNON ; qLIN * g(qNON) ],
%   where g(·) stacks slave amplitudes in a compact SVD basis. This function
%   prepares Φ_lin, Φ_non,mst/slv and the spline data so that τ(q), ∂τ/∂q and
%   ∂²τ/∂q² can be evaluated consistently (e.g., by tauFUN_1paramDAMAGE).
%
% WORKFLOW
%   1) Choose inner product G:
%        • 'StiffnessMatrix' → G = K_ll,
%        • 'GeometricMassMatrix' → G = M_geo,ll (if provided).
%   2) Linear mode (elastic):
%        Φ_lin = WSVDT( SNAPdisp(DOFl,ind_linear), G ), enforced rank = 1.
%   3) Nonlinear (damage) subspace:
%        Extract Φ_non in the G-orthogonal complement of Φ_lin
%        (NonlinearModesNMODESgivenM, truncated by DATAoffline settings).
%   4) Plotting pack:
%        Store bases mapped to full DOFs (periodic/constraint map A if used)
%        in OTHER_output.Phi_To_Plot for GiD visualization.
%   5) Master/slave split and spline learning (default FIRST_MODE):
%        • Set nREDcoor = 1 (linear) + 1 (nonlinear master).
%        • Call Bsplines_EP_1stmodeDMG(Φ_lin, Φ_non, G, ...):
%            – returns Φ_non,mst, Φ_non,slv,
%              qMASTER_nonl (≡ qNON) and qSLAVE_nonl,
%              DATA_evaluateTAU_and_DER with spline handles for g, g′, g″.
%        Alternative latent-variable options are retained only for reference
%        and are marked “not viable” for encoder construction.
%
% INPUTS
%   SNAPdisp     : [nDOF×nsnap] displacement snapshots.
%   DOFl         : Indices of independent/free DOFs.
%   ind_linear   : Snapshot indices for (nearly) elastic states.
%   ind_nonlinear: Snapshot indices for inelastic/damage evolution.
%   DATAoffline  : Offline settings (norm selection, nmodes limits, tolerances).
%   DATA_interp  : Spline/learning options; includes method selector for the
%                  nonlinear master ('FIRST_MODE' / 'FIRST_SVD_MODE' default).
%   OTHER_output : Struct holding FE operators (e.g., K) and optional mapping A.
%   OPERFE       : Optional FE data (e.g., Mgeo for geometric mass inner product).
%   INTV_GLO     : (Reserved) internal-variable snapshots for alternative masters.
%
% OUTPUTS
%   PhiMaster_lin  : [nDOF×1] G-orthonormal elastic/linear mode Φ_lin.
%   PhiMaster_nonl : [nDOF×1] nonlinear *master* damage mode Φ_non,mst.
%   PhiSlave_nonl  : [nDOF×(r_nl−1)] G-orthonormal nonlinear *slave* modes.
%   qMASTER_nonl   : [1×n_snap] master coordinate along Φ_non,mst.
%   qSLAVE_nonl    : [(r_nl−1)×n_snap] slave coordinates along Φ_non,slv.
%   OTHER_output   : Updated (adds Phi_To_Plot, messages on mode counts, etc.).
%   DATA_evaluateTAU_and_DER :
%                    Spline evaluators/metadata to assemble τ, τ′, τ″.
%   nREDcoor       : Reduced dimension used by the damage manifold (usually 2).
%   DATA_interp    : Echoed back with inferred fields (e.g., index of qNON).
%   GmatN          : The inner-product matrix actually used (K_ll or M_geo,ll).
%
% NOTES & DIAGNOSTICS
%   • All bases and coordinates are computed in the G-inner product. When
%     constraints are active (d = A d_L), the code works with “periodic”
%     forms K_ll = Aᵀ K A (and analogously for M_geo).
%   • Mode counts and basic diagnostics are printed to console.
%   • DATAoffline.nmodes_PLASTIC (or analogous) caps the nonlinear basis size.
%   • This routine is agnostic to hyper-reduction; it only prepares the
%     manifold ingredients required by subsequent τ-evaluation and tangents.
%
% VERSION / AUTHOR
%   23-Oct-2025 — J.A. Hernández (Barcelona). Damage-adapted basis builder
%   with K_ll (or M_geo,ll) inner product and first-mode nonlinear master.
% -------------------------------------------------------------------------

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

 DATAoffline = DefaultField(DATAoffline,'NormToUseForSVD_truncation_Displacements','StiffnessMatrix') ; %GeometricMassMatrix
% 

switch DATAoffline.NormToUseForSVD_truncation_Displacements
    case 'GeometricMassMatrix'
        GmatN = MgeoLL ; MgeoLL = [] ;
        disp('Using geometric mass matrix  (G= MgeoLL)  for orthogonality')
    case 'StiffnessMatrix'
        % Before 18-Oct-2025
        GmatN = Kll ; 
        disp('Using stiffness    matrix (G =Kll) for orthogonality')
     otherwise
        error('Option not implemented')
end


% LET US BEGIN BY DETERMINING the elastic (linear) mode
% For consistency, we shall use the norm induced by the constrained
% stiffness matrix
DATAs.TOL = 1e-10 ;
[PhiMaster_lin,SS,VV,GmatN_chol] = WSVDT(SNAPdisp(DOFl,ind_linear),GmatN,DATAs);
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
 % 
% switch DATAoffline.NormToUseForSVD_truncation_Displacements
%     case 'GeometricMassMatrix'
%         
   [PhiALL,Phi_non] =...
       NonlinearModesNMODESgivenM(PhiMaster_lin,SNAPdisp,DOFl,ind_nonlinear,DATAoffline,GmatN,GmatN_chol) ;
%     case 'StiffnessMatrix'
%         % Before 18-Oct-2025
%         [PhiALL,Phi_non] =PlasticModesNMODESgiven(PhiMaster_lin,SNAPdisp,DOFl,ind_nonlinear,DATAoffline,Kll_chol,Kll) ;
%     otherwise
%         error('Option not implemented')
% end




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

        DATA_interp.IndexLatentNonlinearVariable = size(PhiMaster_lin,2) + 1; 
        DATA_interp.nREDcoor = nREDcoor ;  
%         [PhiMaster_nonl,PhiSlave_nonl,qMASTER_nonl,...
%             qSLAVE_nonl,...
%             DATA_evaluateTAU_and_DER] = ...
%             Bsplines_EP_1stmodeDMG(PhiALL(:,1),Phi_non,GmatN,DATA_interp,SNAPdisp(DOFl,ind_nonlinear),GmatN_chol) ;
        
         [PhiMaster_nonl,PhiSlave_nonl,qMASTER_nonl,...
            qSLAVE_nonl,...
            DATA_evaluateTAU_and_DER] = ...
            Bsplines_EP_1stmodeDMG_v2(PhiALL(:,1),Phi_non,GmatN,DATA_interp,SNAPdisp(DOFl,ind_nonlinear),GmatN_chol) ;
        
        
        
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
        %         %qMASTER_nonl =   sum(INTV_GLO(:,ind_nonlinear).*OPERFE.wSTs,1)/sum(OPERFE.wSTs) ;
        %
        %         %error('Option not implemented')
        %         PhiMaster_nonl = [] ; % Phi_non(:,1) ;  % Master mode, plastic
        %         PhiSlave_nonl =  Phi_non ; % Slave modes, plastic
        %
        %         qSLAVE_nonl = PhiSlave_nonl'*(Kll*SNAPdisp_nonlin_L) ;
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
        %         [DATA_evaluateTAU_and_DER, nREDcoor] = BsplinesLeastSquares_nonlin_al(DATA_interp, qMASTER_nonl, VV', UU, SS,q_firstMODE) ;
        %
        %
        %
        %
    otherwise
        error(['Option not implemented...'])
        
end




% Since we intend to construct qSLAVE_nonl = f(qMASTER_nonl) using

