function  [z_mst, z_slv,w_mst,w_slv,DATA_regress_eta_der] = DiscreteECM_adaptWEIGHTS_ELPLAST2(A,z,w,DATA_interp,wFE)
%--------------------------------------------------------------------------
% DiscreteECM_adaptWEIGHTS_ELPLAST2
%
% PURPOSE
%   Adaptive Empirical Cubature Method (ECM) tailored to elastoplastic
%   problems driven by a *single loading parameter* but evolving on a
%   *two-latent-variable* manifold (e.g., plastic strain + hardening var.).
%   Among the current ECM points, the algorithm selects one “master” point
%   whose integrand value is most explanatory, and reconstructs the
%   remaining “slave” points as a *nonlinear* function of the master via a
%   spline regression built on a weighted SVD (WSVD) of the slave data.
%
% CONTEXT / WHAT CHANGED VS. DiscreteECM_adaptWEIGHTS
%   • Designed for elastoplastic responses split into components associated
%     with first/second derivatives (or, more generally, two latent
%     coordinates). Internally the snapshot matrix restricted to ECM points
%     is partitioned as:
%          A_zGLO = { A_z(:,1:2:end),  A_z(:,2:2:end) }.
%   • Returns a compact *evaluator struct* (DATA_regress_eta_der) per
%     component instead of explicit function handles, enabling fast,
%     vectorized calls to evaluate the slave reconstruction from the master.
%
% METHOD (per component i = 1..2)
%   1) Exact/reference integral:  I_ref = A' * wFE
%      Standard ECM estimate:     I_ecm = A(z,:)' * w
%      Report relative error ||I_ecm − I_ref|| / ||I_ref||.
%   2) Choose the master ECM row index:
%        - If DATA_interp.IndexMasterPointGiven is provided, use it.
%        - Otherwise, call pick_most_monotone_injective(A_z) to favor rows
%          whose evolution across snapshots is (piecewise) monotone and as
%          injective as possible (proxy for a well-posed master map).
%   3) Split ECM into master/slaves:
%        A_mst = A_z(i_mst,:),
%        A_slv = A_z(i_slv,:), with corresponding weights w_mst, w_slv.
%   4) Sort snapshots by A_mst (enforces a well-ordered abscissa for
%      regression) and apply a *weighted* SVD to the slave block:
%        [U_s, S_s, V_s] = WSVDT(A_slv, diag(w_slv), DATA_SVD)
%      DATA_SVD.TOL = DATA_interp.LocalToleranceECMnonlinear_masterPOINTS
%      (relative truncation if DATA_SVD.RELATIVE_SVD = 1).
%   5) Regress V_s^T(U_s S_s) coefficients as smooth functions of A_mst
%      with least-squares *B-splines*:
%        DATA_regress_eta_der{i} = BsplinesLeastSquares_fastECM(...)
%      which produces an evaluator to map A_mst → A_slv (nonlinear η).
%   6) Assemble the adaptive ECM estimate using the reconstructed slaves:
%        A_slv_hat = feval(DATA_regress_eta_der{i}.nameFunctionEvaluate, ...
%                          A_mst, DATA_regress_eta_der{i});
%        I_adapt = A_mst' * w_mst + A_slv_hat * w_slv
%      Report the relative deviation vs. the original ECM component:
%        ||I_adapt − (A_z' * w)(sorted)|| / ||A_z' * w||.
%
% INPUTS
%   A           : Cell array or matrix of integrand snapshots (M × Nsnap).
%                 If cell, it is concatenated with cell2mat(A).
%   z           : Indices (rows) of the selected ECM points (length p).
%   w           : ECM weights aligned with z (length p).
%   DATA_interp : Struct with interpolation/SVD settings. Recognized fields:
%                 - IndexMasterPointGiven (scalar or empty)
%                 - LocalToleranceECMnonlinear_masterPOINTS (default 1e-5)
%                 - NSAMPLES (default propagated to NKNOTS_ECMadaptive)
%                 - NKNOTS_ECMadaptive (no. of spline samples/knots)
%                 - (optional) plotting/label fields used by the helpers.
%   wFE         : Full FE quadrature weights (M × 1), used for I_ref.
%
% OUTPUTS
%   z_mst       : ECM *global* index of the master point (scalar).
%   z_slv       : ECM *global* indices of the slave points (vector length p−1).
%   w_mst       : Weight of the master point (scalar).
%   w_slv       : Weights of the slave points (vector length p−1).
%   DATA_regress_eta_der
%               : 2×1 cell (one per component) with spline-regression
%                 evaluators. Each cell is a struct suitable for:
%                    A_slv_hat = feval( .nameFunctionEvaluate, A_mst, . )
%                 Typical fields (implementation-dependent):
%                    .nameFunctionEvaluate  : function name (char)
%                    .knots, .degree, .coeff: spline data
%                    .scales, .shifts       : normalization metadata
%                    .rSVD                  : retained WSVD rank
%
% ASSUMPTIONS / PRACTICAL NOTES
%   • Single loading parameter (snapshots ordered consistently in “time” or
%     load), but response evolves on a 2-latent manifold (elastoplasticity).
%   • The master trace A_mst should be (piecewise) monotone/injective; if
%     not, slave regression may be ill-posed. Provide
%     DATA_interp.IndexMasterPointGiven to override automatic selection.
%   • Sorting by A_mst is mandatory before spline fitting.
%   • WSVD is essential: weighs samples by w_slv to honor the cubature
%     measure during basis extraction.
%   • For very stiff/non-smooth transitions (yield onset), increase the
%     number of spline knots (NKNOTS_ECMadaptive) and/or relax the SVD
%     tolerance to avoid underfitting/overtruncation.
%
% DIAGNOSTICS
%   • Prints: standard ECM error vs. full FE; adaptive ECM error vs.
%     standard ECM for each component.
%   • Optional plots: (i) raw ECM traces across snapshots; (ii) right
%     singular vectors V_s against A_mst (to gauge smoothness).
%
% LIMITATIONS / EDGE CASES
%   • If A_mst is flat or has strong backtracks (non-injective), automatic
%     selection may fail—set IndexMasterPointGiven manually.
%   • Multi-parameter load paths or more than two latent variables require
%     extending the partitioning and/or using multi-variate regression.
%
% EXAMPLE (schematic)
%   [zm, zs, wm, ws, eta_data] = DiscreteECM_adaptWEIGHTS_ELPLAST2(A, z, w, DATA, wFE);
%   A_mst_new  = A(zm,:);                         % master trace
%   A_slv_hat1 = feval(eta_data{1}.nameFunctionEvaluate, A_mst_new, eta_data{1});
%   A_slv_hat2 = feval(eta_data{2}.nameFunctionEvaluate, A_mst_new, eta_data{2});
%
% AUTHOR / VERSION
%   J.A. Hernández (UPC/CIMNE)
%   26 Aug 2025, Honest Greens – Tuset, Barcelona
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------

%--------------------------------------------------------------------------

if nargin == 0
    load('tmpPOINTS_SEP_el.mat')
    % DATA_interp.METHOD_SELECT_MASTER_ECM_POINTS = 'MAXIMUM_WEIGHT' ;
    DATA_interp.LocalToleranceECMnonlinear_masterPOINTS = 1e-3;
    %  DATA_interp.NSAMPLES = 500 ;
    close all
end

A = cell2mat(A) ;
ExactInt = A'*wFE;
A_z = A(z,:);
Int_ECM_stand = A(z,:)'*w ;

ERRO_ECMstand =  norm(Int_ECM_stand-ExactInt)/norm(ExactInt) ;
disp(['Error Standard ECM = ',num2str(ERRO_ECMstand*100),' %'])


% First step: separate the contributions of the first derivative and the
% second derivative
A_zGLO = {A_z(:,1:2:end),A_z(:,2:2:end)}  ;

% for iii = 1:length(A_zGLO)
%     [UU,SS,ZZ] = SVDT(A_zGLO{iii});
%     rLOCAL = length(SS) ;
%     disp(['Rank submatrix A_',num2str(iii),'=',num2str(rLOCAL)]) ;
%
% end
disp(['Number of ECM points = ',num2str(length(z))])




% NOW WE CONFRONT THE MOST DELICATE PHASE OF THE ENTIRE "ADAPTIVE" ECM: WE
% HAVE TO ANSWER THE QUESTION:  IS THERE ANY SINGLE POINT (THE MASTER
% POINT) such that we can express the remaining variables as nonlinear
% functions of the integrand at this master point ?


%DATA_interp = DefaultField(DATA_interp,'METHOD_SELECT_MASTER_ECM_POINTS','MOST_EXPLANATORY') ;

% switch DATA_interp.METHOD_SELECT_MASTER_ECM_POINTS
%     case 'MOST_EXPLANATORY'
DATA_regress_eta_der = cell(length(A_zGLO),1) ;

%DATA_interp.IndexMasterPointGiven = 35;

DATA_interp = DefaultField(DATA_interp,'IndexMasterPointGiven',[]) ; 

% Selected in /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/08_PLAST_adapECM.mlx
% using /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/08_PLAST_adapECM/INSPECTING_Afint_uSEP.m
DATA_interp = DefaultField(DATA_interp,'NKNOTS_ECMadaptive', DATA_interp.NSAMPLES) ;
DATA_interp = DefaultField(DATA_interp,'LocalToleranceECMnonlinear_masterPOINTS',1e-5);

for iLOC = 1:length(A_zGLO)
    
    
    disp(['-----------------------------------------------'])
    disp(['Analysis component FintRED number = ',num2str(iLOC)])
    disp(['-----------------------------------------------'])
    
    A_z = A_zGLO{iLOC} ;
    
    if ~isempty(DATA_interp.IndexMasterPointGiven)
        best_col_index =DATA_interp.IndexMasterPointGiven ;
        disp(['Master Point for both components = ',num2str(z(best_col_index))])
        disp(['( selected by an auxiliar program) '])
    else
        [best_col_index, scores,report] = pick_most_monotone_injective(A_z) ;
        
        
        
        
        disp(['SELECTED POINT = ',num2str(best_col_index)])
        disp(['IS MONOTONIC?, CHECKING INDEX MR (IF = 1, IT IS MONOTONIC) =  ',num2str(scores(best_col_index).MR)])
        
        if abs(scores(best_col_index).MR - 1) <= 1e-14
            
            disp('It is indeed monotonic !')
            disp(['Is it injective ? '])
            disp(report.is_injective)
            %  disp([])
            
        end
    end
    
    
    compl_Set = setdiff(1:size(A_z,1),best_col_index) ;
    
    A_mst = A_z(best_col_index,:) ;
    A_slv = A_z(compl_Set,:) ;
    z_mst = z(best_col_index) ;
    w_mst = w(best_col_index) ;
    z_slv = z(compl_Set);
    w_slv = w(compl_Set) ;
    
    PLOT_EVOLUTION_Az = 1 ;
    
    if  PLOT_EVOLUTION_Az == 1
        disp('WARNING: The following graphical representation is useless (or misleading) in more than one latent dimension')
        figure(3094 + iLOC)
        hold on
        
        title(['Internal work density, component = ',num2str(iLOC)])
        xlabel('Snapshot ')
        ylabel('Internal work density')
        numberECMpoints = size(A_z,1);
        % numberECMpoints = min(numberECMpoints, size(A_z,1)) ;
        wREL = w/sum(w)*100 ;
        
        
        for ipoints = 1:numberECMpoints
            y = A_z(ipoints,:)  ;
            plot(y,'DisplayName',['z_{',num2str(ipoints),'} =', num2str(z(ipoints)),'  w % = ',num2str(wREL(ipoints))])  ;
        end
        
        legend show
        
        
    end
    
    % 5-08-2025:  Why did we introduce this in the first place?
    % See discussion in /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/05_meta_AECM.mlx
    
    % = DefaultField(DATA_interp,'METHOD_SORT_ENTRIES_MASTER_POINTS','monotonic') ;
    
    % switch DATA_interp.METHOD_SORT_ENTRIES_MASTER_POINTS
    %     case 'monotonic'
    [A_mst,indMST] = sort(A_mst) ;
    A_slv = A_slv(:,indMST) ;
    
    
    
    
    
    
    
    
    
    
    
    
    
    TOL = DATA_interp.LocalToleranceECMnonlinear_masterPOINTS;
    M  =diag(w_slv) ;
    DATA_SVD.RELATIVE_SVD = 1;
    DATA_SVD.TOL = TOL ;
    % This is the weighted SVD,
    [U_s,S_s,V_s] = WSVDT(A_slv,M,DATA_SVD) ;
    
    
    PLOT_BEFORE_SPLINES = 1;
    
    if  PLOT_BEFORE_SPLINES == 1
        figure(463+iLOC)
        hold on
        title(['Right singular vectors SVD Amaster as function of Aslave, component = ',num2str(iLOC)]) ;
        xlabel('A_{mst}')
        ylabel('V_{slv}')
        for iii  =1:size(V_s,2)
            vs = S_s(iii)/S_s(1) ;
            plot(A_mst,V_s(:,iii),'DisplayName',['\lambda = ',num2str(vs)] )
        end
    end
    legend show
    
    DATA_interpLOC = DATA_interp ;
    
    DATA_interpLOC.Nfigures_base = 450+iLOC ;
    
    DATA_interpLOC = DefaultField(DATA_interpLOC,'Legend_var_x','A_{mst}') ;
    DATA_interpLOC = DefaultField(DATA_interpLOC,'Legend_var_y','\eta  (A_{slv})') ;
    DATA_interpLOC = DefaultField(DATA_interpLOC,'TITLE_additional',[' COMP =',num2str(iLOC)]) ;
    DATA_interpLOC = DefaultField(DATA_interpLOC,'NSAMPLES', DATA_interpLOC.NKNOTS_ECMadaptive) ;
    
    % [etaNON, etaNONder,etaNONder2] = BsplinesLeastSquares_ECM(DATA_interp, A_mst, V_s', U_s, S_s) ;
    [DATA_regress_eta_der{iLOC}, nREDcoor] = BsplinesLeastSquares_fastECM(DATA_interpLOC, A_mst,  V_s', U_s, S_s) ;
    
    %  [DATA_evaluateTAU_and_DER, nREDcoor] = BsplinesLeastSquares_fast(DATA_interp, qINF, VrightT, UU, SS) ;
    
    
    %[etaNON, etaNONder,etaNONder2] = BsplinesLeastSquares(DATA_interp, A_mst, V_s', U_s, S_s) ;
    
    %[etaNON, etaNONder,nREDcoor] = SplineInterpLOC_ECM(DATA_interp, A_mst, V_s', U_s, S_s) ;
    
    
    %Int_AdaptiveECM =  etaNON(A_mst)'*[w_mst; w_slv] ;
    
    A_slv = feval(DATA_regress_eta_der{iLOC}.nameFunctionEvaluate,A_mst,DATA_regress_eta_der{iLOC}) ;
    
    Int_AdaptiveECM =A_mst'*w_mst +  A_slv*w_slv ;
    
    ExactInt_ECM = A_zGLO{iLOC}'*w ;
    
    ERRO_ECMnew =  norm(Int_AdaptiveECM-ExactInt_ECM(indMST))/norm(ExactInt_ECM) ;
    disp(['Error adaptive ECM with respect ECM = ',num2str(ERRO_ECMnew*100),' %'])
    
    
    
end
