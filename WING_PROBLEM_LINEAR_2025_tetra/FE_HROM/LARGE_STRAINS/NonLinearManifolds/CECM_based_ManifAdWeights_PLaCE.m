function  ECMdata = CECM_based_ManifAdWeights_PLaCE(A,DATA,Wfe,DATA_interp,zINI,wINI,...
    qLATENT,LOCALdata,Uel,Qall_w,DATAoffline)
%==========================================================================
%==========================================================================
% CECM_based_ManifAdWeights_PLaCE
%
% PURPOSE
% -------
% Manifold-Adaptive Weights ECM for inelastic problems (plasticity/damage)
% with **explicit nonnegativity (w ≥ 0)** and **exhaustive single-point
% elimination** per iteration. Starting from an ECM rule (zINI, wINI), the
% algorithm repeatedly tests removing each currently active integration
% point, re-solves the constrained weights for all clusters, and accepts the
% trial that **minimizes the worst number of zeros across clusters**, which
% promotes balanced sparsity. It returns the reduced rule and a regressor
% mapping qLATENT → weights. Optionally records a GIF of weight evolution.
%
% CONTEXT & LINEAGE
% -----------------
% • CECM  : Continuous Empirical Cubature Method.
% • ECM   : Classical sparse cubature with constant (state-independent) w.
% • MAW-ECM: State-dependent weights (functions of a latent variable q).
%   This “PLaCE” variant enforces w ≥ 0 explicitly and uses exhaustive
%   elimination to avoid unbalanced per-cluster sparsity.
%
% PROBLEM STATEMENT (per cluster c)
% ---------------------------------
%   Find w(:,c) ≥ 0  such that  U{c}' * w(:,c) = b{c}.
% U{c} includes a **fixed elastic scaffold** (Uel and the constant function)
% and plastic/damage snapshots, optionally denoised by an orthonormal Q in
% the FE inner product defined by Wfe.
%
% INPUTS
% ------
% A            : 1×nC cell. A{c} ∈ ℝ^(nIP×m_c), plastic/damage snapshots for
%                cluster c (columns grouped by internal-force components).
% DATA         : Mesh/meta info. Uses DATA.MESH.ngaus to map points→elements.
% Wfe          : (nIP×1) finite-element integration weights/measures.
% DATA_interp  : Options for regression q → weights (order, extrapolation…).
% zINI, wINI   : Initial ECM candidate indices and associated weights.
% qLATENT      : Latent parameters (e.g., plastic multiplier). Code uses the
%                component given by DATA_interp.IndexLatentPlasticVariable
%                (defaults to 2 for plasticity; see “Plastic/Damage range”).
% LOCALdata    : Struct with logical selection of snapshots:
%                  • ind_plastic   – indices of plastic range   (preferred),
%                  • ind_nonlinear – indices for damage range   (alternative).
%                Also provides/assumes availability of Uel (see below).
% Uel          : Elastic basis to be injected in every cluster basis.
% Qall_w       : Columns spanning modes to build Q orthonormal in Wfe
%                (stabilization/denoising of plastic/damage part).
% DATAoffline  : Optional knobs (overlap flags, tolerances, active-set caps,
%                maximum number of ECM points, etc.).
%
% OPTIONAL LOCAL CONTROLS (internal; created if absent)
% -----------------------------------------------------
% LOCAL.MAKE_GIF     (false)  : Save a GIF of weight evolution.
% LOCAL.GIF_FILE     ('weights_evolution.gif')
% LOCAL.DELAY_TIME   (0.25)   : Seconds per frame in the GIF.
% LOCAL.DPI          (120)    : Export resolution.
% LOCAL.SHOW_LEGEND  (true)
%
% KEY OPTIONS (auto-defaulted if missing)
% ---------------------------------------
% DATAoffline.NumberOfOverlappingClustersBothDirections : 0 (no overlap) or 1.
% DATAoffline.MaxNumberZeros_Active_Set_Loop_MAW_ECM    : Max zeros allowed
%                                                         in explicit w≥0 loop (default 1).
% DATAoffline.MaximumNumberOfECMpoints                   : 0 ⇒ ignored; >0 caps rule size.
%
% OUTPUTS
% -------
% ECMdata.setPoints          : Surviving global quadrature-point indices.
% ECMdata.setElements        : Corresponding element ids (large2smallINCLUDEREP).
% ECMdata.wRED.DATA_regress  : Regressor mapping qLATENT → reduced weights.
% ECMdata.wRED.IndexDOFl_q   : Index of latent component used in regression (=2 by default).
%
% HIGH-LEVEL ALGORITHM
% --------------------
% 0) Build Q orthonormal in Wfe; build fixed scaffold Ufixed = SVDT([Uel, ones]).
% 1) For each cluster c:
%      • Optionally project plastics/damage: A_c ← Q(QᵀA_c); column-scale.
%      • Build basis U{c} = SVDT([Ufixed, A_c]) and target b{c}.
%      • Restrict to zINI; initialize wADAPT(:,c) = wINI.
% 2) While #active points > max #modes across clusters:
%      • Rank active points by total weight R = ∑_c w(:,c); sort ascending.
%      • **Exhaustive try** (in that order): remove one point, then for each
%        cluster solve     minimize ‖w − w_prev‖₂  s.t. U' w = b,  w ≥ 0
%        via an active-set NNLS with equalities.
%      • Discard infeasible trials; among feasible trials commit the one
%        that minimizes the **maximum number of zeros** across clusters.
%      • Store wADAPT snapshot (for optional GIF).
% 3) Regress weights vs qLATENT: DATA_regress = MAW_ECM_regression(…).
% 4) Quality check of the rule: CheckAccurayMAT_ECMl2(…).
% 5) Optionally export a GIF of weight evolution.
%
% NUMERICAL NOTES
% ---------------
% • SVDT orthonormalizes columns under the Euclidean inner product; adding
%   the constant function preserves the total weight (volume) constraint.
% • The redistribution step is a **least-change** update anchored at the
%   previous w to improve stability across iterations and clusters.
% • The selection metric (minimize worst #zeros) prevents over-sparsifying
%   a single cluster and typically yields more robust rules.
%
% COMPLEXITY (rule-of-thumb)
% --------------------------
% Let nP = #active points, r = #exactness modes per cluster, nC = #clusters.
% One outer iteration may attempt up to nP removals; each attempt solves nC
% NNLS-with-equalities problems of size (r × (nP−1)). Using Q-projection and
% SVDT on [Uel, A_c] keeps r modest.
%
% PRACTICAL TIPS / LIMITATIONS
% ----------------------------
% • Exhaustive search is costlier than greedy; prefer it when cluster balance
%   matters or when greedy elimination stalls early.
% • Many infeasible trials? Relax/skip Q-projection, rescale columns of A_c,
%   or enrich Ufixed (e.g., add low-energy elastic/plastic modes).
% • Some viewers loop GIFs regardless of “play once”; export MP4 if needed.
%
% SIDE EFFECTS
% ------------
% • Generates diagnostic figures (internal force evolution).
% • May write a GIF if LOCAL.MAKE_GIF = true.
%
% DEPENDENCIES
% ------------
% SVDT, DefaultField, MAWecmBasisMATRIX_LOCAL_noverlap / _overlap,
% Loop_MAWecmNOENF_cl, Loop_MAWecmNOENFe, nn_update_active_set,
% MAW_ECM_regression, CheckAccurayMAT_ECMl2, large2smallINCLUDEREP.
%
% REFERENCES (CECM / MAW-ECM)
% ---------------------------
% Hernández, J.A., “Continuous Empirical Cubature Method (CECM),”
%   Computer Methods in Applied Mechanics and Engineering (CMAME), (year).
% See internal notes and demos:
%   .../112_NonLIN_ROM_RBF/11_MAW_ECM_plast/LATEX/MAW_ECM_v1.pdf
%   .../112_NonLIN_ROM_RBF/11_MAW_ECM_plast.mlx
%
% AUTHOR / PLACE / DATE
% ---------------------
% J.A. Hernández Ortega (JAHO) — Barcelona — 21-Sep-2025
% (Comments updated by ChatGPT-5 Thinking)
%==========================================================================

%==========================================================================

%==========================================================================


if nargin == 0
    load('tmp1.mat')
    %  DATAoffline.MaxNumberZeros_Active_Set_Loop_MAW_ECM = 0 ;
    close all
    % DATAoffline.MaximumNumberOfECMpoints = 15;
    %  DATAoffline.MaxNumberZeros_Active_Set_Loop_MAW_ECM = 200 ;
    %    DATAoffline.Project_SNAP_withQ__MAW_ECM = 1 ;
    %     DATAoffline.NumberOfOverlappingClustersBothDirections = 1;
    %     DATAoffline.ToleranceOverlappingClusters = 1e-2;
    
    
    %      DATA_interp.NSAMPLES  = 20 ;
    
    
    %   DATAoffline.errorFINT_linearpart =1e-6;
    % DATAoffline.Project_SNAP_withQ__MAW_ECM = 1;
    %   DATA_interp.PortionExtrapolation_plot_WEIGHTS = 0.3;
    %  DATA_interp.Extrapolation_Method_WEIGHTS_projection_active_set =  1 ;
end

if nargin < 12 || isempty(LOCAL)
    LOCAL.MAKE_GIF   = false;                 % if true, build GIF after loop
    LOCAL.GIF_FILE   = 'weights_evolution.gif';
    LOCAL.DELAY_TIME = 0.25;                  % seconds per frame
    LOCAL.DPI        = 120;                   % export resolution
    LOCAL.SHOW_LEGEND = true;                 % legend in the GIF frames
end

% SEPARATING INTO COMPONENTS IN THE LINEAR RANGE AND THE PLASTIC RANGE
%Ael = A(:,LOCALdata.ind_elastic) ;
Aall = cell2mat(A);

nCOMP = size(A{1},2) ;


PLOT_ini = 0;

if PLOT_ini == 1
    figure(156)
    hold on
    xlabel('Snapshot')
    ylabel('Internal force ')
    title(['ELASTIC/PLASTIC  Internal force, comp'])
    
    for icmpo = 1:nCOMP
        
        Indexes_loc = icmpo:nCOMP:size(Aall,2) ;
        IntLoc_exact = Aall(:,Indexes_loc)'*Wfe  ;
        IntLoc_approx= Aall(zINI,Indexes_loc)'*wINI  ;
        plot(IntLoc_exact,'DisplayName',['Exact, comp = ',num2str(icmpo)]) ;
        plot(IntLoc_approx,'DisplayName',['ECM, comp = ',num2str(icmpo)]) ;
    end
    legend show
    
end
IS_DAMAGE_PROBLEM  =false;
if ~isempty(LOCALdata)
    if isfield(LOCALdata,'ind_plastic')
        A = A(LOCALdata.ind_plastic) ; % Plastic range
        
    elseif isfield(LOCALdata,'ind_nonlinear')
        A = A(LOCALdata.ind_nonlinear) ; % Damge range
        IS_DAMAGE_PROBLEM  =true;
    else
        error('Option not implemented')
    end
end
% The elastic component should be included in all local basis
% Accordingly, we apply the SVD to obtain a basis matrix
% TOLloc = 1e-6;
% DATAsvdLOC.HIDE_OUTPUT =  1;
% [Uel,SSS,VVV]  = SRSVD(Ael,TOLloc,DATAsvdLOC)  ;
% % We want also to include the constant function
ONES = ones(size(A{1},1),1) ; % This is a vector of ones. It represents the constant function
% Uel is the elastic basis
Ufixed = SVDT([Uel,ONES]) ;
%Ufixed_w = bsxfun(@times,Ufixed,sqW) ;  % Now Q'*WfeD*Q = I


% PLASTIC (or DAMAGE) PART
% ---------------
% We are going to include only plastic/damage snapshots
% Besides, the latent variable is the second component of qLATENT
if ~isempty(LOCALdata)
    
    %
    if ~IS_DAMAGE_PROBLEM
        % Plasticity problem
        iLAT = 2;
        DATA_interp = DefaultField(DATA_interp,'IndexLatentPlasticVariable',iLAT) ;
        
        
        qLATENT = qLATENT(DATA_interp.IndexLatentPlasticVariable,LOCALdata.ind_plastic)  ;
    else
        % Damage problems
        indNONlinear = 2; indLINEAR = 1;
        % qLATENT =qNONrel
        qLATENT_all = qLATENT ;
        qLATENT = qLATENT(indNONlinear,LOCALdata.ind_nonlinear)./ qLATENT(indLINEAR,LOCALdata.ind_nonlinear) ;
        
        DATA_regress.IndexLinear_damageMODEL = 1;
        DATA_regress.IndexNonlinear_damageMODEL = 2;
        
    end
    
else
    DATA_interp = DefaultField(DATA_interp,'IndexLatentPlasticVariable',1) ;
    
end

[qLATENT,ii] = sort(qLATENT) ;


[~,bbb] = find(abs(qLATENT)<1e-6*max(abs(qLATENT))) ;
if ~isempty(bbb)
    qLATENT(bbb) = [] ;
    ii(bbb) = [] ;
end
A = A(ii) ;



Delta_q = diff(qLATENT) ;
% % Norm of A
% rA = zeros(size(A)) ;
% for iCL = 1:length(A)
%     rA(iCL) = norm(A{iCL},'fro') ;
% end
%
%
% % STEP 1
% % Determine the basis matrices of each cluster using A = A{1},A{2} ...
%
% [qLATENT,ii] = sort(qLATENT) ;
% [aab,bbb] = find(abs(qLATENT)<1e-6*max(abs(qLATENT))) ;
% if ~isempty(bbb)
%     qLATENT(bbb) = [] ;
%     ii(bbb) = [] ;
% end
% A = A(:,ii) ;
%


sqW = sqrt(Wfe) ;
%WfeD = diag(sparse(Wfe)) ;
Q = bsxfun(@times,Qall_w,1./sqW) ;  % Now Q'*WfeD*Q = I
Q  = SVDT(Q) ;
%  disp('Check integration error')
%  Am = cell2mat(A) ;
%
% ApproxERR = Am - Q*(Q'*Am);
% ApproxERR = norm(ApproxERR,'fro')/norm(Am,'fro') ;

% Basis matrix for

% DATAsvdLOC.ISRELATIVE = 1;
% TOL_LOC = 1e-10 ;

DATAoffline = DefaultField(DATAoffline,'NumberOfOverlappingClustersBothDirections',0) ; % = 1;

% Basis matrices
if DATAoffline.NumberOfOverlappingClustersBothDirections == 0
    % No overlapping (before 23-Sept-2502)
    [U,b,wADAPT,DATAoffline,numberMODES] = MAWecmBasisMATRIX_LOCAL_noverlap(Ufixed,zINI,qLATENT,A,DATAoffline,wINI,Wfe,Q) ;
else
    warning('No maintenance given to this function since September 2025')
    %   overlapping
    [U,b,wADAPT,DATAoffline,numberMODES] = MAWecmBasisMATRIX_LOCAL_overlap(Ufixed,zINI,qLATENT,A,DATAoffline,wINI,Wfe,Q) ;
end


DATAoffline = DefaultField(DATAoffline,'MaxNumberZeros_Active_Set_Loop_MAW_ECM',1) ; % = 1;
MaxZerosMAW = DATAoffline.MaxNumberZeros_Active_Set_Loop_MAW_ECM ;

%
VOL  = sum(wINI) ;
% [ hhh,htitle]= PlotAUX_mawecm(LOCAL,qLATENT,wADAPT,VOL) ;
% ---- NEW: history storage (no plotting in-loop) -------------------------
historyW = {};              % cell array; each entry is wADAPT at an iteration
historyTitle = {};          % optional: titles per iteration
% Capture the initial state before reduction starts (as "iteration 0")
historyW{end+1} = wADAPT;
historyTitle{end+1} = sprintf('Weights vs qLATENT — initial (iter 0)');


% Initialization
iCAND = 1:length(wINI);  % Candidate points to be eliminated (local indexes)
jELIM = [] ;
%USE_LEAST_NORM = 1;
k  =1;

DATAoffline= DefaultField(DATAoffline,'MaximumNumberOfECMpoints',0) ; % = 15;

max_nmodes= max(numberMODES) ; % Maximum number of modes
max_nmodes = max(max_nmodes,DATAoffline.MaximumNumberOfECMpoints) ;
%maxIter = 20 ;
while   length(iCAND) > max_nmodes
    disp(['iter = ',num2str(k),', Length initial candidate set = ',num2str(length(iCAND))])
    %     if k == 31
    %         disp('Borrar esto')
    %     end
    % 1) Deciding which points to eliminate
    RRR= sum(wADAPT,2) ;
    ilocPOS = find(RRR>0) ;
    % ilocPOS are the indexes of those points which have not eliminated yet
    % Next we order such points according to which has, on average, the
    % lower weights (ascending order)
    [RRRloc_sort,indSORT] = sort(RRR(ilocPOS));
    
    [wADAPT_tent,iCAND_tent,jELIM_tent,kLOC] = Loop_MAWecmNOENF_cl(U,ilocPOS,iCAND,wADAPT,b,indSORT,jELIM);
    
    if kLOC > length(ilocPOS) && DATAoffline.MaxNumberZeros_Active_Set_Loop_MAW_ECM >0
        disp('The algorithm with no enforcement of w>0 cannot reduce anymore the number of points... ')
        disp('We switch to explicit enforcement')
        [wADAPT_tent,iCAND_tent,jELIM_tent,SOLUTION_FOUND] = ...
            Loop_MAWecmNOENFe(U,ilocPOS,iCAND,wADAPT,b,indSORT,jELIM,MaxZerosMAW,Delta_q);
        if  ~SOLUTION_FOUND
            disp('The algorithm cannot reduce the number of points anymore... ')
            break
        else
            wADAPT = wADAPT_tent ;  iCAND = iCAND_tent ; jELIM = jELIM_tent;
        end
    elseif kLOC <= length(ilocPOS)
        wADAPT = wADAPT_tent ;  iCAND = iCAND_tent ; jELIM = jELIM_tent;
    else
        disp('The algorithm cannot reduce the number of points anymore... ')
        break
    end
    
    % ---- NEW: store snapshot for this iteration (no plotting here) ------
    historyW{end+1} = wADAPT; %#ok<AGROW>
    historyTitle{end+1} = sprintf('Weights vs qLATENT — npoints %d',length(iCAND));
    
    %    hhh = SAWECMlocalPLOT_1(LOCAL,hhh,wADAPT,VOL) ;
    k = k+1;
    
end
RRR= sum(wADAPT,2) ;
ilocPOS = find(RRR>0) ;

disp(['Final number  of points = ',num2str(length(ilocPOS))])
wALL = wADAPT(ilocPOS,:) ;
setPointsALL = zINI(ilocPOS) ;

setElements_cand = large2smallINCLUDEREP(setPointsALL,DATA.MESH.ngaus) ;

DATA_interp.setElements =setElements_cand ;
% CONSTRUCT THE REGRESSION MAPPING FROM qLATENT to wALL
% if IS_DAMAGE_PROBLEM
%     DATA_regress = MAW_ECM_regressionDMG(qLATENT,wALL,DATA_interp,qLATENT_all(:,LOCALdata.ind_nonlinear) ) ;
% else
DATA_regress = MAW_ECM_regression(qLATENT,wALL,DATA_interp ) ;

if IS_DAMAGE_PROBLEM
    DATA_regress.IndexLinear_damageMODEL = 1;
    DATA_regress.IndexNonlinear_damageMODEL = 2;
end


%end


ECMdata.wRED.DATA_regress = DATA_regress;
ECMdata.setPoints = setPointsALL ;
% ECMdata.wRED.Values = wALL ;
% ECMdata.wRED.q = qLATENT ;
ECMdata.setElements = setElements_cand ;


DATAoffline = DefaultField(DATAoffline,'Index_q_used_for_regression_WEIGHTS_SAW_ECM',DATA_interp.IndexLatentPlasticVariable) ;


ECMdata.wRED.IndexDOFl_q = DATAoffline.Index_q_used_for_regression_WEIGHTS_SAW_ECM;
disp(['Selected elements = ',num2str(setElements_cand(:)')])


for   icmp = 1:size(A{1},2)
    Acomp = zeros(length(setPointsALL),length(A)) ;
    for icluster = 1:length(A)
        Acomp(:,icluster) = A{icluster}(setPointsALL,icmp) ;
    end
    
    
    figure(3043+icmp)
    hold on
    title(['Evolution internal forces MAW-ECM points, comp ',num2str(icmp)])
    xlabel('qPLAST')
    ylabel('Internal force density')
    for ipoints = 1:length(setPointsALL)
        plot(qLATENT,Acomp(ipoints,:),'DisplayName',['Point = ',num2str(setPointsALL(ipoints)),' elem = ',num2str(setElements_cand(ipoints))])
    end
    legend show
    
    
    
end

disp('Cheching accuracy MAW-ECM')
CheckAccurayMAT_ECMl2(ECMdata,A,Wfe,qLATENT,wALL) ;


%if LOCAL.MAKE_GIF
make_weights_evolution_gif(historyW, qLATENT, VOL, LOCAL, historyTitle);
%end


%
%
% function hhh = SAWECMlocalPLOT_1(LOCAL,hhh,wADAPT,VOL)
%
%
% if ~LOCAL.MAKE_GIF
%     disp('PLotting new distribution')
%     % Suppose you already have hhh from your plotting loop
%     for iii = 1:length(hhh)
%         % Update the Y data for each line
%
%         if all(wADAPT(iii,:) == 0)
%             set(hhh(iii), 'Visible', 'off', 'HandleVisibility', 'off');   % hide curve
%         else
%             set(hhh(iii), 'YData', wADAPT(iii,:)/VOL*100);
%         end
%
%         % Optional: update the label with the current max
%         if iii ==1
%             [mmm,iMMM]= max(wADAPT(iii,:)) ;
%             maxval = mmm/VOL*100;
%         else
%             mmm = wADAPT(iii,iMMM) ;
%             maxval = mmm/VOL*100;
%         end
%         if maxval > 0
%             set(hhh(iii), 'DisplayName', ...
%                 ['w_{',num2str(iii),'}',' max = ', num2str(maxval, '%.3f')]);
%             %     else
%             %          set(hhh(iii), 'DisplayName', ...
%             %         ['w_{',num2str(iii),'}, eliminated'] );
%         end
%     end
%
%     % Refresh the legend if needed
%     legend show
%     set(htitle,'String',sprintf('Weights versus latent variable  — npoints %d',length(iCAND_new)));
% end
%
% end
%
%
% function  [ hhh,htitle]= PlotAUX_mawecm(LOCAL,qLATENT,wADAPT,VOL)
%
%
% if ~LOCAL.MAKE_GIF
%     figure(153)
%     hold on
%     xlabel('qLATENT')
%     ylabel('weights (non-zero)/VOL*100')
%     htitle = title(['Weights versus latent variable (initial)'])
%     hhh = zeros(length(wINI),1) ;
%     for iii  =1:length(hhh)
%         hhh(iii) = plot(qLATENT,wADAPT(iii,:),'DisplayName',['w_{',num2str(iii),'}',' max/VOL*100 = ',num2str(wADAPT(iii,1)/VOL*100)])
%     end
% else
%     hhh= [] ;
%     htitle = [] ;
% end
%
% end
