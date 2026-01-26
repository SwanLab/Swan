function  ECMdata = CECM_based_ManifAdWeights(A,DATA,Wfe,DATA_interp,zINI,wINI,...
    qLATENT,LOCAL)
%======================================================================
% CECM_based_ManifAdWeights
%
% Manifold-adaptive weights Empirical Cubature Method (MAW-ECM)
% inspired by the Continuous Empirical Cubature Method (CECM).
%
% SYNOPSIS
%   ECMdata = CECM_based_ManifAdWeights(A, DATA, Wfe, DATA_interp, ...
%                                       zINI, wINI, qLATENT, LOCAL)
%
% PURPOSE
%   Build a *manifold-adaptive* cubature rule by sparsifying an initial
%   ECM rule (zINI, wINI) across a set of parameter/latent states qLATENT.
%   Starting from *constant* weights over the manifold (standard ECM),
%   the algorithm *eliminates integration points one at a time* while
%   rebalancing the remaining weights per latent state so that a chosen
%   set of “exactness modes” remains exactly integrated.
%
%   Conceptually, this is a *structured sparsification* problem on the
%   weight trajectories w(q) over the latent manifold, with the constant
%   function included in the exactness set to preserve total volume.
%
% ORIGIN & CONTEXT
%   Methodology builds on the author’s CECM idea (CMAME) and adapts it to
%   learn *q-dependent* weights (“continuous” in the manifold sense) via a
%   greedy elimination + least-norm redistribution procedure.
%
% INPUTS
%   A            : cell array of size [nClusters x 1] (or [1 x nClusters])
%                  Each A{icluster} is an (nIP x nModes_icluster) matrix
%                  of integrand samples at the FE quadrature points for
%                  latent state q = qLATENT(icluster).
%                  NOTE: The constant function is *prepended internally*
%                  (vector of ones) to enforce exactness of total volume.
%
%   DATA         : struct with FE/mesh metadata (only fields needed here
%                  are used downstream to map points to elements).
%                  Required: DATA.MESH.ngaus (for element mapping).
%
%   Wfe          : (nIP x 1) vector of *reference* FE weights used to
%                  compute exact inner products (acts as mass measure).
%
%   DATA_interp  : struct with options for post-regression mapping from
%                  qLATENT to weights (used by MAW_ECM_regression).
%
%   zINI         : (nIP0 x 1) indices of the *initial* candidate point set
%                  (subset of FE quadrature points selected by ECM).
%
%   wINI         : (nIP0 x 1) initial *constant* weights (ECM solution).
%                  These populate all columns of wADAPT at iteration 0.
%
%   qLATENT      : (1 x nClusters) vector of latent/parameter samples
%                  (sorted internally; exact zero entries are removed to
%                  avoid duplicated/ill-posed cluster).
%
%   LOCAL        : (optional) struct with visualization/export options
%                  .MAKE_GIF     (logical) default = true
%                  .GIF_FILE     (char)    default = 'weights_evolution.gif'
%                  .DELAY_TIME   (scalar)  seconds per frame, default = 1
%                  .DPI          (scalar)  export resolution, default = 120
%                  .SHOW_LEGEND  (logical) show legend in frames, default = true
%
% OUTPUTS
%   ECMdata                 : struct with the learned manifold rule
%     .setPoints            : indices of *surviving* integration points
%     .setElements          : mapped element indices (via large2small…)
%     .wRED.DATA_regress    : regression model mapping qLATENT -> weights
%     .wRED.IndexDOFl_q     : index convention for the regressor (=1)
%
% ALGORITHM (high level)
%   0) (Preparation) For each cluster icluster (i.e., each qLATENT sample):
%        • Build an orthonormal “exactness basis” U{icluster} by SVDT of
%          [ones, A{icluster}], evaluated on the *current* point set.
%        • Form the exact right-hand side b{icluster} = U{icluster}'*Wfe.
%        • Initialize wADAPT(:,icluster) = wINI (constant across manifold).
%      Define VOL = sum(wINI) (used for normalized plotting).
%
%   1) (Greedy elimination) While the number of candidates > max #modes:
%        • Compute pointwise total weight RRR = sum over clusters of w.
%        • Sort currently positive-weight points by ascending RRR and try
%          to eliminate the smallest one (ties broken by order).
%        • Tentatively remove point p and for each cluster solve:
%              find Δw on remaining points s.t.
%                 U{icluster}' * Δw = b{icluster} - U{icluster}' * w_before
%              using *least-norm* solution   Δw = lsqminnorm(U', RHS).
%              (This preserves exactness on the basis U while minimizing
%               the ℓ2 change of weights.)
%          Reject the elimination if any updated weight becomes negative.
%        • If accepted, commit the removal (set eliminated row to zero).
%      Stopping criterion: cannot remove any further point or the number
%      of points reaches the maximum #modes required across clusters.
%
%   2) (Regression over q) After sparsification, collect the surviving
%      weights wALL over the manifold and build a regressor
%      DATA_regress = MAW_ECM_regression(qLATENT, wALL, DATA_interp),
%      which will predict weights at unseen q.
%
%   3) (Accuracy check) Call CheckAccurayMAT_ECMl to assess the final rule.
%
% MATHEMATICAL NOTES
%   • Exactness set per cluster is span{ 1, A{icluster} }, orthonormalized
%     by SVDT. The constant function ensures conservation of VOL.
%   • The redistribution subproblem is posed on U' (tall) with least-norm
%     update Δw = argmin ||Δw||_2  s.t.  U'Δw = RHS,   solved by lsqminnorm.
%   • Nonnegativity is *enforced by rejection*: if any updated weight is
%     negative, the candidate removal is discarded and the next candidate
%     is tried. (Active-set variants could be plugged in if desired.)
%
% COMPLEXITY (rough guideline)
%   Let nP be current #points, r the #exactness modes per cluster,
%   and nC the #clusters. Each tentative removal solves nC least-norm
%   linear systems of size (r x (nP-1)); overall cost is dominated by the
%   number of accepted/attempted removals times these solves.
%
% NUMERICAL/IMPLEMENTATION DETAILS
%   • qLATENT is sorted; clusters with |q| ≈ 0 are removed to avoid
%     duplication (tolerance: 1e-6*max|q|).
%   • SVDT is assumed to return an orthonormal basis (left singular vecs).
%   • The constant column (ONES) is included in every cluster basis.
%   • USE_LEAST_NORM = 1 uses lsqminnorm for robustness on tall/ill-cond.
%   • VOL = sum(wINI) is used only for normalization in plots/GIF.
%   • History of weights per iteration is stored and optionally rendered
%     to a GIF by make_weights_evolution_gif.
%
% GIF EXPORT (play once vs. loop & “static at end”)
%   • To *play once and stop*, set (in make_weights_evolution_gif):
%         imwrite(imind, cm, LOCAL.GIF_FILE, 'gif', 'LoopCount', 1, ...
%                 'DelayTime', LOCAL.DELAY_TIME);
%     IMPORTANT: Some viewers (e.g., certain web apps or OS viewers) may
%     ignore the GIF loop flag and *restart anyway*. That’s a player issue,
%     not the file. Workarounds:
%       – Export as video (e.g., .mp4) for guaranteed “stop at end”.
%       – In PowerPoint/Keynote, set animation to “Play once”.
%   • LOCAL.SHOW_LEGEND controls legend visibility in frames.
%
% LIMITATIONS / TIPS
%   • Greedy elimination with rejection can stall if many candidates lead
%     to negative weights; consider warm-starts or active-set NNLS updates
%     if you need deeper sparsity.
%   • The choice/size of the exactness basis per cluster determines the
%     minimal feasible #points; increasing modes raises the lower bound.
%   • Pre-scaling A columns (energy-based) can improve conditioning.
%
% DEPENDENCIES (called by this function)
%   SVDT, MAW_ECM_regression, CheckAccurayMAT_ECMl,
%   large2smallINCLUDEREP (points→elements mapping).
%
% EXAMPLE
%   % Build A, Wfe, zINI, wINI, qLATENT, DATA, DATA_interp beforehand…
%   LOCAL = struct('MAKE_GIF',true,'GIF_FILE','maw_evol.gif', ...
%                  'DELAY_TIME',0.8,'DPI',120,'SHOW_LEGEND',true);
%   ECMdata = CECM_based_ManifAdWeights(A, DATA, Wfe, DATA_interp, ...
%                                       zINI, wINI, qLATENT, LOCAL);
%
% SEE ALSO
%   ECM (standard), CECM (continuous), lsqminnorm, SVD/POD,
%   MAW_ECM_regression, CheckAccurayMAT_ECMl.
%
% REFERENCE
%   J. A. Hernández, “Continuous Empirical Cubature Method,”
%   Computer Methods in Applied Mechanics and Engineering, CMAME, (year).
%   (This implementation is an author’s extension for manifold-adaptive
%    weights and greedy sparsification.)
%
% AUTHOR / VERSION
%   Joaquín A. Hernández (JAHO), UPC, Terrassa.
%   v1.0 — 17-Sep-2025: initial public header & GIF history export.
%======================================================================

if nargin == 0
    load('tmp2.mat')
    close all
    %   DATA_interp.PortionExtrapolation_plot_WEIGHTS = 0.3;
    %  DATA_interp.Extrapolation_Method_WEIGHTS_projection_active_set =  1 ;
end

if nargin < 8 || isempty(LOCAL)
    LOCAL.MAKE_GIF   = true;                 % if true, build GIF after loop
    LOCAL.GIF_FILE   = 'weights_evolution.gif';
    LOCAL.DELAY_TIME = 1;                  % seconds per frame
    LOCAL.DPI        = 120;                   % export resolution
    LOCAL.SHOW_LEGEND = true;                 % legend in the GIF frames
end


% Derivation of cubature rules with Manifold-adaptive weights
% JAHO, 17th-September-2025, Wednesday, Terrassa UPC
% Here A is a cell array containing the integrand functions for each value
% of qLATENT. zINI and wINI are the initial cubature rule, while


% STEP 1
% Determine the basis matrices of each cluster using A = A{1},A{2} ...

[qLATENT,ii] = sort(qLATENT) ;
[aab,bbb] = find(abs(qLATENT)<1e-6*max(abs(qLATENT))) ;
if ~isempty(bbb)
    qLATENT(bbb) = [] ;
    ii(bbb) = [] ;
end
A = A(:,ii) ;



%sqW = sqrt(Wfe) ;
%WfeD = diag(sparse(Wfe)) ;
%Q = bsxfun(@times,Qw,1./sqW) ;  % Now Q'*WfeD*Q = I
U = cell(size(A));
b = cell(size(A,1),1) ;
% Basis matrix for
ONES = ones(size(A{1},1),1) ; % This is a vector of ones. It represents the constant function
wADAPT = zeros(length(zINI),length(qLATENT)) ;
numberMODES = zeros(length(qLATENT),1) ;
for icluster = 1:length(A)
    [U{icluster},SSSS]  = SVDT([ONES,A{icluster}]) ;
    
    b{icluster}  = U{icluster}'*Wfe ;   % This is the EXACT integral of this vector
    
    U{icluster}  = U{icluster}(zINI,:) ;
    %   b{icluster}  = U{icluster}'*wINI ; % This is the ECM integral of this vector
    
    
    wADAPT(:,icluster) = wINI ;
    numberMODES(icluster) = length(b{icluster}) ;
end

%
VOL  = sum(wINI) ;

if ~LOCAL.MAKE_GIF
    
    figure(153)
    hold on
    xlabel('qLATENT')
    ylabel('weights (non-zero)/VOL*100')
    htitle = title(['Weights versus latent variable (initial)'])
    hhh = zeros(length(wINI),1) ;
    for iii  =1:length(hhh)
        hhh(iii) = plot(qLATENT,wADAPT(iii,:),'DisplayName',['w_{',num2str(iii),'}',' max/VOL*100 = ',num2str(wADAPT(iii,1)/VOL*100)])
    end
    
end
% ---- NEW: history storage (no plotting in-loop) -------------------------
historyW = {};              % cell array; each entry is wADAPT at an iteration
historyTitle = {};          % optional: titles per iteration
% Capture the initial state before reduction starts (as "iteration 0")
historyW{end+1} = wADAPT;
historyTitle{end+1} = sprintf('Weights vs qLATENT — initial (iter 0)');


% Initialization
iCAND = 1:length(wINI);  % Candidate points to be eliminated (local indexes)
jELIM = [] ;
USE_LEAST_NORM = 1;
k  =1;
max_nmodes= max(numberMODES) ; % Maximum number of modes
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
    
    kLOC = 1;
    while kLOC <= length(ilocPOS)
        indMIN = indSORT(kLOC) ;
        indMINglo = ilocPOS(indMIN) ;
        iCAND_new = setdiff(iCAND,indMINglo) ; % = [] ;
        
        %
        %     [~,indMIN] = min(RRR(ilocPOS)) ;
        %     indMINglo = ilocPOS(indMIN) ;
        %
        %
        %     % Now we loop over all clusters
        [EXIT_WHILE_cluster,wADAPT_new] = LoopCheckClusterMAWecm(wADAPT,U,iCAND_new,indMINglo,b,USE_LEAST_NORM) ; 
       
        
        
        if EXIT_WHILE_cluster == 1
            kLOC = kLOC +1 ;
        else
            iCAND = iCAND_new ;
            jELIM = [jELIM; indMINglo] ;
            wADAPT_new(jELIM,:) = 0  ;
            wADAPT = wADAPT_new ;
            break
        end
        
    end
    
    
    if kLOC > length(ilocPOS)
        disp('The algorithm cannot reduce any more the number of points...Exiting')
        break
    end
    
    
    % ---- NEW: store snapshot for this iteration (no plotting here) ------
    historyW{end+1} = wADAPT; %#ok<AGROW>
    historyTitle{end+1} = sprintf('Weights vs qLATENT — npoints %d',length(iCAND_new));
    
    
    if ~LOCAL.MAKE_GIF
        
        disp('PLotting new distribution')
        
        % Suppose you already have hhh from your plotting loop
        for iii = 1:length(hhh)
            % Update the Y data for each line
            
            if all(wADAPT(iii,:) == 0)
                set(hhh(iii), 'Visible', 'off', 'HandleVisibility', 'off');   % hide curve
            else
                set(hhh(iii), 'YData', wADAPT(iii,:)/VOL*100);
            end
            
            % Optional: update the label with the current max
            if iii ==1
                [mmm,iMMM]= max(wADAPT(iii,:)) ;
                maxval = mmm/VOL*100;
            else
                mmm = wADAPT(iii,iMMM) ;
                maxval = mmm/VOL*100;
            end
            if maxval > 0
                set(hhh(iii), 'DisplayName', ...
                    ['w_{',num2str(iii),'}',' max = ', num2str(maxval, '%.3f')]);
                %     else
                %          set(hhh(iii), 'DisplayName', ...
                %         ['w_{',num2str(iii),'}, eliminated'] );
            end
        end
        
        % Refresh the legend if needed
        legend show
        set(htitle,'String',sprintf('Weights versus latent variable  — npoints %d',length(iCAND_new)));
    end
    k = k+1;
    
end
RRR= sum(wADAPT,2) ;
ilocPOS = find(RRR>0) ;

disp(['Final number  of points = ',num2str(length(ilocPOS))])
wALL = wADAPT(ilocPOS,:) ;
setPointsALL = zINI(ilocPOS) ;


% CONSTRUCT THE REGRESSION MAPPING FROM qLATENT to wALL
DATA_regress = MAW_ECM_regression(qLATENT,wALL,DATA_interp ) ;


ECMdata.wRED.DATA_regress = DATA_regress;
ECMdata.setPoints = setPointsALL ;
% ECMdata.wRED.Values = wALL ;
% ECMdata.wRED.q = qLATENT ;
setElements_cand = large2smallINCLUDEREP(setPointsALL,DATA.MESH.ngaus) ;
ECMdata.setElements = setElements_cand ;
ECMdata.wRED.IndexDOFl_q = 1;
disp(['Selected elements = ',num2str(setElements_cand(:)')])

disp('Cheching accuracy MAW-ECM')
CheckAccurayMAT_ECMl(ECMdata,A,Wfe,qLATENT,wALL) ;


if LOCAL.MAKE_GIF
    make_weights_evolution_gif(historyW, qLATENT, VOL, LOCAL, historyTitle);
end

end

function make_weights_evolution_gif(historyW, qLATENT, VOL, LOCAL, historyTitle)
    % Create a GIF that shows the evolution of weights per integration point
    % historyW: cell, each cell is [nPoints x nClusters] weights at an iteration
    % qLATENT: 1 x nClusters
    % VOL: scalar to normalize (% of volume)
    % LOCAL: settings (GIF_FILE, DELAY_TIME, DPI, SHOW_LEGEND)
    % historyTitle: cell of strings (titles per frame)
    
    [nPoints, nClusters] = size(historyW{1});
    
    % Prepare figure
    f = figure('Color','w','Visible','off');
    ax = axes('Parent',f); hold(ax,'on'); box(ax,'on');
    xlabel(ax,'qLATENT'); ylabel(ax,'weights (non-zero) / VOL * 100');
    
    % Pre-create line handles for consistency across frames
    hLines = gobjects(nPoints,1);
    for iii = 1:nPoints
        hLines(iii) = plot(ax, qLATENT, nan(1,nClusters), 'DisplayName', sprintf('w_{%d}',iii));
    end
    if LOCAL.SHOW_LEGEND
        legend(ax,'show','Location','eastoutside');
    end
    
    % axis limits based on all frames to avoid resizing
    maxVal = 0;
    for it = 1:numel(historyW)
        maxVal = max(maxVal, max(historyW{it}(:)));
    end
    ylim(ax, [0, max(1e-12, maxVal/VOL*100)*1.05]);
    
    % Write GIF
    for it = 1:numel(historyW)
        W = historyW{it};
        for iii = 1:nPoints
            y = W(iii,:)/VOL*100;
            if all(W(iii,:)==0)
                set(hLines(iii), 'YData', nan(1,nClusters), 'Visible','off', 'HandleVisibility','off');
            else
                set(hLines(iii), 'YData', y, 'Visible','on', 'HandleVisibility','on', ...
                    'DisplayName', sprintf('w_{%d}', iii));
            end
        end
        if LOCAL.SHOW_LEGEND
            legend(ax,'show');
        end
        if it <= numel(historyTitle)
            title(ax, historyTitle{it});
        else
            title(ax, sprintf('Weights vs qLATENT — frame %d', it));
        end
        
        drawnow;
        
        % Capture frame and write/append to GIF
        f.Position(3:4) = [1000 600]; % make it reasonably wide
        frame = getframe(f);
        [imind, cm] = rgb2ind(frame2im(frame), 256);
        if it == 1
%             imwrite(imind, cm, LOCAL.GIF_FILE, 'gif', ...
%                 'LoopCount', Inf, 'DelayTime', LOCAL.DELAY_TIME);
       %     imwrite(imind, cm, LOCAL.GIF_FILE, 'gif', ...
   % 'LoopCount', 1, 'DelayTime', LOCAL.DELAY_TIME);
            
            imwrite(imind, cm, LOCAL.GIF_FILE, 'gif', 'DelayTime', LOCAL.DELAY_TIME);
        else
            imwrite(imind, cm, LOCAL.GIF_FILE, 'gif', ...
                'WriteMode', 'append', 'DelayTime', LOCAL.DELAY_TIME);
        end
    end
    close(f);
    fprintf('GIF written to: %s\n', LOCAL.GIF_FILE);
end

