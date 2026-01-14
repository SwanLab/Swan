function [z,w,ERROR_GLO,DATAOUT]= EmpiricalCubatureMethod_CANDcomplV(BasisF,SingVal_F,W,DATA,wOLD,zOLD)
% Modification of EmpiricalCubatureMethod_CANDcompl, described below 
% The goal of the modification is to explore whether different choices of
% candidate smooths out the transition in multicluster problems
% 6-Sept-2025, Balmes 185, Barcelona 
% 
%==========================================================================
% EmpiricalCubatureMethod_CANDcompl
%--------------------------------------------------------------------------
% Purpose
%   Discrete Empirical Cubature Method (ECM) point/weight selection with a
%   "complementary candidates" safeguard. If the active candidate set is
%   not rich enough to reach the target tolerance, the algorithm
%   automatically enlarges the pool by adding the complementary set of
%   points and resumes the greedy search. This makes the procedure
%   (theoretically) fail-safe with respect to premature stagnation.
%
% Context
%   This routine implements the ECM greedy selection used in the 2024 work
%   on Subspace Adaptive Weights ECM (Hernández, Bravo, Ares de Parga,
%   Rossi, 2024), preserving the classic OMP-like structure (select, update
%   LS weights, prune nonpositive weights).
%
% Inputs
%   BasisF    : (n_f x M) basis/evaluations of the integrand snapshots or
%               reduced basis functions at the M quadrature points. We use
%               its transpose as G = BasisF'.
%   SingVal_F : (optional) singular values associated with BasisF rows (for
%               error measurement/scaling). If empty or inconsistent,
%               defaults to ones.
%   W         : (M x 1) underlying (full) positive quadrature weights
%               (e.g., Gauss weights) at each point.
%   DATA      : struct with options (all optional). Relevant fields:
%       .IncludeSingularValuesF (0)   -> If 1, scales columns of G and b
%                                        by sqrt(SingVal_F). *Not
%                                        recommended* (kept for backward
%                                        compatibility; warns if used).
%       .TOL (0)                      -> Target relative error tolerance.
%       .TOLFilterCandidatePoints(1e-6)-> Filters very small columns of G.
%       .RemoveColumnsWithNegativeProjection (0) -> (legacy; unused here).
%       .IND_POINTS_CANDIDATES ([])   -> Restrict initial candidate set.
%       .USE_SINGULAR_VALUES_INTF_FOR_MEASURING_ECM_ERROR (1)
%                                     -> Measure error using SingVal_F.
%       .npoints (length(b))          -> Max number of selected points.
%
% Outputs
%   z         : indices of selected (reduced) integration points (subset of
%               {1,...,M}).
%   w         : reduced positive weights associated with z, scaled to the
%               original measure (alpha .* sqrt(W(z))).
%   ERROR_GLO : history of the “actual” relative error at each iteration
%               (with or without singular-value weighting depending on
%               DATA.USE_SINGULAR_VALUES_INTF_FOR_MEASURING_ECM_ERROR).
%   DATAOUT   : structure with diagnostics; currently:
%       .kiteraciones -> number of iterations performed.
%
% Algorithm (high level)
%   1) Build G = BasisF' and the “exact integral” vector b = G * sqrt(W).
%   2) Initialize candidate set y (optionally filtered/pruned).
%   3) Greedy loop:
%        - Pick new index i that maximizes normalized correlation
%          (ObjFun = (G(:,y)'*r) ./ ||G(:,y)||).
%        - Update LS weights alpha and inverse normal matrix H efficiently.
%        - Enforce positivity: if any weights <= 0, move those points back
%          to the candidate set and downdate H.
%        - Compute residual r and error; store histories.
%        - If stuck (no progress) or out of candidates, enlarge y by
%          appending the complementary set.
%      Stop when error <= TOL, or when reaching limits on iterations/points.
%
% Notes / Tips
%   * The positivity enforcement ensures a positive-weight cubature rule.
%   * The complementary-set enlargement is the main difference w.r.t. the
%     “orig” version; it prevents the search from stalling due to an
%     impoverished candidate pool.
%   * Set DATA.npoints to cap sparsity (upper bound on |z|).
%   * If DATA.IncludeSingularValuesF=1 is set, the routine warns because
%     that path proved unreliable in practice; it is preserved only for
%     reproducibility of historical runs.
%
% References
%   * J. A. Hernández, J. R. Bravo, S. Ares de Parga, R. Rossi (2024).
%     Subspace Adaptive Weights ECM (SAW-ECM), CMAME.
%   * Classic ECM / OMP-style greedy cubature literature.
%==========================================================================

% --- Developer test hook: load workspace if called with no args ----------
if nargin == 0
    load('tmp.mat')
end

% --- Robust default for singular values ----------------------------------
if isempty(SingVal_F) || length(SingVal_F) ~= size(BasisF,1)
    SingVal_F  = ones(size(BasisF,2),1) ;
end

% --- Option: legacy singular-value scaling of G and b (not recommended) ---
DATA = DefaultField(DATA,'IncludeSingularValuesF',0) ;
if DATA.IncludeSingularValuesF == 1
    warning('This option has proved unreliable...disable it')
    % Historical variant: scale columns of G by sqrt(SingVal_F)
    % NOTE: Earlier versions used times(SingVal_F) (pre-2019)
    G = bsxfun(@times,BasisF',sqrt(SingVal_F));   % G is (M x n_f)
    b = G*sqrt(W) ;                               % “exact” integral vector
    bEXACT = b ;
else
    % Standard path: no scaling of G/b (recommended)
    G = BasisF' ;
    clear BasisF                                  % free memory
    b = G*sqrt(W) ;
    % For “actual error” measurement we may weight by SingVal_F (optional)
    if ~isempty(SingVal_F)
        bEXACT = b.*SingVal_F ;
    else
        bEXACT = b ;
    end
end
nbEXACT = norm(bEXACT) ;                          % norm used in error metric

% --- Basic sizes and tolerances ------------------------------------------
Gnorm = sqrt(sum(G.*G,1)) ;   % column 2-norms of G (1 x M)
M     = size(G,2) ;           % number of candidate points
DATA  = DefaultField(DATA,'TOL',0) ;
TOL   = DATA.TOL ;

% --- INITIALIZATION of sets/state ----------------------------------------
z = [] ;                       % selected indices (support)
y = (1:M) ;                    % candidate pool (start: all points)

% Optionally filter columns with tiny norm to avoid numerical issues
DATA = DefaultField(DATA,'TOLFilterCandidatePoints',1e-6) ;

% IMPORTANT: In early implementations, G included an extra row (sqrt(W));
% the current code uses the true G only (fix from 18-DEC-2022).
GnormNOONE = sqrt(sum(G.*G,1)) ;

if DATA.TOLFilterCandidatePoints > 0
    TOL_REMOVE = DATA.TOLFilterCandidatePoints*norm(b) ;
    rmvpin     = find(GnormNOONE(y) < TOL_REMOVE) ;
    y(rmvpin)  = [] ;
end

% (Legacy flag kept for API compatibility; not used below)
DATA = DefaultField(DATA,'RemoveColumnsWithNegativeProjection',0);

% If a user-provided subset of candidates is given, intersect with y
DATA = DefaultField(DATA,'IND_POINTS_CANDIDATES',[]) ;
if ~isempty(DATA.IND_POINTS_CANDIDATES)
    y = intersect(y,DATA.IND_POINTS_CANDIDATES) ;
end

% Change 6-Sept-2025
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/09_SAW_ECMplast.mlx
% We try to enforce the same set of weights as in the previous Cluster 

if ~isempty(zOLD)
    yCAND_2 = setdiff(y,zOLD) ;
    if ~isempty(yCAND_2)
        yGLO = {zOLD,yCAND_2} ; 
    else
        yGLO = {y}  ; 
    end
else
    yGLO = {y} ; 
end



yORIG = y ;                     % remember original candidate set

% Error measurement option (with or without singular values)
DATA  = DefaultField(DATA,'USE_SINGULAR_VALUES_INTF_FOR_MEASURING_ECM_ERROR',1) ;
USEsingvERR = DATA.USE_SINGULAR_VALUES_INTF_FOR_MEASURING_ECM_ERROR ;

% Greedy loop state
alpha   = [] ;                  % current LS weights for selected columns
mPOS    = 0 ;                   % number of positive (active) weights
r       = b ;                   % residual vector (start at exact b)
k       = 1 ;                   % iteration counter
errorGLO   = [] ;               % (unused local; kept for compatibility)
ERROR_GLO  = [] ;               % store actual error per iteration
NPOINTS    = [] ;               % store |z| per iteration

% Cap on number of selected points
DATA = DefaultField(DATA,'npoints',length(b)) ;
m    = min(DATA.npoints,length(b)) ;

% Precompute norms and loop limits
normB          = norm(b) ;
nerror         = norm(r)/normB ;
nerrorACTUAL   = nerror ;
H              = [] ;           % inverse of (Gz'*Gz) (Sherman-Morrison updates)
y              = y(:) ;
NITERACIONES   = 50*m ;         % generous upper bound on iterations
iOLD           = [] ;           % last selected index (only for debugging)
yCOMPL         = (1:length(W))';
yCOMPL(y)      = [] ;           % complementary set (used if we get stuck)

% Heuristic: if |z| does not increase for several iters or y is empty, we
% enlarge the pool with yCOMPL.
NITERATIONS_NO_MORE_POINTS = 10 ;
ITER_MINIMUM               = 0 ;
MAX_NUMBER_POINTS          = 0 ;
iLY = 1; 
    y = yGLO{iLY} ; 
%============================== MAIN LOOP =================================
% Stop if: (a) error <= TOL, (b) |z| reached m, or (c) iteration cap.
while  (nerrorACTUAL > TOL) && (mPOS < m) && (k < NITERACIONES)
    
    %-------------------- STEP 1: greedy selection ------------------------
    % Correlate residual with all candidate columns and normalize by column
    % norm to avoid bias towards large-norm columns (akin to OMP). 
    
   % while  iLY <= length(y)
    
    ObjFun = G(:,y)' * r ;             % correlations
    ObjFun = ObjFun ./ Gnorm(y)';      % normalized scores
    %[maxLOC, indSORT] = max(ObjFun) ;  %#ok<ASGLU> pick best index
    [ObjFunSorted,INDEXESall]  = sort(ObjFun,'descend') ; 
    indSELECT=  min(2,length(ObjFunSorted)); 
    indSORT = INDEXESall(indSELECT) ; 
    i = y(indSORT) ;
    
    %-------------------- STEP 2: LS update (unrestricted) ----------------
    % Update alpha and the inverse normal matrix H efficiently.
    if k == 1
        % First point: closed-form LS
        alpha = G(:,i) \ b ;
        H     = 1 / (G(:,i)'*G(:,i)) ;
    else
        % Rank-1 update using previously selected set and inverse (fast)
        [H, alpha] = UpdateWeightsInverse(G(:,z), H, G(:,i), alpha, r) ;
    end
    
    %-------------------- STEP 3: move i from candidates to selected -------
    z = [z; i] ;
    y(indSORT(1)) = [] ;
    
    %-------------------- STEP 4–5: positivity enforcement ----------------
    % If any weights became nonpositive, drop those points back to y and
    % downdate H accordingly, then recompute the LS solution.
    n = find(alpha <= 0) ;
    if ~isempty(n)
        y = [y; z(n)];                     % return offending points
        z(n) = [] ;                         % remove from support
        H = MultiUpdateInverseHermitian(H,n); % consistent downdate
        alpha = H * (G(:,z)' * b);         % recompute LS on the pruned set
    end
    
    
    %-------------------- STEP 6: residual and error ----------------------
    r      = b - G(:,z)*alpha ;
    nerror = norm(r)/norm(b) ;            % nominal relative error
    
    % “Actual” error optionally includes singular-value weighting
    if DATA.IncludeSingularValuesF == 0 && USEsingvERR == 1
        nerrorACTUAL = SingVal_F .* r ;
        nerrorACTUAL = norm(nerrorACTUAL / nbEXACT);
    else
        nerrorACTUAL = nerror ;
    end
    
    %-------------------- Stagnation handling (pool enrichment) -----------
    % If we are not adding new points for too long or we ran out of y,
    % enlarge y by appending the complementary set and continue.
    if ITER_MINIMUM > NITERATIONS_NO_MORE_POINTS || (isempty(y) && ( length(z) < m  ) )
        disp('--------------------------------------------------------------')
      
        if iLY < length(yGLO)
              disp('The algorithm cannot proceed with the current set of points(the same as the previous cluster) ')
        disp('We enlarge the set of candidates with the first complementary set')
            iLY = iLY  +1 ; 
                y = yGLO{iLY} ; 
        else
             y = [y; yCOMPL] ;
              disp('We enlarge the set of candidates with the remaining points')
        end
       
        ITER_MINIMUM = 0 ;
    end
    
    
    % Track whether we just increased |z|
    if length(z) > MAX_NUMBER_POINTS
        ITER_MINIMUM = 0 ;
    else
        ITER_MINIMUM = ITER_MINIMUM + 1 ;
    end
    iOLD = i ; %#ok<NASGU> (kept for debugging/inspection)
    

    
    %-------------------- STEP 7: logging ---------------------------------
    disp(['k = ',num2str(k),', m=',num2str(length(z)),' ,','; error n(res)/n(b) (%) = ',...
        num2str(nerror*100,10),';  Actual error % =',num2str(nerrorACTUAL*100,10)]) ;
    ERROR_GLO(k) = nerrorACTUAL ;
    NPOINTS(k)   = length(z) ;
    
    %-------------------- Book-keeping for loop control -------------------
    mPOS = length(z) ;
    MAX_NUMBER_POINTS = max(mPOS,MAX_NUMBER_POINTS) ;
    k = k + 1 ;
end
%=========================== END MAIN LOOP ================================

% If we hit the iteration cap, advise the user to enlarge the initial pool.
if k >= NITERACIONES
    PROPORpoints = length(yORIG)/length(W)*100;
    error(['NO CONVERGENCE. ENLARGE THE SET OF CANDIDATE POINTS  (NOW IT CONTAINS ',...
        num2str(PROPORpoints),' % of the total number of Gauss points)'])
end

% Map LS weights (alpha) to final reduced weights on the original measure
w = alpha .* sqrt(W(z)) ;

disp(['Total number of iterations =',num2str(k)])

% Optional quick plot (disabled by default to keep function side-effect free)
PLOT_error = 0 ;
if PLOT_error == 1
    figure(500); hold on
    xlabel('Number of points'); ylabel('Error (%)')
    plot(NPOINTS, ERROR_GLO*100, 'k')
end

% Diagnostics
DATAOUT.kiteraciones = k ;

% % (Optional) debug plot of |z| vs iterations:
% figure(501); hold on
% xlabel('Number of iterations'); ylabel('Number of points')
% plot(NPOINTS,'k')
%==========================================================================

end






% BELOW IS THE VERSION BEFORE CHATGPT COMMENTS, 6-SEPT-2025

% function [z,w,ERROR_GLO,DATAOUT]= EmpiricalCubatureMethod_CANDcompl(BasisF,SingVal_F,W,DATA)
% % Version created 23-FEB-2021 to cope with clustering
% % The only difference with respect to
% % ../SVDlibrary/ECM/EmpiricalCubatureMethod_orig.m
% % is that this version never fails to converge (in theory). If the
% % set of candidates does not contain any solution, and the algorithm fails
% % to converge in a given number of iterations, then the set of candidates is enlarged
% % with the complementary, initial set, and then the search goes on
% % ----------------------------------------------------
%
% if nargin == 0
%     load('tmp.mat')
% end
%
% if isempty(SingVal_F) || length(SingVal_F) ~= size(BasisF,1)
%     SingVal_F  =ones(size(BasisF,2),1) ;
% end
%
% DATA = DefaultField(DATA,'IncludeSingularValuesF',0) ; %  ; .IncludeSingularValuesF = 1
% if DATA.IncludeSingularValuesF == 1
%     warning('This option has proved unreliable...disable it')
%     %   G = bsxfun(@times,BasisF', (SingVal_F));  %  Version before  5th Dec-2019...
%     G = bsxfun(@times,BasisF',sqrt(SingVal_F));  %
%     b = G*sqrt(W) ;  % b Vector (exact integral)
%     bEXACT = b ;
% else
%     G = BasisF' ;
%     clear BasisF
%     b = G*sqrt(W) ;  % b Vector (exact integral)
%     if ~isempty(SingVal_F)
%         bEXACT = b.*SingVal_F ;
%     else
%         bEXACT =b ;
%     end
% end
% nbEXACT = norm(bEXACT) ;
%
%
% Gnorm =sqrt(sum(G.*G,1)) ; % Norm of Columns of G
% M = size(G,2) ;  % Number of FE points
% DATA = DefaultField(DATA,'TOL',0) ; % Default tolerance for convergence
% TOL = DATA.TOL ;
% % INITIALIZATIONS
% % ------------------------
% z = [] ; % Set of integration points
% % Set of candidate points (those whose associated column has low norm are removed)
%
% %  PointsWithZero =  find(sum(G(1:end-1,:),1)==0) ;
%
% y=1:M ;
% %y(PointsWithZero) = []  ;
% DATA = DefaultField(DATA,'TOLFilterCandidatePoints',1e-6) ;
% %GnormNOONE =sqrt(sum(G(1:end-1,:).*G(1:end-1,:),1)) ; % Norm of Columns of G  % change 18-DEC-2022
% % The above was a mistake originated from the first implementations in
% % which G matrix was "expanded" with a row (sqrt(W))
% GnormNOONE =sqrt(sum(G.*G,1)) ; % Norm of Columns of G  % change 18-DEC-2022
%
%
% if DATA.TOLFilterCandidatePoints >0
%     TOL_REMOVE = DATA.TOLFilterCandidatePoints*norm(b) ;
%     rmvpin = find(GnormNOONE(y)<TOL_REMOVE) ;
%     y(rmvpin) = [] ;
% end
%
% DATA = DefaultField(DATA,'RemoveColumnsWithNegativeProjection',0); % 4-Dec-2019
%
%
%
% DATA = DefaultField(DATA,'IND_POINTS_CANDIDATES',[]) ;
%
% if ~isempty(DATA.IND_POINTS_CANDIDATES)
%     y = intersect(y,DATA.IND_POINTS_CANDIDATES) ;
% end
% yORIG = y ;
%
% DATA  = DefaultField(DATA,'USE_SINGULAR_VALUES_INTF_FOR_MEASURING_ECM_ERROR',1) ;  % 27th-april-2020
% USEsingvERR = DATA.USE_SINGULAR_VALUES_INTF_FOR_MEASURING_ECM_ERROR ;
%
%
% alpha = [] ; % Vector of weights
% mPOS = 0 ; % Number of nonzero weights
% r = b ; % Residual vector
% k = 1;  % Number of iterations
% errorGLO = [] ; %(for storing error)
% % Default number of points
% DATA = DefaultField(DATA,'npoints',length(b)) ;
% m = min(DATA.npoints,length(b)) ;
% % END INITIALIZATIONS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normB = norm(b) ;
% nerror = norm(r)/normB  ;
% H = [] ; % Inverse of (Gz'*Gz)
% ERROR_GLO = [] ;
% NPOINTS =[] ;
% nerrorACTUAL = nerror;
% y = y(:) ;
% NITERACIONES = 50*m ;
% iOLD = [] ;
% yCOMPL = (1:length(W))';
% yCOMPL(y) = [] ; % Complementary set of the candidate set (in case this set is not sufficiently rich)
%
%
% NITERATIONS_NO_MORE_POINTS = 10 ; %    NUMBER OF ITERATIONS IN WHICH THE CODE IS ALLOWED TO "ITERATE" IN ORDER TO FIND
% %  THE MINIMUM
%
% ITER_MINIMUM = 0;
% MAX_NUMBER_POINTS = 0 ;
% % while  nerrorACTUAL >TOL && mPOS <m   && ~isempty(y) && k<NITERACIONES
% % ...Before Nov. 4th -2021, see EmpiricalCubatureMethod_CANDcompl_aux.mlx
% while  nerrorACTUAL >TOL && mPOS <m     && k<NITERACIONES
%
%     % TIMELOC_k = tic ;
%     % STEP 1. Compute new point
%     ObjFun = G(:,y)'*r ;
%     ObjFun = ObjFun./Gnorm(y)';
%     [maxLOC, indSORT] = max(ObjFun)  ;
%     i = y(indSORT(1)) ;
%     % STEP 2.  Update alpha and H  (unrestricted least-squares)
%     if k==1
%         alpha =  G(:,i)\b ;
%         H = 1/(G(:,i)'*G(:,i)) ;
%     else
%         [H alpha] = UpdateWeightsInverse(G(:,z),H,G(:,i),alpha,r) ;
%     end
%
%     % STEP 3. Move i from set y to set z
%     z = [z;i] ;     y(indSORT(1)) = [] ;
%     % STEP 4. Find possible negative weights
%     n = find(alpha<=0) ;
%     if  ~isempty(n)
%         % STEP 5
%         y = [y; z(n)];  z(n)=[] ;
%         H = MultiUpdateInverseHermitian(H,n) ;
%         % Recomputing alpha
%         alpha = H*(G(:,z)'*b );
%         % else
%         % ITER_MINIMUM = 0 ;
%     end
%
%     if ITER_MINIMUM > NITERATIONS_NO_MORE_POINTS || isempty(y)
%         disp('--------------------------------------------------------------')
%         disp(['The algorithm cannot proceed with the current set of points'])
%         disp(['We enlarge the set of candidates with the complementary set'])
%         y = [y;yCOMPL] ;
%         ITER_MINIMUM = 0 ;
%     end
%
%
%
%     if length(z) > MAX_NUMBER_POINTS
%         ITER_MINIMUM = 0 ;
%     else
%         ITER_MINIMUM = ITER_MINIMUM + 1;
%     end
%
%     iOLD = i ;
%     % STEP 6
%     r = b-G(:,z)*alpha ;
%     nerror = norm(r)/norm(b) ; % Relative error (using r and b)
%     if DATA.IncludeSingularValuesF == 0 && USEsingvERR == 1
%         nerrorACTUAL = SingVal_F.*r ;
%         nerrorACTUAL = norm(nerrorACTUAL/nbEXACT);
%     else
%         nerrorACTUAL = nerror ;
%     end
%     % STEP 7
%     disp(['k = ',num2str(k),', m=',num2str(length(z)),' ,','; error n(res)/n(b) (%) = ',...
%         num2str(nerror*100,10),';  Actual error % =',num2str(nerrorACTUAL*100,10)]) ;
%     ERROR_GLO(k) = nerrorACTUAL ;
%     NPOINTS(k) =  length(z) ;
%
%     mPOS = length(z) ;
%     MAX_NUMBER_POINTS = max(mPOS,MAX_NUMBER_POINTS) ;
%     k = k + 1 ;
%
%     %     if length(z) == m
%     %         dbstop('88')
%     %         disp('')
%     %     end
%
% end
%
%
% if  k>= NITERACIONES
%     PROPORpoints = length(yORIG)/length(W)*100;
%     error(['NO CONVERGENCE. ENLARGE THE SET OF CANDIDATE POINTS  (NOW IT CONTAINS ',num2str(PROPORpoints),' % of the total number of Gauss points)'])
% end
%
%
% w = alpha.*sqrt(W(z)) ;
%
% disp(['Total number of iterations =',num2str(k)])
%
% PLOT_error = 0 ;
%
% if PLOT_error == 1
%
%     figure(500)
%     hold on
%     xlabel('Number of points')
%     ylabel('Error (%)')
%     plot(NPOINTS,ERROR_GLO*100,'k')
% end
%
% DATAOUT.kiteraciones = k ;
%
% % figure(501)
% % hold on
% % xlabel('Number of iterations')
% % ylabel('Number of points')
% % plot(NPOINTS,'k')
%
%
