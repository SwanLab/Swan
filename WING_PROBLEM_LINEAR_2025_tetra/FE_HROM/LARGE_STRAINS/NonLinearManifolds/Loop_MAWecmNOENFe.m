function [wADAPT,iCAND,jELIM,FEASIBLE_SOLUTION_FOUND] = ...
    Loop_MAWecmNOENFe(U,ilocPOS,iCAND,wADAPT,b,indSORT,jELIM,MaxZerosMAW,Delta_q)
%==========================================================================
% Loop_MAWecmNOENFe
%
% PURPOSE
% -------
% Exhaustive single-point elimination **with explicit nonnegativity** for
% MAW-ECM across multiple clusters. In one pass, the routine evaluates the
% removal of each currently active candidate (row) and, for every cluster,
% recomputes a feasible weight vector w ≥ 0 that satisfies U{c}'*w = b{c}.
% Among all **feasible** removals in the pass, it commits the one that
% **minimizes the worst (maximum) number of zero weights across clusters**,
% promoting balanced sparsity.
%
% PROBLEM SETTING (per cluster c)
% -------------------------------
% Given active set iCAND, try iCAND_new = iCAND \ {p}. For cluster c:
%   Find w_after ≥ 0 such that  U{c}(iCAND_new,:)' * w_after = b{c},
% starting from w_before = wADAPT(iCAND_new,c). The solve is performed with
% an active-set NNLS with equality constraints (nn_update_active_set).
%
% INPUTS
% ------
% U          : 1×nC cell; U{c} ∈ ℝ^(|zINI|×r_c) restricted bases (rows align
%              with global candidate indices).
% ilocPOS    : Indices (in local numbering) of still-positive rows. Combined
%              with indSORT to map to the global row to test.
% iCAND      : Current active set of global candidate indices.
% wADAPT     : |zINI|×nC matrix of current weights.
% b          : 1×nC cell; b{c} ∈ ℝ^(r_c×1) target integrals.
% indSORT    : Permutation of 1:length(ilocPOS) giving the **trial order**
%              (typically ascending total weight).
% jELIM      : Column vector of already eliminated global indices.
% MaxZerosMAW: Upper bound on the accepted **max zero-count** across clusters.
% Delta_q    : Spacing of latent coordinate (used only in optional gradient
%              diagnostics; currently not active).
%
% OUTPUTS
% -------
% wADAPT                 : Updated weights after committing the selected
%                          removal (rows in jELIM are explicitly zeroed).
% iCAND                  : Updated active set (one index removed if success).
% jELIM                  : Updated eliminated list with the committed index.
% FEASIBLE_SOLUTION_FOUND: true if at least one candidate removal was feasible
%                          on **all** clusters and passed the zero-count test.
%
% ALGORITHM (single exhaustive pass)
% ----------------------------------
% For each trial kLOC in the order indSORT:
%   1) Map local position → global index: indMINglo = ilocPOS(indSORT(kLOC)).
%      Define iCAND_new = setdiff(iCAND, indMINglo).
%   2) For every cluster c:
%        Uloc = U{c}(iCAND_new,:);
%        Check full column rank of Uloc (via SVDT); if rank-deficient → reject.
%        w_before = wADAPT(iCAND_new,c).
%        [w_after, ok] = nn_update_active_set(Uloc, b{c}, w_before).
%        If ~ok → reject this trial; continue with next candidate.
%        Record w_after and its zero count for cluster c.
%   3) If all clusters feasible → store:
%        • wADAPT_new (assembled with w_after per cluster),
%        • iCAND_new, the removed index indMINglo,
%        • NumberZeroWeigth = max_c (#zeros in w_after_c).
% After trying all candidates:
%   • If no feasible trial or min(NumberZeroWeigth) > MaxZerosMAW:
%       FEASIBLE_SOLUTION_FOUND = false (no changes).
%   • Else select the feasible trial with minimal NumberZeroWeigth,
%       commit (update iCAND, append jELIM, set wADAPT = chosen wADAPT_new,
%       and zero rows at jELIM), FEASIBLE_SOLUTION_FOUND = true.
%
% KEY IDEAS
% ---------
% • Explicit feasibility: nonnegativity and equalities enforced per cluster.
% • Selection metric: minimize the **worst** cluster sparsity to keep a
%   balanced support across clusters (avoid over-sparsifying one cluster).
% • A hard cap MaxZerosMAW prevents accepting excessively sparse outcomes.
%
% NUMERICAL NOTES
% ---------------
% • Rank check uses SVDT(Uloc); requirement: size(Uloc,2) singular values.
% • nn_update_active_set should return the feasible solution closest (in a
%   least-change sense) to w_before; noisy “exact” alternatives (e.g., lsqlin)
%   are available but discouraged (see USE_exact_method in code).
% • Gradient smoothness diagnostics w.r.t. q (using Delta_q) are present as
%   commented code and can be re-enabled if needed.
%
% SIDE ARTIFACTS (stored for debugging/analysis)
% ----------------------------------------------
% • WeightsHistory      : all feasible wADAPT_new matrices.
% • CandidatesHistory   : corresponding iCAND_new sets.
% • NumberZeroWeigth    : max zero-count per feasible trial.
% • indMINglo_hst       : globally removed index per feasible trial.
%
% BEHAVIOR
% --------
% • Exactly one removal can be committed per pass. Outer loops (caller) can
%   invoke this routine repeatedly until no further feasible removal exists.
%
% COMPLEXITY
% ----------
% Up to |ilocPOS| trials; each feasible/infeasible trial attempts nC solves
% of NNLS-with-equalities over dimensions (|iCAND_new| × r_c).
%
% DEPENDENCIES
% ------------
% SVDT, nn_update_active_set  (optionally nn_update_lsqlin for testing).
%
% AUTHOR / PLACE / DATE
% ---------------------
% J.A. Hernández Ortega (JAHO) — Barcelona — 21-Sep-2025
% (Comments updated by ChatGPT-5 Thinking)
%==========================================================================

%==========================================================================

if nargin == 0
    load('tmpE1.mat')
end

%kLOC = 1;
WeightsHistory = {} ;
CandidatesHistory = {} ;
NumberZeroWeigth = []  ;
indMINglo_hst = [] ;
% MaxGradient = [] ;
% Norm2Gradient = [] ;

USE_exact_method = 0;

for kLOC =1: length(ilocPOS) % Here we test all the possible combinations

%for kLOC =2:2 % Here we test all the possible combinations
%    warning('borrar esto')
    indMIN = indSORT(kLOC) ;
    indMINglo = ilocPOS(indMIN) ;
    iCAND_new = setdiff(iCAND,indMINglo) ; % = [] ;
    
    EXIT_WHILE_cluster = 0 ;
    wADAPT_new = wADAPT ;
    NzerosCLUSTER= zeros(length(U),1) ;
    for icluster = 1:length(U)
        Uloc = U{icluster}(iCAND_new,:) ;
        [~,SSS] = SVDT(Uloc) ;
        if   length(SSS) == size(Uloc,2)
            w_before= wADAPT(iCAND_new,icluster) ;
            if USE_exact_method == 0
                [ w_after,SOLUTION_fOUND] = nn_update_active_set(Uloc, b{icluster}, w_before) ;
            else
                % Not recommended, it produces noisy solutions
                [w_after,SOLUTION_fOUND] = nn_update_lsqlin(Uloc, b{icluster}, w_before) ;
                
            end
            NzerosCLUSTER(icluster) = length(find(w_after==0)) ;
            if ~SOLUTION_fOUND
                % disp(['No feasible solution found, exiting, icluster = ',num2str(icluster)])
                EXIT_WHILE_cluster = 1 ;
                break
                
            end
            wADAPT_new(iCAND_new,icluster) = w_after ;
        else
            disp('Uloc is not full rank, exiting')
            EXIT_WHILE_cluster = 1 ;
            break
            
        end
        
        
        
    end
    
    
    if EXIT_WHILE_cluster == 0
        indMINglo_hst(end+1) = indMINglo ;
        NumberZeroWeigth(end+1) = max(NzerosCLUSTER) ;
        WeightsHistory{end+1} = wADAPT_new; %(iCAND_new,:) ;
        CandidatesHistory{end+1} = iCAND_new ;
        % Compute gradient
        %         dw = diff(wADAPT_new')';
        %         dw_dq = bsxfun(@times,dw',1./Delta_q')';
        %         max_grad = max(max(abs(dw_dq))) ;
        %         MaxGradient(end+1) = max_grad ;
        %          ndw_dq = sum(sum(dw_dq.^2)) ;
        %          Norm2Gradient(end+1) = sqrt(ndw_dq)  ;
    end
    
    
end


minNZ = min(NumberZeroWeigth) ;

if isempty(CandidatesHistory)  || minNZ > MaxZerosMAW
    FEASIBLE_SOLUTION_FOUND = false;
    disp('Feasible solution not found')
else
    FEASIBLE_SOLUTION_FOUND = true ;
    disp('Feasible solution    found')
    % Max criterion
    %  [nZZ,indZZ_max] = min(MaxGradient)  ;
    [nZZ,indZZ] = min(NumberZeroWeigth)  ;
    %  CandiLOC =  find(NumberZeroWeigth == nZZ) ;
    %  [~,iLOC] = min(MaxGradient(CandiLOC))  ;
    %  indZZ = CandiLOC(iLOC) ;
    %       [nZZ,indZZ_norm] = min(Norm2Gradient)  ;
    %     indZZ = indZZ_norm ;
    
    disp(['Number of zeros in the set of weights = ',num2str(nZZ)])
    iCAND = CandidatesHistory{indZZ} ;
    indMINglo = indMINglo_hst(indZZ) ;
    jELIM = [jELIM; indMINglo] ;
    wADAPT_new = WeightsHistory{indZZ} ;
    wADAPT_new(jELIM,:) = 0  ;
    wADAPT = wADAPT_new ;
    %   break
end