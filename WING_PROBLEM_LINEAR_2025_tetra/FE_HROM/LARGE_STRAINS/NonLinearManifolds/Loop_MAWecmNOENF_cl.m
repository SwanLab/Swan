function [wADAPT,iCAND,jELIM,kLOC] = Loop_MAWecmNOENF_cl(U,ilocPOS,iCAND,wADAPT,b,indSORT,jELIM)
%==========================================================================
% Loop_MAWecmNOENF_cl
%
% PURPOSE
% -------
% Single-cluster **candidate elimination trial** without explicit w ≥ 0
% enforcement (NO ENForcement = “NOENF”). Given:
%   • current candidate set iCAND and active-positive indices ilocPOS,
%   • per-cluster restricted bases U{c} and targets b{c},
%   • current weights wADAPT,
% try removing ONE point (chosen by the provided ascending order indSORT),
% and redistribute its weight across the remaining iCAND points by a
% least-change update that preserves exactness (U' w = b) **per cluster**.
% If any cluster becomes infeasible (negatives or rank deficiency), reject
% the trial and try the next candidate in indSORT.
%
% INPUTS
% ------
% U         : 1×nC cell; U{c} ∈ ℝ^(|zINI|×r_c), restricted cluster bases.
% ilocPOS   : Vector of indices (in local/candidate numbering) whose current
%             total weight is > 0 (i.e., not yet eliminated).
% iCAND     : Current set of active candidate indices (subset of 1:|zINI|).
% wADAPT    : |zINI|×nC current weights (columns = clusters).
% b         : 1×nC cell; b{c} ∈ ℝ^(r_c×1), target integrals for each cluster.
% indSORT   : Permutation of 1:length(ilocPOS) giving the **ascending** order
%             of candidates to try (e.g., by small total weight).
% jELIM     : Column vector of indices already eliminated (for bookkeeping).
%
% OUTPUTS
% -------
% wADAPT    : Updated weights if a feasible elimination is committed;
%             otherwise unchanged.
% iCAND     : Updated active set after committing the elimination.
% jELIM     : Updated list with the newly removed index appended.
% kLOC      : Local trial counter; on exit, either points to the candidate
%             that succeeded (commit happened) or equals length(ilocPOS)+1
%             if all trials failed (no elimination committed).
%
% METHOD (one outer attempt across clusters)
% -----------------------------------------
% For kLOC = 1 .. length(ilocPOS):
%   1) Pick the kLOC-th candidate (in ascending indSORT), map to global idx
%      indMINglo = ilocPOS(indSORT(kLOC)), and form the tentative set
%      iCAND_new = iCAND \ {indMINglo}.
%   2) For each cluster c:
%        Uloc = U{c}(iCAND_new,:)    % restrict rows to tentative set
%        If rank(Uloc) = r_c:
%           w_before = wADAPT(iCAND_new,c)
%           Solve **least-change** update that preserves exactness:
%               RHS     = b{c} − Ulocᵀ w_before
%               w_incre = lsqminnorm(Ulocᵀ, RHS)   % minimal ‖Δw‖₂ solution
%               w_after = w_before + w_incre
%           If any(w_after < 0): reject this candidate for this cluster.
%        Else: reject (rank-deficient).
%      If any cluster rejects → try next kLOC (continue).
%   3) If all clusters accept:
%        - Commit: iCAND ← iCAND_new; append indMINglo to jELIM
%        - Set wADAPT(jELIM,:) = 0 to mark eliminated point(s)
%        - Return updated (wADAPT, iCAND, jELIM) and break.
%
% NUMERICAL NOTES
% ---------------
% • The update uses a **minimal 2-norm increment** (lsqminnorm) to preserve
%   Uᵀ w = b after removing one row in U (i.e., one integration point).
% • No explicit nonnegativity enforcement is performed here; any negative
%   entry in w_after causes the trial to be rejected.
% • Full-column rank of Uloc is required (r_c ≤ |iCAND_new| and independent).
%
% BEHAVIOR & RETURN CONVENTION
% ----------------------------
% • If no candidate passes across all clusters, the loop exits with
%   kLOC = length(ilocPOS)+1, and (wADAPT, iCAND, jELIM) are unchanged.
% • On success, the function returns immediately after the first feasible
%   elimination is committed.
%
% COMPLEXITY
% ----------
% Up to length(ilocPOS) trials; each trial solves nC least-squares systems
% of size r_c × |iCAND_new| (via lsqminnorm on Ulocᵀ).
%
% DEPENDENCIES
% ------------
% SVDT (for rank check via singular values), lsqminnorm (MATLAB).
% J.A. Hernández Ortega (JAHO) — Barcelona — 7-NOV-2025
% (Comments updated by ChatGPT-5 Thinking)
%==========================================================================


if nargin == 0
    load('tmp.mat')
end



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
    
    EXIT_WHILE_cluster = 0 ;
    wADAPT_new = wADAPT ;
    for icluster = 1:length(U)
        %  disp( ['Cluster = ',num2str(icluster)])
        % Right-hand side   b-U(:,z)'*w_before_update
        
        
        
        
        Uloc = U{icluster}(iCAND_new,:) ;
        [~,SSS] = SVDT(Uloc) ;
        if length(SSS) == size(Uloc,2)
            
            w_before= wADAPT(iCAND_new,icluster) ;
            %w_eliminated = wADAPT(indMINglo,icluster)  ;
            RHS = b{icluster}-Uloc'*w_before ;
            %             USE_CHT = 1;
            %             if USE_CHT == 0
            %                 if USE_LEAST_NORM == 0
            %                     w_incre = Uloc'\RHS;
            %                 else
            w_incre   = lsqminnorm(Uloc',RHS) ;
            %                 end
            %                 CHECK_w = sum(w_incre) -w_eliminated ;
            w_after = w_before + w_incre ;
            %             else
            
            % [ w_after,SOLUTION_fOUND] = nn_update_active_set(Uloc, b{icluster}, w_before) ;
            %    errorLOC = Uloc'*w_after-b{icluster}   ;
            %   errorLOC = norm(errorLOC)/norm(b{icluster}) ;
            %             if errorLOC >1e-10
            %                  disp(['error_appox, rel, over 1 =',num2str(errorLOC)])
            %             end
            %   end
            
            
            
            
            % disp(['Difference betwee sum(Delta w) and eliminated weight =',num2str(CHECK_w)])
            if  any(w_after <0)
                % disp(['NEGATIVE WEIGHTS, exiting, cluster = ',num2str(icluster)])
                EXIT_WHILE_cluster = 1 ;
                break
            end
            % Updating
            wADAPT_new(iCAND_new,icluster) = w_after ;
            
        else
            disp('Uloc is not full rank, exiting')
            EXIT_WHILE_cluster = 1 ;
            break
        end
        
        
    end
    
    
    
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