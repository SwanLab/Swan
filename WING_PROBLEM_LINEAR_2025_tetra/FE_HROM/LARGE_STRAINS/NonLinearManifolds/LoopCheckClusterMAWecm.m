function [EXIT_WHILE_cluster,wADAPT_new] = LoopCheckClusterMAWecm(wADAPT,U,iCAND_new,indMINglo,b,USE_LEAST_NORM)
%==========================================================================
% LoopCheckClusterMAWecm
%
% PURPOSE
%   Attempt elimination of a candidate integration point (indMINglo) from
%   the current adaptive weights wADAPT across all clusters, by
%   redistributing its contribution to the remaining points (iCAND_new).
%   For each cluster, the redistribution is solved in least-norm sense so
%   that exactness on the cluster basis U{icluster} is preserved.
%
% INPUTS
%   wADAPT        : [nPoints × nClusters] current weight matrix.
%                   Each column corresponds to one latent cluster.
%   U             : cell{nClusters}, each entry is [nActive × r_c] matrix
%                   of basis functions (restricted to active candidate
%                   points) for cluster c.
%   iCAND_new     : vector of indices of remaining candidate points after
%                   tentative elimination of indMINglo.
%   indMINglo     : scalar, global index of the candidate point to be
%                   eliminated from the active set.
%   b             : cell{nClusters}, each entry is RHS vector
%                   b{c} = U{c}' * Wfe (exact integrals for cluster c).
%   USE_LEAST_NORM: logical flag; if true (1), use lsqminnorm for robust
%                   least-norm redistribution; if false (0), use backslash
%                   (less robust, exact if square).
%
% OUTPUTS
%   EXIT_WHILE_cluster : flag (0 or 1)
%                        0 → elimination accepted (no negatives found)
%                        1 → elimination rejected (negative weights arose)
%   wADAPT_new          : updated weight matrix after redistribution if
%                        successful; unchanged if elimination rejected.
%
% METHOD
%   For each cluster:
%     – Restrict basis U to surviving candidates.
%     – Compute residual RHS = b - U'*w_before.
%     – Solve U'^Δw = RHS with least-norm solution.
%     – Update weights w_after = w_before + Δw.
%     – If any weight is negative, abort elimination and return flag=1.
%
% NOTES
%   • The eliminated point’s weights are not directly reassigned; instead,
%     redistribution is computed to preserve exactness on span{1, A}.
%   • Enforces nonnegativity via *rejection* of point removal.
%   • Designed as a subroutine inside the greedy loop of
%     CECM_based_ManifAdWeights.
%
% AUTHOR
%   Joaquín A. Hernández (JAHO), UPC, Sept 2025.
%==========================================================================

if nargin == 0
    load('tmp.mat')
end

EXIT_WHILE_cluster = 0 ;
wADAPT_new = wADAPT ;
for icluster = 1:length(U)
    %  disp( ['Cluster = ',num2str(icluster)])
    % Right-hand side   b-U(:,z)'*w_before_update
    
    Uloc = U{icluster}(iCAND_new,:) ;
    
    [~,SSS] = SVDT(Uloc) ;
    if length(SSS) == size(Uloc,2)
        w_before= wADAPT(iCAND_new,icluster) ;
        w_eliminated = wADAPT(indMINglo,icluster)  ;
        RHS = b{icluster}-Uloc'*w_before ;
        if USE_LEAST_NORM == 0
            w_incre = Uloc'\RHS;
        else
            w_incre   = lsqminnorm(Uloc',RHS) ;
        end
        %
        %     r = Uloc' * w_incre - RHS;
        % res_norm = norm(r);                  % residual of the linear system
        % rel_res  = res_norm / max(1,norm(RHS));
        % CHECK_w = sum(w_incre) -w_eliminated ;
        w_after = w_before + w_incre ;
        % disp(['Difference betwee sum(Delta w) and eliminated weight =',num2str(CHECK_w)])
        if any(w_after <0)
            disp('NEGATIVE WEIGHTS, exiting')
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
