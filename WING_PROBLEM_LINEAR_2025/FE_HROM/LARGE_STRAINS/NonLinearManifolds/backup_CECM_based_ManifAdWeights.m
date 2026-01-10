function  ECMdata = backup_CECM_based_ManifAdWeights(A,DATA,Wfe,DATAoffline,zINI,wINI,...
    Qw,qLATENT)
if nargin == 0
    load('tmp3.mat')
    close all
end
% Derivation of cubature rules with Manifold-adaptive weights
% JAHO, 17th-September-2025, Wednesday, Terrassa UPC
% Here A is a cell array containing the integrand functions for each value
% of qLATENT. zINI and wINI are the initial cubature rule, while
% Qw is the "weighted" basis matrix used by the ECM to determine zINI and
% wINI

% STEP 1
% Determine the basis matrices of each cluster using A = A{1},A{2} ... and
% Qw
% First we determine the actual basis matrix Q

[qLATENT,ii] = sort(qLATENT) ;
[aab,bbb] = find(abs(qLATENT)<1e-6*max(abs(qLATENT))) ;
if ~isempty(bbb)
    qLATENT(bbb) = [] ;
    ii(bbb) = [] ;
end
A = A(:,ii) ;



sqW = sqrt(Wfe) ;
WfeD = diag(sparse(Wfe)) ;
Q = bsxfun(@times,Qw,1./sqW) ;  % Now Q'*WfeD*Q = I
U = cell(size(A));
b = cell(size(A,1),1) ;
% Basis matrix for
ONES = ones(size(Q,1),1) ; % This is a vector of ones. It represents the constant function
wADAPT = zeros(length(zINI),length(qLATENT)) ;
for icluster = 1:length(A)
    U{icluster}  = SVDT([ONES,A{icluster}]) ;
    U{icluster}  = U{icluster}(zINI,:) ;
    b{icluster}  = U{icluster}'*wINI ; % This is the ECM integral of this vector
    wADAPT(:,icluster) = wINI ;
end

%
VOL  = sum(wINI) ;

figure(153)
hold on
xlabel('qLATENT')
ylabel('weights (non-zero)/VOL*100')
htitle = title(['Weights versus latent variable (initial)'])
hhh = zeros(length(wINI),1) ;
for iii  =1:length(hhh)
    hhh(iii) = plot(qLATENT,wADAPT(iii,:),'DisplayName',['w_{',num2str(iii),'}',' max/VOL*100 = ',num2str(wADAPT(iii,1)/VOL*100)])
end
% Initialization
iCAND = 1:length(wINI);  % Candidate points to be eliminated (local indexes)
jELIM = [] ;
USE_LEAST_NORM = 1;
k  =1;
while   length(iCAND) >= 2
    
    disp(['iter = ',num2str(k),', Length initial candidate set = ',num2str(length(iCAND))])
%     if k == 31
%         disp('Borrar esto')
%     end
    % 1) Deciding which points to eliminate
    RRR= sum(wADAPT,2) ;
    ilocPOS = find(RRR>0) ;    
    [~,indMIN] = min(RRR(ilocPOS)) ;
    indMINglo = ilocPOS(indMIN) ;
    jELIM = [jELIM; indMINglo] ;
    iCAND = setdiff(iCAND,indMINglo) ; % = [] ;
    % Now we loop over all clusters
    
    EXIT_WHILE = 0 ;
    
    for icluster = 1:length(U)
        %  disp( ['Cluster = ',num2str(icluster)])
        % Right-hand side   b-U(:,z)'*w_before_update
        Uloc = U{icluster}(iCAND,:) ;
        w_before= wADAPT(iCAND,icluster) ;
        w_eliminated = wADAPT(indMINglo,icluster)  ;
        RHS = b{icluster}-Uloc'*w_before ;
        if USE_LEAST_NORM == 0
            w_incre = Uloc'\RHS;
        else
            w_incre   = lsqminnorm(Uloc',RHS) ;
        end
        CHECK_w = sum(w_incre) -w_eliminated ;
        w_after = w_before + w_incre ;
        % disp(['Difference betwee sum(Delta w) and eliminated weight =',num2str(CHECK_w)])
        if any(w_after <0)
            warning('NEGATIVE WEIGHTS, exiting')
            EXIT_WHILE = 1 ;
            break
        end
        % Updating
        wADAPT(iCAND,icluster) = w_after ;
    end
    
    
    
    if EXIT_WHILE == 1
        break
    else
        wADAPT(jELIM,:) = 0  ;
    end
    
    
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
        maxval = max(wADAPT(iii,:))/VOL*100;
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
    set(htitle,'String',sprintf('Weights versus latent variable  — iteration %d',k));
    k = k+1;
    
end
RRR= sum(wADAPT,2) ;
ilocPOS = find(RRR>0) ;

disp(['Final number  of points = ',num2str(length(ilocPOS))])



