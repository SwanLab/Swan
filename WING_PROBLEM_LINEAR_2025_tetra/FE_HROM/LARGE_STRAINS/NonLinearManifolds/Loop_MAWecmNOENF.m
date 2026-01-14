function [wADAPT,iCAND,jELIM,kLOC] = Loop_MAWecmNOENF(U,ilocPOS,iCAND,wADAPT,b,icluster,indSORT,jELIM)

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
        w_before= wADAPT(iCAND_new,icluster) ;
        %w_eliminated = wADAPT(indMINglo,icluster)  ;
       % RHS = b{icluster}-Uloc'*w_before ;
        %             USE_CHT = 1;
        %             if USE_CHT == 0
        %                 if USE_LEAST_NORM == 0
        %                     w_incre = Uloc'\RHS;
        %                 else
        %                     w_incre   = lsqminnorm(Uloc',RHS) ;
        %                 end
        %                 CHECK_w = sum(w_incre) -w_eliminated ;
        %                 w_after = w_before + w_incre ;
        %             else
        
        [ w_after,SOLUTION_fOUND] = nn_update_active_set(Uloc, b{icluster}, w_before) ;
        %    errorLOC = Uloc'*w_after-b{icluster}   ;
        %   errorLOC = norm(errorLOC)/norm(b{icluster}) ;
        %             if errorLOC >1e-10
        %                  disp(['error_appox, rel, over 1 =',num2str(errorLOC)])
        %             end
        %   end
        
        
        
        
        % disp(['Difference betwee sum(Delta w) and eliminated weight =',num2str(CHECK_w)])
        if SOLUTION_fOUND == 0
            disp(['No feasible solution found, exiting, icluster = ',num2str(icluster)])
            EXIT_WHILE_cluster = 1 ;
            break
        end
        % Updating
        wADAPT_new(iCAND_new,icluster) = w_after ;
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