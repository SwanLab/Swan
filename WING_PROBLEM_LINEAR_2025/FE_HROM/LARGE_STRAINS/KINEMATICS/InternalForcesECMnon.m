function  [Fint,fINTredNON_mst_DOFl,Bst_non] = InternalForcesECMnon(OPERFE,PoneST,PK2STRESS,DATA,VAR,tauNONallDER_q) ;  
% Internal forces, large strains, in terms of the Bst, W operators, as well
% as the 1st PK stress (stacked) 
% 
if nargin == 0
    load('tmp3.mat')
end


if ~isempty(PoneST)
    
    
    nF = DATA.MESH.ndim^2 ;
    % PoneST is the 1st Piola-Kirkhoff stresses for all the ECM points (master ECM points)
    % PoneST = [PoneST_1; PoneST_2 ; \vdots PoneST_mECM_master]
    % What we need is the internal work density matrix   \fINTredNON_{mst}
    % such that 
    % (\FintREDnon{}{})_{mst} =   \fINTredNON_{mst} w_{mst}
    % Here w_{mst} is OPERFE.wSTs
    % The size of \fINTredNON_{mst} is:  as many columns as master ECM
    % points
    % What about the rows ?  Number of  generalized coordinates, including the constrained ones. 
         
    Bst_non = OPERFE.Bst*tauNONallDER_q ; 
    
    nREDall = size(Bst_non,2) ; 
    mMST =length(OPERFE.wSTs); 
    fINTredNON_mst = zeros(nREDall,mMST);
    for icomp = 1:nF
        icol = icomp:nF:size(PoneST,1) ;
        
        BstLOC = Bst_non(icol,:) ; 
        % The size of the above matrix is
        %  mMST x nmodesUall 
        PoneSTloc = PoneST(icol) ; 
        % The size of the above matrix mMSTx1
        % Thus, we have to compute 
        fINTloc = bsxfun(@times,BstLOC,PoneSTloc); 
        % fINTloc is now  a mMST x nmodesUall matrix 
        
        % Now we just add the contribution
        fINTredNON_mst = fINTredNON_mst + fINTloc'; 
    end    
        % The end result 
        %  (\FintREDnon{}{})_{mst} =   \fINTredNON_{mst} w_{mst}
       Fint  = fINTredNON_mst*OPERFE.wSTs ; 
       
       % NExt step: we use just the unconstrained component of FintREDnon_mst
       fINTredNON_mst_DOFl = fINTredNON_mst(OPERFE.DOFl,:) ; 
       % as an input for the nonlinear mapping called "etanon"         
       % \fINTredNON_{slv,\DOFl} =\etaNON(\fINTredNON_{mst,\DOFl})
       fINTredNON_slv_DOFl = OPERFE.etaNON(fINTredNON_mst_DOFl)' ; 
       
       
       
       % And then multiplied by their weights to get the contribution
       % to the reduced internal forces 
       FintREDnon_slv_DOFl = fINTredNON_slv_DOFl*OPERFE.wRED_slv ; 
       % Now we add the slave contribution to the one from the master
       % integration points
       Fint(OPERFE.DOFl) = Fint(OPERFE.DOFl) + FintREDnon_slv_DOFl ; 
       
    
    
    
else
   error('Option not implemented yet, 21-July-2025')
    
    if  DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES == 0
        nF = DATA.MESH.nstrain ;
        for icomp = 1:nF
            icol = icomp:nF:length(PK2STRESS) ;
            PK2STRESS(icol,:) = PK2STRESS(icol,:).*OPERFE.wSTs;
        end
        Fint = OPERFE.Bst'*PK2STRESS ;
    else
        % Special implementation for EIFEM
        % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx
        nF = DATA.MESH.nstrain ;
        for icomp = 1:nF
            icol = icomp:nF:length(PK2STRESS) ;
            PK2STRESS(icol,:) = VAR.PK2STRESS_incre(icol,:).*OPERFE.wSTs;
        end
        Fint = OPERFE.KstiffLINEAR*VAR.DISP + OPERFE.Bst'*PK2STRESS ;
    end
end

