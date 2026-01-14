function  Fint = InternalForces_COROT(OPERFE,PoneST,PK2STRESS,DATA,VAR,dCqLOC,KcLINq,BstQ,D_QrotALL) ;
% Adaptation of InternalForces.m to the corotational approach
% ------------------------------------------------------------
% % % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/02_PURE_ROTATION.mlx
% JAHO, 29-Oct-2024, UPC, CAMPUS NORD, Barcelona.

%
if nargin == 0
    load('tmp1.mat')
end


if ~isempty(PoneST)
    nF = DATA.MESH.ndim^2 ;
    if  DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES == 0
        error('Option not implemented')
        for icomp = 1:nF
            icol = icomp:nF:length(PoneST) ;
            PoneST(icol,:) = PoneST(icol,:).*OPERFE.wSTs;
        end
        
        Fint = OPERFE.Bst'*PoneST ;
    else        
        % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/108_EIFEM_metamat/04_EIFEM_NONL.mlx
        %          nF = DATA.MESH.ndim^2 ;
        %         for icomp = 1:nF
        %             icol = icomp:nF:length(PoneST) ;
        %             PoneST(icol,:) = VAR.PK1STRESS_incre(icol,:).*OPERFE.wSTs;
        %         end
        %         Fint = OPERFE.KstiffLINEAR*VAR.DISP + OPERFE.Bst'*PoneST ;
        % -------------------------------------------------------------------------------------------------------------------------          
        % \FintC  & =  \KcLINq \dC   -  {\LboolCall^T  \DiagC{\QrotALL}   \DiagC{\KcLIN}      \dCqLOC  }
        %   \\&  + \BstQ^T \DiagC{\Wecm} (\PoneFlocST - \PoneFlocLINst)
        %
        Fint = KcLINq*VAR.DISP ;
        Fint = Fint - OPERFE.LboolCall'*(D_QrotALL*(OPERFE.D_KcLINloc*dCqLOC));
        Fint = Fint +    BstQ'*OPERFE.D_Wecm*VAR.PK1STRESS_incre ;
    end
    
    
    
    
else
    
    error('Options not implemented')
    
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

