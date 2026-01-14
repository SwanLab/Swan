function  Fint = InternalForces(OPERFE,PoneST,PK2STRESS,DATA,VAR) ;  
% Internal forces, large strains, in terms of the Bst, W operators, as well
% as the 1st PK stress (stacked) 
% 
if nargin == 0
    load('tmp.mat')
end


if ~isempty(PoneST)
    
    
    
    nF = DATA.MESH.ndim^2 ;
    if  DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES == 0
        for icomp = 1:nF
            icol = icomp:nF:size(PoneST,1) ;
            PoneST(icol,:) = PoneST(icol,:).*OPERFE.wSTs;
        end
        
        Fint = OPERFE.Bst'*PoneST ;
    else
       
        % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/108_EIFEM_metamat/04_EIFEM_NONL.mlx
         nF = DATA.MESH.ndim^2 ;
        for icomp = 1:nF
            icol = icomp:nF:length(PoneST) ;
            PoneST(icol,:) = VAR.PK1STRESS_incre(icol,:).*OPERFE.wSTs;
        end
        Fint = OPERFE.KstiffLINEAR*VAR.DISP + OPERFE.Bst'*PoneST ;
         
    end
    
    
    
    
else
   
    
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

