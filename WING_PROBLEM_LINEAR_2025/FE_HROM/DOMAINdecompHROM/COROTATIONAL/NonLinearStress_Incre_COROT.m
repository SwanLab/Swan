function VAR = NonLinearStress_Incre_COROT(DATA,OPERFE,VAR,GgradFlocST)
% Adaptation of NonLinearStress.m to the corotational approach
% % % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/02_PURE_ROTATION.mlx
% JAHO, 29-Oct-2024, Green's Aribau, Barcelona.  
% -----------------------------------------------------------------------------------------------------------------------

if DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES == 1
    
    
    if DATA.SMALL_STRAIN_KINEMATICS==1 
        error('Option not implemented yet, 29-Oct-2024')
        % THIS IS FOR SMALL STRAINS
        % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx
        %   celastST_iniD = ConvertBlockDiag(OPERFE.celastST_ini) ; % Diagonal block matrix
        % OPERFE.celastINI_Bst= celastST_iniD*OPERFE.Bst;
        
        %         OPTION = 0;
        %
        %         if OPTION == 1
        STRESS_elastic = zeros(size(VAR.PK2STRESS)) ;
        nstrain = size(OPERFE.celastST_ini,2) ;
        for istrain = 1:nstrain
            istrainCOL = istrain:nstrain:size(OPERFE.celastST_ini,1) ;
            for jstrain = 1:nstrain
                jstrainCOL = jstrain:nstrain:size(OPERFE.celastST_ini,1) ;
                STRESS_elastic(istrainCOL) =   STRESS_elastic(istrainCOL) + OPERFE.celastST_ini(istrainCOL,jstrain).*VAR.GLSTRAINS(jstrainCOL) ;
            end
        end
        VAR.PK2STRESS_incre = VAR.PK2STRESS - STRESS_elastic ;
        
        %         else
        %
        %             VAR.PK2STRESS_incre = VAR.PK2STRESS  -OPERFE.celastINI_Bst*VAR.DISP  ;
        %         end
        
        
    else
        
          % PoneFlocLINst =   \DiagC{\CtangFlinST} \GgradFlocST
          PoneFlocLINst = OPERFE.D_CtangFlinST*GgradFlocST; 
        
        VAR.PK1STRESS_incre = VAR.PK1STRESS - PoneFlocLINst ;
        
        
        
        
        
        %   large strains
        %
        %         METHOD_loop = 1;
        %
        %         if METHOD_loop == 1
        
%         STRESS_elastic = zeros(size(VAR.PK1STRESS)) ;
%         nstrain = size(OPERFE.celastST_ini,2) ;
%         for istrain = 1:nstrain
%             istrainCOL = istrain:nstrain:size(OPERFE.celastST_ini,1) ;
%             for jstrain = 1:nstrain
%                 jstrainCOL = jstrain:nstrain:size(OPERFE.celastST_ini,1) ;
%                 
%                 
%                 STRESS_elastic(istrainCOL) =   STRESS_elastic(istrainCOL) ...
%                     + OPERFE.celastST_ini(istrainCOL,jstrain).*(FgradST(jstrainCOL)-OPERFE.IDENTITY_F(jstrainCOL)) ;
%                 % 9-oct-2024...SHOULDN T WE SYMMETRIZED THIS? ECM WORKS
%                 % WITH SYMMETRIZED GRADIENTS, SEE
%                 % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/DOMAINdecompHROM/MultiHROM/MultiSnapStressFromDispNECMlarg.m
%                 
%             end
%         end
       
       
        
        %         else
        %             % SEcond method
        %             celastST_ini = ConvertBlockDiag(OPERFE.celastST_ini) ;
        %
        %             VAR.PK1STRESS_incre = VAR.PK1STRESS - celastST_ini*(OPERFE.Bst*VAR.DISP) ;
        %
        %         end
        
        
        
        
        
        
        
        
    end
    
    
end
