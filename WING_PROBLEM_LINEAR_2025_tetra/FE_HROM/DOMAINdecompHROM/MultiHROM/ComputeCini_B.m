function [ OPERFE]= ComputeCini_B(DATA,OPERFE,MATPRO,VAR)

if nargin == 1
    load('tmp3.mat')
end



if DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES == 1
    
    
    % Internal variables (value at previous time step)
    % ----------------------
    VARint_n = [];
    if ~isempty(DATA.ListFieldInternalVariables)
        for iINTVAR = 1:length(DATA.ListFieldInternalVariables)
            NameIntVar= DATA.ListFieldInternalVariables{iINTVAR} ;
            VARint_n.(NameIntVar) = VAR.(NameIntVar) ;
        end
    end
    SMALL_STRAIN_KINEMATICS = DATA.SMALL_STRAIN_KINEMATICS ;  % =
    
    if SMALL_STRAIN_KINEMATICS == 1
         DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES = 0 ;
        [~, OPERFE.celastST_ini,~,~] = StressesFromDisplacementsVAR(OPERFE,VAR,MATPRO,DATA,VARint_n) ;
        
    else
        % See  /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/108_EIFEM_metamat/04_EIFEM_NONL.mlx
        % We are interested in the linear part of the deformation
        % This is obtained by setting
        DATA.SMALL_STRAIN_KINEMATICS = 1;
        DATA.NO_USE_Deformation_gradient_in_Small_Strains = 0 ;
        DATA.MESH.ndim = DATA.MESH.ndimFINE ;
        DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES = 0 ;
        [~, celastST_ini,FgradST,~] = StressesFromDisplacementsVAR(OPERFE,VAR,MATPRO,DATA,VARint_n) ;
        
        % In contrast to the small-strains case, here we will turn
        % celastST_ini into a diagonal, block matrix, see more details in
        % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/108_EIFEM_metamat/04_EIFEM_NONL.mlx
        % See also /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/LARGE_STRAINS/KINEMATICS/KstiffLargeStrains.m
        
        OPERFE.celastST_ini = CelasLARGEmat_allgauss(celastST_ini,FgradST,DATA.MESH.ndim) ;
        
        % OPERFE.celastST_ini = ConvertBlockDiag(celasLARGE) ; % Diagonal block matrix
        
        
        
    end
    
    
end



%