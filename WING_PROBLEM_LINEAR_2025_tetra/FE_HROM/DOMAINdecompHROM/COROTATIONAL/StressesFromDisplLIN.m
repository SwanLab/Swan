function [VAR,celastST,FgradST,detFgrad] = StressesFromDisplLIN(OPERFE,VAR,MATPRO,DATA,VARint_n)
% Copy of StressesFromDisplacementsVAR.m
% Adapted for co-rotational approach. Evaluation of linear constitutive tensors,  
 % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/02_PURE_ROTATION.mlx
 % 28-Oct-2024, Balmes 185, Barcelona 
if nargin == 0
    load('tmp1.mat')
end

% Gradient of displacements
%  \GgradFlocST  =  \BstQ \dC
%  \BstQ \defeq \DiagC{\BmatIst}  \DiagC{\QrotALL^T}  \LboolCall; 
% Displacement at element level
% Rotation matrix 
%dCeST =   QrotINIall_GLOBAL(OPERFE.QrotINI,OPERFE.INDEXsparseROTmat); 

%ndof = size(OPERFE.Bst,2) ;
FgradST =  zeros(DATA.MESH.ndofSTRESS,1)  ; %+ OPERFE.IDENTITY_F ;
if DATA.SMALL_STRAIN_KINEMATICS ==0
    for idim = 1:DATA.MESH.ndim
        LOCROWS = idim:DATA.MESH.ndim^2:length(FgradST) ;
        FgradST(LOCROWS,:) = 1+FgradST(LOCROWS,:) ;
    end
    % 3. Green-Lagrante strains at all Gauss points
    VAR.GLSTRAINS = StrainGreenLagrange(FgradST,DATA.MESH.ndim) ;
else
    if DATA.NO_USE_Deformation_gradient_in_Small_Strains ==0
        % In small strain kinematics, FgradST == gradU (because we haven't added the identity matrix )
        [VAR.GLSTRAINS,FgradST] = StrainGreenLagrange_small(FgradST,DATA.MESH.ndim) ;
    else
        % 8-Feb-2022. See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/PLASTICITY_HROM/README_plasticity.mlx
        VAR.GLSTRAINS = FgradST ; % Small strains are used, from the beginning
        FgradST = [] ;
    end
end
% 4. 2nd Piola-Kirchhoff stresses at all Gauss Points
DATA.CALC_CTANG = 1 ;
if size(VAR.GLSTRAINS,2) == 1
    %  [PK2STRESS,celastST,detFgrad ]= PK2stress_Constitutive_Model(GLSTRAINS,MATPRO,DATA,FgradST) ;
    [VAR,celastST,detFgrad ]= PK2stress_Constitutive_ModelVAR(VAR,MATPRO,DATA,FgradST,VARint_n,OPERFE) ;
else
    % This is used only for post-process purposes (not in the actual Newton-Raphson algorithm)
    % error('ADapt this option to the new framework (9-Feb-2022)')
    GLSTRAINS  = VAR.GLSTRAINS ;
    PK2STRESS = zeros(size(GLSTRAINS)) ;   celastST = [] ; detFgrad = [] ;
    if ~isempty(VARint_n)
        FFF = fieldnames(VARint_n)   ;
        %   VARINTtime = [] ;
        for ivar = 1:length(FFF)
            FLOC = FFF{ivar} ;
            VARint_n.(FLOC) = VAR.(FLOC) ;
            % VARINTtime.(FLOC) = zeros(size(VAR.(FLOC),1),size(GLSTRAINS,2)) ;
        end
    end
    % VAR.GLSTRAINS
    DATA.CALC_CTANG = 0 ;
    DATA.kiter = 2;
    for  itime = 1:size(GLSTRAINS,2)
        VAR.GLSTRAINS = GLSTRAINS(:,itime) ;
        if isempty(FgradST)
            FgradSTloc = [] ;
        else
            FgradSTloc = FgradST(:,itime) ;
        end
        [VAR,celastST,detFgrad ]= PK2stress_Constitutive_ModelVAR(VAR,MATPRO,DATA,FgradSTloc,VARint_n) ;
        PK2STRESS(:,itime) = VAR.PK2STRESS ;
        
        % Internal variables
        if ~isempty(VARint_n)
            FFF = fieldnames(VARint_n)   ;
            for ivar = 1:length(FFF)
                FLOC = FFF{ivar} ;
                VARint_n.(FLOC) = VAR.(FLOC) ;
                
            end
        end
        
        
    end
    VAR.PK2STRESS = PK2STRESS ;
end

if isempty(VAR.PK2STRESS)
    VAR.PK1STRESS = [] ;
else
    % 5. 1st Piola-Kirchhoff stresses at all Gauss Points
    if ~isempty(FgradST)
        VAR.PK1STRESS = PK1stress(VAR.PK2STRESS,FgradST,DATA.MESH.ndim) ;
    else
        % This means that we are in the small strain regime.
        VAR.PK1STRESS = [] ;% PK2STRESS ;
    end
    
end

%  For the case in which DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES == 1
VAR = NonLinearStress_Incre(DATA,OPERFE,VAR,FgradST) ;


