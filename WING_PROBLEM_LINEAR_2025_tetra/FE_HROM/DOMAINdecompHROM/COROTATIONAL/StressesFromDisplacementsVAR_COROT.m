function [VAR,celastST,FgradST,detFgrad] = StressesFromDisplacementsVAR_COROT(OPERFE,VAR,MATPRO,DATA,VARint_n,...
    dCqLOC,BstQ,D_QrotALL)
% Adaptation of StressesFromDisplacementsVAR.m  to the corotational
% approach
% % % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/02_PURE_ROTATION.mlx
% JAHO, 29-Oct-2024, Green's Aribau, Barcelona.


% See
if nargin == 0
    load('tmp1.mat')
end

%
ndof = DATA.MESH.ndof ;
% FgradST =  OPERFE.Bst*VAR.DISP(1:ndof,:)  ; % Why this 1:ndof ? ????
% 29-Oct-2024

% Gradient of incremental displacements
%  \GgradFlocST  =  \BstQ \dC - \DiagC{\BmatIst}   \dCqLOC
GgradFlocST = BstQ*VAR.DISP(1:ndof,:) - OPERFE.D_BmatIst*dCqLOC ;
% Deformation gradient
FgradST = GgradFlocST ;

if DATA.SMALL_STRAIN_KINEMATICS ==0
    for idim = 1:DATA.MESH.ndim
        LOCROWS = idim:DATA.MESH.ndim^2:length(GgradFlocST) ;
        % Deformation gradient
        FgradST(LOCROWS,:) = 1+GgradFlocST(LOCROWS,:) ;
    end
    % 3. Green-Lagrante strains at all Gauss points
    VAR.GLSTRAINS = StrainGreenLagrange(FgradST,DATA.MESH.ndim) ;
else
    error('Option not implemented (29-Oct-2024)')
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

if DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES == 1
    PoneFlocLINst = OPERFE.D_CtangFlinST*GgradFlocST;
    VAR.PK1STRESS_incre = VAR.PK1STRESS - PoneFlocLINst ;
else
    error('Option not implemented')
end

%
% VAR = NonLinearStress_Incre_COROT(DATA,OPERFE,VAR,GgradFlocST) ;
%

