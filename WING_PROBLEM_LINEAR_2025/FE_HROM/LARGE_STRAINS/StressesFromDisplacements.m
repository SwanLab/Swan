function [GLSTRAINS,PK2STRESS,celastST,FgradST,PoneST,detFgrad] = StressesFromDisplacements(OPERFE,d,MATPRO,DATA)
% 2. Deformation gradient at all Gauss points
if nargin == 0
    load('tmp.mat')
end


ndof = size(OPERFE.Bst,2) ;
FgradST =  OPERFE.Bst*d(1:ndof,:)  ; %+ OPERFE.IDENTITY_F ;
if DATA.SMALL_STRAIN_KINEMATICS ==0
    for idim = 1:DATA.MESH.ndim
        LOCROWS = idim:DATA.MESH.ndim^2:length(FgradST) ;
        FgradST(LOCROWS,:) = 1+FgradST(LOCROWS,:) ;
    end
    % 3. Green-Lagrante strains at all Gauss points
    GLSTRAINS = StrainGreenLagrange(FgradST,DATA.MESH.ndim) ;
else
    if DATA.NO_USE_Deformation_gradient_in_Small_Strains ==0
        % In small strain kinematics, FgradST == gradU (because we haven't added the identity matrix )
        [GLSTRAINS,FgradST] = StrainGreenLagrange_small(FgradST,DATA.MESH.ndim) ;
    else
        % 8-Feb-2022. See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/PLASTICITY_HROM/README_plasticity.mlx
        GLSTRAINS = FgradST ; % Small strains are used, from the beginning
        FgradST = [] ;
    end
end
% 4. 2nd Piola-Kirchhoff stresses at all Gauss Points
if size(GLSTRAINS,2) == 1
    [PK2STRESS,celastST,detFgrad ]= PK2stress_Constitutive_Model(GLSTRAINS,MATPRO,DATA,FgradST) ;
else
    % This is used only for post-process purposes (not in the actual Newton-Raphson algorithm)
    PK2STRESS = zeros(size(GLSTRAINS)) ;   celastST = [] ; detFgrad = [] ;
    for  itime = 1:size(GLSTRAINS,2)
        if ~isempty(FgradST)
            Fgrad_inp = FgradST(:,itime)  ;
        else
            Fgrad_inp = [] ;
        end
        [PK2STRESS(:,itime)]= PK2stress_Constitutive_Model(GLSTRAINS(:,itime),MATPRO,DATA,Fgrad_inp) ;
    end
end

if isempty(PK2STRESS)
    PoneST = [] ;
else
    % 5. 1st Piola-Kirchhoff stresses at all Gauss Points
    if ~isempty(FgradST)
        PoneST = PK1stress(PK2STRESS,FgradST,DATA.MESH.ndim) ;
    else
        % This means that we are in the small strain regime.
        PoneST = [] ;% PK2STRESS ;
    end
    
end


