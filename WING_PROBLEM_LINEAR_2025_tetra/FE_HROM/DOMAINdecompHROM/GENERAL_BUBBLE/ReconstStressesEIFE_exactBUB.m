function   [STRESSES_FINE,INTERNAL_VAR,VMSTRESS] = ReconstStressesEIFE_exactBUB(EIFE_prop,dCOARSE,DATA,...
    SCALE_FACTOR_ALLS,PROPMAT,Vrot,ROTATIONmat)
% Reconstruction of stresses and internal variables by solving integrating
% the constitutive equations along time
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/05_VECTORIZ/02_PLASTICITY.mlx
% See  /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/04_GeneralTheory.mlx
% JAHO, 23-March-2023, 16-Oct-2023 (adaptation to Bubble modes)
if nargin == 0
    load('tmp1.mat')
end

% DETERMINATION OF VARIABLE "MATPRO" FOR ALL THE GAUSS POINTS OF THE EIF
% ELEMENT (in the parent configuration)
DATA.MESH.ngaus = size(EIFE_prop.MESH.posgp,2) ;
[MATPRO,MESH,DATA,INICOND] = MaterialPropertiesIntegPoints(EIFE_prop.MESH,DATA,PROPMAT) ;


if isempty(DATA.ListFieldInternalVariables)
    VARint_n = [] ;
else
    DATA = DefaultField(DATA,'NumberInternalVariableToPrint',[length(DATA.ListFieldInternalVariables)]);
    for ivar = 1:length(DATA.ListFieldInternalVariables)
        FLOC = DATA.ListFieldInternalVariables{ivar};
        VARint_n.(FLOC) = INICOND.(FLOC)  ;
        DATA.STORE.VAR.(FLOC) = 1;
    end
end
% STRESSES FROM HISTORY OF DISPLACEMENTS
% See
% Bst = Bmat_PhiDEF*EIFE_prop.OPER.HdefINV_PsiDEFfT*Vrot/SCALE_FACTOR_ALLS;
% EIFEoper.RECONSTRUCTION.STRAINS.BASIS = OPERFE.Bst*MODES.PhiDEF ;

% --------------------
% OLD IMPLEMENTATION -
% --------------------
%OPERFE.Bst = (EIFE_prop.RECONSTRUCTION.STRAINS.BASIS)*(EIFE_prop.OPER.HdefINV_PsiDEFfT*Vrot)/SCALE_FACTOR_ALLS ;

% --------------------
% NEW IMPLEMENTATION -
% --------------------
IndexesStrainModes= EIFE_prop.INFO.IndPhiDEF;
IndexesBubbleModes= EIFE_prop.INFO.IndGammaBUB;

BstB=   (EIFE_prop.RECONSTRUCTION.STRAINS.BASIS(:,IndexesStrainModes))*(EIFE_prop.OPER.HdefINV_PsiDEFfT*Vrot)/SCALE_FACTOR_ALLS ;
BstBUB=   (EIFE_prop.RECONSTRUCTION.STRAINS.BASIS(:,IndexesBubbleModes))/SCALE_FACTOR_ALLS ; % INCLUDE SCALE_FACTOR_ALLS IN 6-nOV-2024
OPERFE.Bst = [BstB,BstBUB] ;



VAR.DISP = dCOARSE ;

[VAR,VARINTtime,~,~,~] = StressesFromDisplacementsVAR_EIFE(OPERFE,VAR,MATPRO,DATA,VARint_n) ;
STRESSES_FINE_REF = VAR.PK2STRESS;

ndim = size(MESH.COOR,2) ;
VMSTRESS = VonMisesStressCOMP(STRESSES_FINE_REF,ndim,DATA) ;

if ~isempty(DATA.ListFieldInternalVariables)
    NameInternalVariable = DATA.ListFieldInternalVariables{DATA.NumberInternalVariableToPrint} ;
    INTERNAL_VAR = VARINTtime.(NameInternalVariable)  ;
else
    INTERNAL_VAR = [] ;
end


if EIFE_prop.MESH.nstrain == 4
    
    STRESSES_FINE =zeros(size(STRESSES_FINE_REF)) ;
    if length(ROTATIONmat) == 1
        ROTATION_STRESSES = eye(3) ;
    else
        ROTATION_STRESSES =    RotateStress2Dplanestrain(ROTATIONmat(1,1),ROTATIONmat(2,1)) ;
    end
    
    for itime = 1:size(STRESSES_FINE_REF,2)
        STRESSES_FINE_loc =  reshape(STRESSES_FINE_REF(:,itime),EIFE_prop.MESH.nstrain,[])  ;
        
        STRESSES_FINE_loc(1:3,:) = ROTATION_STRESSES*STRESSES_FINE_loc(1:3,:) ;
        
        % GID CONVENTION FOR PRINTING
        %   STRESSES_FINE = STRESSES_FINE([1 2 4 3],:) ;
        STRESSES_FINE(:,itime) = STRESSES_FINE_loc(:);
        %  ngaus =size(EIFE_prop.MESH.posgp,2) ;
        %  nstrain = (EIFE_prop.MESH.nstrain) ;
        
        %  STRESSES_FINE_loc =  reshape(STRESSES_FINE_loc,nstrain*ngaus,[])  ;
        
    end
    
    
else
    % ROTATION OF FINE-SCALE STRESSES NOT IMPLEMENTED YET
    STRESSES_FINE = STRESSES_FINE_REF ;
    disp('Rotation of stresses not implemented yet in 3D')
end



