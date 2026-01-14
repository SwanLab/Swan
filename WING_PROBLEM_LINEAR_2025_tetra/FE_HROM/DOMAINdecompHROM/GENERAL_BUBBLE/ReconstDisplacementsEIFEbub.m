function [DISP,strainC] = ReconstDisplacementsEIFEbub(dCOARSEloc,EIFE_prop,Vrot,DATA,SCALE_FACTOR_ALLS,ROTATIONmat,ndim)
% JAHO, 16-oct-2023, Barcelona 
% Original file explained in: /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/07_PostProcess_2D.mlx
% Adaptation to bubble modes: /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/04_GeneralTheory.mlx
if nargin == 0
    load('tmp.mat')
end
% deformational part (amplitude modes)

% BOUNDARY DOFS AND COARSE-SCALE DOFS 
% Defined in /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/DOMAINdecompHROM/GENERAL_BUBBLE/EIFE_operatorsBUB.m
% EIFEoper.INFO.DOFsBUB = DOFsBUB; % BUBBLE DOFS
% EIFEoper.INFO.DOFsB= DOFsB; % boundary DOFS 
DOFsBUB= EIFE_prop.INFO.DOFsBUB; 
DOFsB= EIFE_prop.INFO.DOFsB; 
IndexesStrainModes= EIFE_prop.INFO.IndPhiDEF; 
IndexesBubbleModes= EIFE_prop.INFO.IndGammaBUB; 


DISP = zeros(size(EIFE_prop.RECONSTRUCTION.DEF_DISP.BASIS,1),size(dCOARSEloc,2)) ; 
% 1) DEFORMATIONAL PART, BOUNDARY NODES % ---------------------------------------------------------------------
strainC = EIFE_prop.OPER.HdefINV_PsiDEFfT*Vrot*dCOARSEloc(DOFsB,:) ;   %%  OPER.HdefINV_PsiDEFfT = Ldef
DISP = DISP + EIFE_prop.RECONSTRUCTION.DEF_DISP.BASIS(:,IndexesStrainModes)*strainC ; 
% 2) DEFORMATIONAL PART, bubble DOFS % ---------------------------------------------------------------------
 
DISP = DISP + EIFE_prop.RECONSTRUCTION.DEF_DISP.BASIS(:,IndexesBubbleModes)*dCOARSEloc(DOFsBUB,:) ; 
% 3) rigid body  PART, BOUNDARY NODES % ---------------------------------------------------------------------
aRBb = EIFE_prop.OPER.Ginv_PhiRBfT_m_HrbHdefINV_PsiDEFfT*Vrot*dCOARSEloc(DOFsB,:) ;   %%  OPER.Ginv_PhiRBfT_m_HrbHdefINV_PsiDEFfT = Lrb
DISP = DISP + EIFE_prop.RECONSTRUCTION.RB_DISP.BASIS*aRBb ; 
% 4) rigid body  PART, BUBBLE DOFS % ---------------------------------------------------------------------
aRBbub = EIFE_prop.OPER.LrbBUB*dCOARSEloc(DOFsBUB,:) ;   %%  OPER.Ginv_PhiRBfT_m_HrbHdefINV_PsiDEFfT = Lrb
if ~isempty(aRBbub)
DISP = DISP + EIFE_prop.RECONSTRUCTION.RB_DISP.BASIS*aRBbub ; 
end

%  

if DATA.UNIFORM_SCALING_REFERENCE_ELEMENT == 1
    
    %    dREFrot = SCALE_FACTOR_ALLS*ROTATIONmat*reshape(dREF,ndim,[]) ;
    
    
     
    for itime = 1:size(DISP,2)
        
        dREFrot = ROTATIONmat*reshape(DISP(:,itime),ndim,[]) ;
        DISP(:,itime) = dREFrot(:) ;
        
    end
    
    
    
else
    error('Option not implemented')
end