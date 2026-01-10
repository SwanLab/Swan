function [ PsiDEFf_aligned,Mintf_inv ]= SelfEquilibratedModes_EIFE(MESH,Vall,PhiRB,PsiRBf,Mintf,SNAPreact,DATA,INFOPLOTMESHBND)
% METHOD FOR CALCULATING THE SELF-EQUILIBRATED MODES FOR THE EIFE METHOD 
%JAHO, 2-APRIL-2023
% sEE /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/06_3D_27_HEXAHEDRA/04_FinalExam.mlx
% ------------------------------------------------------------------------
if nargin == 0
    load('tmp.mat')
end

% STEP 3: DEFORMATIONAL PART OF THE DISPLACEMENT INTERFACE MODES Vall 
% -------------------------------------------------------------------

Vall_deform = DefComponentInterfaceModes(MESH,PhiRB,Vall,Mintf,INFOPLOTMESHBND,DATA) ; 
disp('------------------------')
disp(['Number of deformational interface modes = ',num2str(size(Vall_deform,2))])
disp('------------------------')


% STEP 4:  Self-equilibrated modes  (standard calculation)
% -----------------------------------------
disp('------------------------')
disp('Self-equilibrated modes...')
[PsiDEFf,Mintf_inv,~,Mintf_inv_chol] =   SelfEquilibratedModes_1dom(MESH.faceDOFSall,PsiRBf,Mintf,SNAPreact,DATA)  ;
disp(['Number of self-equilibrated modes (before computing   conjugate to deformational interface modes) : ',num2str(size(PsiDEFf,2))])


% STEP 5: SELF-EQUILIBRATED MODES MORE ALIGNED WITH Vall_deform% 
% ---------------------------------------------------------------
% PROJECTION Vall_deform ONTO THE SPAN OF PsiDEFf
COEFF = PsiDEFf'*Vall_deform ; 
[uu,ss,vv] = SVDT(COEFF) ; 
% SINCE PsiDEFf is orthogonal in Mff^-1, and Vall_deform to Mff, then ss
% are all ranged between 0 and 1, and can be interpeted as 
%  the cosines of the principle angles formed by the subspaces
% Vall_deform and PsiDEFf
% SUSPICIOUSLY LOW VALUES OF ss may imply ill-conditioning 
TOL_LIMIT_WELL_POSEDNESS  = 1e-2  ;
if any(ss < TOL_LIMIT_WELL_POSEDNESS)
    
     NameLoc =     'NON_EXCITED_MODES' ;
    NAME_MODES_FOLDER = [cd,filesep,'MODES',filesep];
    NAME_MODES_DISP = [NAME_MODES_FOLDER,DATA.NAME_BASE,NameLoc] ;
    NameFileMesh = [NAME_MODES_DISP,'.msh'] ;
    NameFile_res = [NAME_MODES_DISP,'.res'] ;
    DATALOC = [];
    
    indices=  find(ss < TOL_LIMIT_WELL_POSEDNESS) ;
    
    Vall_deformPLOT = Vall_deform*vv(:,indices) ; 
    
    
%     nv = sqrt(sum(Vall_deform.^2,1)) ;
%     Vall_deformPLOT = bsxfun(@times,Vall_deformPLOT',1./nv')' ; 
    
    GidPostProcessModesDOML(INFOPLOTMESHBND.COOR, INFOPLOTMESHBND.CN ,INFOPLOTMESHBND.TypeElement,Vall_deformPLOT,[],NameFileMesh,NameFile_res,[],DATALOC) ;
    
    
    error('TRAINING SET NOT APPROPRIATE. THERE ARE ONE OR MORE DEFORMATIONAL INTERFACE MODES THAT DOES NOT CONTRIBUTE TO THE WORK DONE BY THE REACTIVE MODES')
end
% Now we have to determine A BASIS MATRIX FOR  the subspace of SPAN(PsiDEFf) which is more
% aligned to Vall_deform. This is given by:
DATALOCw.Mchol = Mintf_inv_chol ; 
PsiDEFf_aligned = WSVDT(PsiDEFf*COEFF,[],DATALOCw) ;   %  