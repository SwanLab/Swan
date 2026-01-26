function [PhiDEF,PsiDEFf,Vall,Mintf,Mintf_inv,Mdom,PhiRB,PsiRBf,MESH,DATA] ...
    = MODESandROMoper_EIFE(OPERFE,MESH,DATA,SNAPreact,BasisStwo,...
    SNAPdisp,DATAcommon,DATAoffline)
% DETERMINATION OF MODES FOR A GIVEN SUBDOMAIN of the corresponding deformational, reactive and interface modes
% Modification of MODESandROMoper_1dom 
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/06_3D_27_HEXAHEDRA/04_FinalExam.mlx
     % WE FOUND 
     % that the version used in giuliudory2023multiscale.pdf (MODESandROMoper_1dom) is prone to
     % instabilities. The version implemented here attempts to eliminate
     % such instabilities
% JAHO, 1-APRIL-2023
% ---------------------------------------------------------------------------------
if nargin == 0
    load('tmp.mat')
    
end

% STEP 1
% -------
% Nodes and DOFs of the interface boundaries. 
[faceDOFS,PhiRB,PsiRBf,Mdom,Mintf,MESH] =   GeometricVarDOMAINS(OPERFE,MESH,DATA,DATAcommon);
% PhiRB = Rigid body modes
% PsiRBf  =Resultant modes 
% Mdom --> Geometric mass matris for the domain 
% MdomFF --> Geometric mass matrix interface 
% MESH > COOR,CN,... of the PARENT DOMAIN 

% STEP 2
% ------------
% Interface modes 
PhiDEF = [] ; 
PsiDEFf = [] ; 
[Vall,MESH,DATA,INFOPLOTMESHBND] = FictInterfaceDISP(PhiDEF,PsiDEFf,MESH,DATAcommon,DATA,Mintf,PhiRB,PsiRBf,DATAoffline);


% STEP 3 AND 4) SELF-EQUILIBRATED MODES 

 [PsiDEFf,Mintf_inv]= SelfEquilibratedModes_EIFE(MESH,Vall,PhiRB,PsiRBf,Mintf,SNAPreact,DATA,INFOPLOTMESHBND) ; 

 
% STEP 4
% ------------
% Deformational modes
% ---------------------------------------------------
disp('Deformational modes...')

[PhiDEF,PsiDEFf] =   DeformationalModes_EIFE(MESH.faceDOFSall,SNAPdisp,PsiDEFf,PhiRB,Mdom)  ; ;

disp(['Number of deform/self-equilibrated modes (final) : ',num2str(size(PsiDEFf,2))])
disp('------------------------')

 


% PLOTTING MODES
 % --------------------------------------------------------------------- 
    PlotModesDEF_SubdomainLevel(DATA,[PhiRB,PhiDEF],MESH)
  
% -------------------------------------------------
% POST-PROCESSING SELF-EQUILIBRATED MODES
% -------------------------------------------------
 PlotModesSE_SubdomainLevel(DATA,PsiRBf,PsiDEFf,MESH)


% STEP 5
% ------------
% Displacement interface modes: rigid-body/Shape functions and candidates for fluctuation
% modes
% disp('Displacement interface modes...')
% 


%  DATAinputLOC = DefaultField(DATAinputLOC,'USE_ONLY_INFO_SLICE',[] ) ; % = [50] ;
% 
% save(DATAinputLOC.NAME_STORE_BASIS_INFORMATION_FOR_OTHER_PROJECTS,'PsiDEFf','MdomFFinv','PhiDEF','Vall','Vrb','Mintf') ;
% 
% 
% 
% % Taking
% PhiDEF_all = cell(size(DATAROM)) ;
% PsiDEFf_all = cell(size(DATAROM)) ;
% 
% 
% PhiDEF_all(:) = {PhiDEF}  ;
% PsiDEFf_all(:) = {PsiDEFf}  ;
% nINTF = length(unique(cnINTF(:))) ;
% Vall_all = cell(1,nINTF) ;
% Mintf_all = cell(1,nINTF) ;
% Vrb_all = cell(1,nINTF) ;
% 
% Vall_all(:) = {Vall}  ;
% Mintf_all(:) = {Mintf{1}}  ;
% Vrb_all(:) = {Vrb{1}}  ;
% 
% 
