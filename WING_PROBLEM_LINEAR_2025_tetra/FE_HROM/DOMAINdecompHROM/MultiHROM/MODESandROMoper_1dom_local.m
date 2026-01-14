function [PhiDEF,PsiDEFf,Vall,Mintf,Mdom,PhiRB,PsiRBf,MESH,DATA] ...
    = MODESandROMoper_1dom_local(OPERFE,MESH,DATA,SNAPreact,BasisStwo,...
    SNAPdisp,DATAcommon,DATAoffline)
% DETERMINATION OF MODES FOR A GIVEN SUBDOMAIN
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/README_HOMOG.mlx
% JAHO, 10-FEB-2023

if nargin == 0
    load('tmp.mat')
    
end

% STEP 1
% -------
% Nodes and DOFs of the interface boundaries
[faceDOFS,PhiRB,PsiRBf,Mdom,Mintf,MESH] =   GeometricVarDOMAINS(OPERFE,MESH,DATA,DATAcommon);

% STEP 2
% ------------
% Self-equilibrated modes


disp('------------------------')
disp('Self-equilibrated modes...')
DATA.TOL = DATAoffline.errorSVD_SELF_EQUILIBRATED ; 

DATAoffline = DefaultField(DATAoffline,'NumberOfDeformationalModes',[]) ; 

if ~isempty(DATAoffline.NumberOfDeformationalModes)
    DATA.TOL = 0 ; 
end

 

[PsiDEFf,MintfINV] =   SelfEquilibratedModes_1dom(MESH.faceDOFSall,PsiRBf,Mintf,SNAPreact,DATA)  ;

if ~isempty(DATAoffline.NumberOfDeformationalModes)
     nmodes= min(DATAoffline.NumberOfDeformationalModes,size(PsiDEFf,2)) ; 
     PsiDEFf = PsiDEFf(:,1:nmodes); 
end



disp(['Number of self-equilibrated modes (before computing deformational conjugate) : ',num2str(size(PsiDEFf,2))])


 
% STEP 4
% ------------
% Deformational modes
disp('Deformational modes...')

[PhiDEF,PsiDEFf] =   DeformationalModes_1dom(MESH.faceDOFSall,PsiRBf,Mintf,SNAPdisp,PsiDEFf,MintfINV,PhiRB,Mdom)  ; ;

disp(['Number of deform/self-equilibrated modes (final) : ',num2str(size(PsiDEFf,2))])
disp('------------------------')

% if  DATAinputLOC.GetBasisMatricesFromOldCode ==1
%     disp('Retrieving basis matrix from old code...')
%     % for testing purposes
%     for e = 1:length(PhiDEF)
%         PhiDEF{e} = WSVDT(DATAROM{e}.BasisUdef,Mdom{e}) ;
%         f = [faceDOFS{e}{1};faceDOFS{e}{2}] ;
%         PsiDEFf{e} = WSVDT(DATAROM{e}.BasisRdef(f,:),cell2mat(MdomFFinv{e})) ;
%     end
%
% elseif DATAinputLOC.GetBasisMatricesFromOldCode ==2
% %     % Comparing subspaces
% %     disp('')
% %      for e = 1:length(PhiDEF)
% %          TOLcos_intersection = 0.01;
% %     [COSINE_ANGLES,WdefCAND]= PRANGLES(PhiDEF{e},PhiDEF_old{e},Mdom{e},TOLcos_intersection) ;
% %      end
% end



% PLOTTING MODES
 % --------------------------------------------------------------------- 
PlotModesDEF_SubdomainLevel(DATA,[PhiRB,PhiDEF],MESH);
  
% -------------------------------------------------
% POST-PROCESSING SELF-EQUILIBRATED MODES
% -------------------------------------------------
 PlotModesSE_SubdomainLevel(DATA,PsiRBf,PsiDEFf,MESH) ; 


% STEP 5
% ------------
% Displacement interface modes: rigid-body/Shape functions and candidates for fluctuation
% modes
% disp('Displacement interface modes...')
% 
[Vall,MESH,DATA] = FictInterfaceDISP(PhiDEF,PsiDEFf,MESH,DATAcommon,DATA,Mintf,PhiRB,PsiRBf,DATAoffline);


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
