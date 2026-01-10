function [PhiDEF,PsiDEFf,Vall,Mintf,Mdom,PhiRB,PsiRBf,MESH,DATA] ...
    = MODESandROMoper_elastic(OPERFE,MESH,DATA,SNAPreact,BasisStwo,...
    SNAPdisp,DATAcommon,DATAoffline,MATPRO)
% DETERMINATION OF MODES FOR A GIVEN SUBDOMAIN
% 1) Elastic range
% 2) Trained without external forces
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/103_EIFEM_PLATES/02_CONV_SIZE_TRAC.mlx
% JAHO, 23-AUG-2023
% ----------------------------

if nargin == 0
    load('tmp.mat')
    
end
% ------------------------------------------------------------------------------------------------
% STEP 1
% -------
% Nodes and DOFs of the interface boundaries
[faceDOFS,PhiRB,PsiRBf,Mdom,Mintf,MESH] =   GeometricVarDOMAINS(OPERFE,MESH,DATA,DATAcommon);

if ~isempty(MATPRO.celasglo)
    %     disp('Checking accuracy CECM points ')
    %     disp('-------------------------------------------')
    %   disp('Coarse-scale stiffness matrix ')
    celastST = MATPRO.celasglo ;
    nF = size(MATPRO.celasglo,2) ;
    for icomp = 1:nF
        icol = icomp:nF:size(celastST,1) ;
        celastST(icol,:) = bsxfun(@times,celastST(icol,:),OPERFE.wSTs) ;
    end
    celastST = ConvertBlockDiag(celastST) ; % Diagonal block matrix
    Kstiff = OPERFE.Bst'*(celastST*OPERFE.Bst);
else
    error('Option not viable')
end

% Deformational modes
disp('Deformational modes (COMPUTED DIRECTLY)...')

DATAoffline= DefaultField(DATAoffline,'errorSVD_DEFORMATIONAL_MODES',0) ; %errorSVD_DEFORMATIONAL_MODES=  1e-4;
AsnapDISPdef =  SprojDEF_operator(PhiRB,Mdom,SNAPdisp) ;
DATALOC.TOL = DATAoffline.errorSVD_DEFORMATIONAL_MODES ;

DATAoffline = DefaultField(DATAoffline,'NumberOfDeformationalModes',[] ) ; 

if ~isempty(DATAoffline.NumberOfDeformationalModes)
    DATALOC.TOL = 0 ; 
end 



[PhiDEF,S,V,~ ]= WSVDT(AsnapDISPdef,Mdom,DATALOC) ;

if ~isempty(DATAoffline.NumberOfDeformationalModes)
   nmodes = min(size(PhiDEF,2),DATAoffline.NumberOfDeformationalModes) ; 
   PhiDEF = PhiDEF(:,1:nmodes) ; 

end 


%%%% 
% Uncomment the following lines to check that the approach is indeed valid (already line) 
% [AsnapREACse,MdomFFinv] =  HprojSEf_operator(PsiRBf,Mintf,SNAPreact(MESH.faceDOFSall,:)) ;
% [PsiDEFf_fromSNAP,S,V,MdomFFinv_chol ]= WSVDT(AsnapREACse,MdomFFinv,DATALOC) ;

MintfINV = inv(Mintf) ; 


PsiDEFf = Kstiff*PhiDEF ; 
PsiDEFf = PsiDEFf(MESH.faceDOFSall,:) ; 
[PsiDEFf,S,V,~ ]= WSVDT(PsiDEFf,MintfINV) ;
%  
% [P1,~,~] =SVDT(PsiDEFf_fromSNAP)  ; 
% [P2,~,~] =SVDT(PsiDEFf)  ; 
% [UU,SS,VV] = SVDT(P1'*P2) ; 



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
