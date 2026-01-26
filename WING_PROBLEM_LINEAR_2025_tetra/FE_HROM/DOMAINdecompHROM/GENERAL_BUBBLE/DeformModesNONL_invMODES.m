function [BASIS_DEF,Mdom,PhiRB,MESH,DATA,PsiRBf,Mintf,DATAoffline] ...
    = DeformModesNONL_invMODES(OPERFE,MESH,DATA,...
    SNAPdisp,DATAcommon,DATAoffline,SNAPdisp_AUX)
% DETERMINATION OF  INVARIANT deformational  MODES FOR A GIVEN SUBDOMAIN
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/02_BEAMQ8_IprofE.mlx
% See archetypical input files 
%  /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/02_BEAMQ8_IprofE/Train_Q8_IprofE_L.m
% and
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/02_BEAMQ8_IprofE/DEF_Q8_Ish_EL_L.m
%% JAHO 1-April-2024, Santa Cruz de Tenerife 
% -------------------------------------------------
if nargin == 0
    load('tmp1.mat')
    %     DATAoffline.BASIC_MODES_DISPLACEMENTS_GROUPS =  {[1:12],[13:60]} ;
    % DATAoffline.TOL_SVD_BASIC_MODES_DISPLACEMENTS_groups =   [0,1e-2] ;
end

disp('------------------------------------------------------------------------------')
disp('Determination of invariant modes')
disp('------------------------------------------------------------------------------')

% STEP 1
% --------------------------------------------------------------------------
% Nodes and DOFs of the interface boundaries
[faceDOFS,PhiRB,PsiRBf,Mdom,Mintf,MESH] =   GeometricVarDOMAINS(OPERFE,MESH,DATA,DATAcommon);

% STEP 2
%%%%%%%%---------------------------
% Basic tests/ Complementary (inelastic tests) tests
% --------------------------------------
ntestsCOMPL = length(DATAoffline.AdditionalTests) ;
ntestsALL = length(SNAPdisp) ;
INDEX_BASIC_TESTS = 1:(ntestsALL-ntestsCOMPL) ;
INDEX_COMPL_TESTS = (INDEX_BASIC_TESTS(end)+1):ntestsALL;
DATAoffline.INDEX_BASIC_TESTS = INDEX_BASIC_TESTS;
DATAoffline.INDEX_COMPL_TESTS = INDEX_COMPL_TESTS;

SNAPbasic = cell2mat(SNAPdisp(INDEX_BASIC_TESTS)) ;
SNAPcompl = SNAPdisp(INDEX_COMPL_TESTS) ;

%%% STEP 2. PURGING RIGID BODY DISPLACEMENTS, basic tests. Basis PhiDEFbs
% -------------------------------------------------------------------
DATAcommon = DefaultField(DATAcommon,'INFO_BASIC_DEFORMATIONAL_MODES',[]) ;
DATAcommon.INFO_BASIC_DEFORMATIONAL_MODES = DefaultField(DATAcommon.INFO_BASIC_DEFORMATIONAL_MODES,'INDEX_FUNDAMENTAL_MODES',[]) ;
DATAoffline.INFO_BASIC_DEFORMATIONAL_MODES = DATAcommon.INFO_BASIC_DEFORMATIONAL_MODES ;
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/01_BEAMQ8_ELAST.mlx
disp('Determining basic modes (elastic range)')
% We purge here the rigid body part

% REMOVING ZERO COLUMNS
%
SNAPbasic_NORM = sum(SNAPbasic.^2,1) ;
nonzerosCOL = find(SNAPbasic_NORM~=0);
SNAPbasic = SNAPbasic(:,nonzerosCOL) ;


DATAoffline = DefaultField(DATAoffline,'THRESHOLD_RATIO_CONSECUTIVE_SINGLE_VALUES_INVARIANT_DEF_MODES',10) ; %
 


SNAPbasic = SprojDEF_operator(PhiRB,Mdom,SNAPbasic) ;
%if isempty(DATAoffline.INFO_BASIC_DEFORMATIONAL_MODES.INDEX_FUNDAMENTAL_MODES)
    
 %   DATALOC.TOL = DATAoffline.TOLSVD_DETERMINE_INVARIANT_DEF_MODES ; 
    [ PhiDEFbs,Sbs,Vbs,MdomCHOL] = WSVDT( SNAPbasic,Mdom) ;
    
    RATIO_SV = Sbs(1:end-1)./Sbs(2:end) ; 
    
  III =   find(RATIO_SV>DATAoffline.THRESHOLD_RATIO_CONSECUTIVE_SINGLE_VALUES_INVARIANT_DEF_MODES) ;
  
nINVmodes = III(1) ;
disp('***********************************************************************************************************')
disp('Basic modes')  
disp(['Total number of  deformational modes = ',num2str(length(Sbs))]) ; 
disp(['Number of  deformational modes that may be deemed INVARIANT = ',num2str(nINVmodes)]) ; 
disp(['Jump in singular values (ratio consecutive singular values) =',num2str(RATIO_SV(nINVmodes))])
disp('***********************************************************************************************************')

PhiDEFbs = PhiDEFbs(:,1:nINVmodes) ; 
Sbs = Sbs(1:nINVmodes) ; 
Vbs = Vbs(:,1:nINVmodes) ; 


 
  
    
%     
%     
%     
% else
%     % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/01_BEAMQ8_ELAST.mlx
%     % FUNDAMENTAL MODES
%     indFUND = DATAoffline.INFO_BASIC_DEFORMATIONAL_MODES.INDEX_FUNDAMENTAL_MODES ;
%     indOTHERS = setdiff(1:size(SNAPbasic,2),indFUND) ;
%     % indFUND are the indices of the snapshots corresponding to
%     % "fundamental" elastic modes (i.e., 6 beam modes for a Q8 element )
%     SNAPfundamental = SNAPbasic(:,indFUND)  ;
%     % Next we apply the wSVD
%     DATAsvdLOC.TOL = 1e-10 ;
%     [PhiDEFbs_fund,SS,VV,MdomCHOL] = WSVDT(SNAPfundamental,Mdom,DATAsvdLOC) ;
%     
%     if size(PhiDEFbs_fund,2) ~= DATAoffline.INFO_BASIC_DEFORMATIONAL_MODES.NUMBER_OF_FUNDAMENTAL_MODES
%         error(['The number of fundamental elastic modes should be equal to ',num2str(DATAoffline.INFO_BASIC_DEFORMATIONAL_MODES.NUMBER_OF_FUNDAMENTAL_MODES)])
%     end
%     % Now we address the case of the other deformational modes
%     % First we purge both the components corresponding to   fundamental deformational modes
%     SNAPadditional= SprojDEF_operator([PhiDEFbs_fund],Mdom, SNAPbasic(:,indOTHERS)) ;
%     
%     DATAsvdLOC.TOL = 1e-10 ;
%     DATAsvdLOC.Mchol = MdomCHOL;
%     [PhiDEFbs_add,SS,VV] = WSVDT(SNAPadditional,[],DATAsvdLOC) ;
%     
%     PhiDEFbs = [PhiDEFbs_fund,PhiDEFbs_add] ;
%     
%     
%     if size(PhiDEFbs,2) ~= DATAoffline.INFO_BASIC_DEFORMATIONAL_MODES.NUMBER_OF_BASIC_DEFORMATIONAL_MODES
%         error(['The number of   elastic modes should be equal to ',num2str(DATAoffline.INFO_BASIC_DEFORMATIONAL_MODES.NUMBER_OF_BASIC_DEFORMATIONAL_MODES)])
%     end
%     
%     
%     
% end

% PLOTTING MODES
% ---------------------------------------------------------------------
disp('Plotting invariant modes')
DATA.NAME_BASE = 'Invariant';
PlotModesDEF_SubdomainLevel(DATA,PhiDEFbs,MESH);


 disp(['Storing information ....PhiDEFbs, Sbs and Vbs....'])
 disp(['------------------------------------------------------'])
[~, nameWSloc,~] = fileparts(DATA.NameFileMeshDATA); 
 
nameWSloc = ['MODES/','INVARIANTS_',nameWSloc,'.mat'] ; 
save(nameWSloc,'PhiDEFbs','Sbs','Vbs') ; 

disp(['Run again the program, setting DATAoffline.METHOD_INVARIANT_DEF_MODES.STAGE_TRAINING =2'])
disp(['Set also:  DATAlauch.TRAINING =1;'])
return 




% 
% %%% STEP 3. PURGING RIGID BODY DISPLACEMENTS, complementary tests. Basis
% %%% PhiDEFcomp
% 
% disp('-----------------------------------------------')
% disp('Determining complementary/inelastic modes ')
% disp('-------------------------------------------------')
% 
% 
% % METHOD_DETERMINING_MODES = 'Orthogonal_interface' ;
% %
% % METHOD_DETERMINING_MODES = 'Standard_Orthogonal_Entire_Domain' ;
% %METHOD_DETERMINING_MODES = 'Orthogonal_interface_mixed' ;
% %METHOD_DETERMINING_MODES = 'All_bubble_modes' ;
% 
% 
% % Purging rigid-body modes
% for iproj = 1:length(SNAPcompl)
%     SNAPcompl{iproj} = SprojDEF_operator(PhiRB,Mdom,SNAPcompl{iproj}) ;
% end
% 
% 
% 
% % Partitioned SVD
% nproj = 1 + length(SNAPcompl) ;
% A  = cell(1,nproj) ;
% A{1} = MdomCHOL*SNAPbasic ;
% for iproj = 1:length(SNAPcompl)
%     A{iproj+1} =  MdomCHOL*SNAPcompl{iproj} ;
% end
% TOL = DATAoffline.TOLSVD_complementary_modes*ones(size(SNAPcompl)) ;
% RELTOL = [1e-10,TOL] ;
% DATAlocS = [] ;
% [U,S_upsilonDEF,V_upsilonDEF] = SRSVD(A,RELTOL,DATAlocS) ;
% 
% UpsilonDEF = MdomCHOL\U ;  % All modes (including basic ones)
% 
% 
% DATAoffline = DefaultField(DATAoffline,'NUMBER_OF_DEFORMATIONAL_MODES',[]) ;
% 
% if  ~isempty(DATAoffline.NUMBER_OF_DEFORMATIONAL_MODES)
%     nmodes = length(S_upsilonDEF) ;
%     nmodes = min(nmodes,DATAoffline.NUMBER_OF_DEFORMATIONAL_MODES) ;
% else
%     nmodes = length(S_upsilonDEF) ;
% end
% 
% 
% BASIS_DEF.UpsilonDEF = UpsilonDEF(:,1:nmodes) ;
% BASIS_DEF.S_upsilonDEF = S_upsilonDEF(1:nmodes) ;
% BASIS_DEF.V_upsilonDEF = V_upsilonDEF(:,1:nmodes) ;
% 
% 
% 
% DATA.NAME_BASE = 'AllDEFORMATIONAL';
% PlotModesDEF_SubdomainLevel(DATA,UpsilonDEF,MESH);
% 
% % %% OUT OF CURIOSITY, HOW DOES THE SVD OF THE BOUNDARY ENTRIES LOOK LIKE ?
% %                 b = MESH.faceDOFSall;
% %
% % Ub  = bsxfun(@times,UpsilonDEF(b,:)',S_upsilonDEF)' ;
% % DATALOCaaa.TOL = 1e-4 ;
% % [ Ub_b,S_b,V_b] = WSVDT( Ub,Mintf,DATALOCaaa) ;
% %
% 
% 
% % ORTHOGONAL COMPLEMENT (visualization purposes)
% 
% 
% ncomp  = size(BASIS_DEF.UpsilonDEF,2)-size(PhiDEFbs,2) ;
% 
% if ncomp >1 && ~isempty(SNAPcompl)
%     
%     PhiDEFcomp = SprojDEF_operator(PhiDEFbs,Mdom, BASIS_DEF.UpsilonDEF) ;
%     DATALOC.TOL = 1e-6 ;
%     DATALOC.Mchol = MdomCHOL ;
%     [ PhiDEFcomp,S,~,~] = WSVDT( PhiDEFcomp,[],DATALOC) ;
%     
%     PhiDEFcomp = PhiDEFcomp(:,1:ncomp)  ;
%     
%     DATA.NAME_BASE = 'Complementary_orth';
%     PlotModesDEF_SubdomainLevel(DATA,PhiDEFcomp,MESH);
%     
%     
% else
%     disp(['There are no complementary (inelastic modes ) here, only elastic'])
%     PhiDEFcomp = [] ;
% end
% 
% 
% 
% BASIS_DEF.PhiDEFbs = PhiDEFbs ;
% BASIS_DEF.PhiDEFcomp = PhiDEFcomp ;
% 
% 
