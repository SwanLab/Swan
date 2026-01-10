function [BASIS_DEF,Mdom,PhiRB,MESH,DATA,PsiRBf,Mintf,DATAoffline] ...
    = DeformModesNONL_invCOMPALL(OPERFE,MESH,DATA,...
    SNAPdisp,DATAcommon,DATAoffline,SNAPdisp_AUX)
%
%% JAHO, 3-Apr-2024
%  Barcelona, Sandwichez, C/Mallorca
% See  /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/02_BEAMQ8_IprofE.mlx
% Determination of deformational modes for the case in which we are
% interesed in computing "invariant" modes, such as "beam" modes
% See training stage for elastic range
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/02_BEAMQ8_IprofE/Train_Q8_48tests.m

if nargin == 0
    load('tmp1.mat')
end

% STEP 1
% --------------------------------------------------------------------------
% Nodes and DOFs of the interface boundaries
%[faceDOFS,PhiRB,PsiRBf,Mdom,Mintf,MESH] =   GeometricVarDOMAINS(OPERFE,MESH,DATA,DATAcommon);

% STEP 2
%%%%%%%%---------------------------
% Basic tests/ Complementary (inelastic tests) tests
% --------------------------------------
ntestsCOMPL = length(DATAoffline.AdditionalTests) ;  % THESE ARE THE INELASTIC TESTS
ntestsALL = length(SNAPdisp) ;
INDEX_BASIC_TESTS = 1:(ntestsALL-ntestsCOMPL) ;   % THESE ARE THE INDICES OF THE ELASTIC MDES
INDEX_COMPL_TESTS = (INDEX_BASIC_TESTS(end)+1):ntestsALL;
DATAoffline.INDEX_BASIC_TESTS = INDEX_BASIC_TESTS;
DATAoffline.INDEX_COMPL_TESTS = INDEX_COMPL_TESTS;

SNAPbasic = cell2mat(SNAPdisp(INDEX_BASIC_TESTS)) ;  % ELASTIC MODES
SNAPcompl = SNAPdisp(INDEX_COMPL_TESTS) ;

% REMOVING ZERO COLUMNS, FOR ELASTIC TESTS
%
SNAPbasic_NORM = sum(SNAPbasic.^2,1) ;
nonzerosCOL = find(SNAPbasic_NORM~=0);
SNAPbasic = SNAPbasic(:,nonzerosCOL) ;

%%% STEP 2.  INVARIANT MODES
% -------------------------------------------------------------------
DATAcommon = DefaultField(DATAcommon,'INFO_BASIC_DEFORMATIONAL_MODES',[]) ;

if ~isfield(DATAcommon.INFO_BASIC_DEFORMATIONAL_MODES,'INDICES_for_invariant_TRAINING_TESTS')
    
    error('You have not specified which are the indices of the tests in which we have to search for the invariant modes')
    %  See for instance /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/LARGE_STRAINS/InputDataFunctions/Q8_for_beamBEHAVIOR_48.m
    
    
end

SNAPinv = SNAPbasic(:,DATAcommon.INFO_BASIC_DEFORMATIONAL_MODES.INDICES_for_invariant_TRAINING_TESTS) ;


% Purging rigid-body displacements
SNAPinv = SprojDEF_operator(PhiRB,Mdom,SNAPbasic(:,DATAcommon.INFO_BASIC_DEFORMATIONAL_MODES.INDICES_for_invariant_TRAINING_TESTS)) ;

DATAlocc.TOL = 1e-10;
[ PhiDEFbs,Sbs,Vbs,MdomCHOL] = WSVDT( SNAPinv,Mdom,DATAlocc) ;

RATIO_SV = Sbs(1:end-1)./Sbs(2:end) ;

DATAoffline = DefaultField(DATAoffline,'THRESHOLD_RATIO_CONSECUTIVE_SINGLE_VALUES_INVARIANT_DEF_MODES',10) ;

III =   find(RATIO_SV>DATAoffline.THRESHOLD_RATIO_CONSECUTIVE_SINGLE_VALUES_INVARIANT_DEF_MODES) ;

nINVmodes = III(1) ;
disp('***********************************************************************************************************')
disp('Searching for invariant modes')
disp(['Total number of  deformational modes (for the first set of training tests) = ',num2str(length(Sbs))]) ;
disp(['Number of  deformational modes that may be deemed INVARIANT = ',num2str(nINVmodes)]) ;
disp(['Jump in singular values (ratio consecutive singular values) =',num2str(RATIO_SV(nINVmodes))])
disp('***********************************************************************************************************')

PhiDEFbs_INV = PhiDEFbs(:,1:nINVmodes) ;

disp(['*****************************************************************************************'])

if ~isfield(DATAcommon.INFO_BASIC_DEFORMATIONAL_MODES,'NUMBER_OF_DEFORMATIONAL_MODES')    
    error('You have not specified THE TOTAL NUMBER OF DEFORMATIONAL MODES FOR THIS TYPE OF ELEMENT')
    %  See for instance /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/LARGE_STRAINS/InputDataFunctions/Q8_for_beamBEHAVIOR_48.m 
end

nDEF_basic = DATAcommon.INFO_BASIC_DEFORMATIONAL_MODES.NUMBER_OF_DEFORMATIONAL_MODES;

if nDEF_basic == nINVmodes
    disp(['No need for determining complementary deformational modes for elastic behavior'])
    PhiDEFbs = PhiDEFbs_INV ; 
else
    disp(['Computing complementary basic tests (aside from those regarded as "invariant")    '])
    nDEF_compb = nDEF_basic-nINVmodes; 
    disp(['Number of complementary deformational modes =',num2str(nDEF_compb)])
end

indCOMPLbasic = setdiff(1:size(SNAPbasic,2),DATAcommon.INFO_BASIC_DEFORMATIONAL_MODES.INDICES_for_invariant_TRAINING_TESTS) ; 
% Purging rigid-body displacements
SNAPbasicCOMPL = SprojDEF_operator(PhiRB,Mdom, SNAPbasic(:,indCOMPLbasic)) ;
DATAlocc.TOL = 1e-10;
[ PhiDEFcomple,Sbs,Vbs,MdomCHOL] = WSVDT( SNAPbasicCOMPL,Mdom,DATAlocc) ;

% 





% 





% PLOTTING MODES
% ---------------------------------------------------------------------
DATA.NAME_BASE = 'Basic';
PlotModesDEF_SubdomainLevel(DATA,PhiDEFbs,MESH);

%%% STEP 3. PURGING RIGID BODY DISPLACEMENTS, complementary tests. Basis
%%% PhiDEFcomp

disp('-----------------------------------------------')
disp('Determining complementary/inelastic modes ')
disp('-------------------------------------------------')


% METHOD_DETERMINING_MODES = 'Orthogonal_interface' ;
%
% METHOD_DETERMINING_MODES = 'Standard_Orthogonal_Entire_Domain' ;
%METHOD_DETERMINING_MODES = 'Orthogonal_interface_mixed' ;
%METHOD_DETERMINING_MODES = 'All_bubble_modes' ;


% Purging rigid-body modes
for iproj = 1:length(SNAPcompl)
    SNAPcompl{iproj} = SprojDEF_operator(PhiRB,Mdom,SNAPcompl{iproj}) ;
end



% Partitioned SVD
nproj = 1 + length(SNAPcompl) ;
A  = cell(1,nproj) ;
A{1} = MdomCHOL*SNAPbasic ;
for iproj = 1:length(SNAPcompl)
    A{iproj+1} =  MdomCHOL*SNAPcompl{iproj} ;
end
TOL = DATAoffline.TOLSVD_complementary_modes*ones(size(SNAPcompl)) ;
RELTOL = [1e-10,TOL] ;
DATAlocS = [] ;
[U,S_upsilonDEF,V_upsilonDEF] = SRSVD(A,RELTOL,DATAlocS) ;

UpsilonDEF = MdomCHOL\U ;  % All modes (including basic ones)


DATAoffline = DefaultField(DATAoffline,'NUMBER_OF_DEFORMATIONAL_MODES',[]) ;

if  ~isempty(DATAoffline.NUMBER_OF_DEFORMATIONAL_MODES)
    nmodes = length(S_upsilonDEF) ;
    nmodes = min(nmodes,DATAoffline.NUMBER_OF_DEFORMATIONAL_MODES) ;
else
    nmodes = length(S_upsilonDEF) ;
end


BASIS_DEF.UpsilonDEF = UpsilonDEF(:,1:nmodes) ;
BASIS_DEF.S_upsilonDEF = S_upsilonDEF(1:nmodes) ;
BASIS_DEF.V_upsilonDEF = V_upsilonDEF(:,1:nmodes) ;



DATA.NAME_BASE = 'AllDEFORMATIONAL';
PlotModesDEF_SubdomainLevel(DATA,UpsilonDEF,MESH);

% %% OUT OF CURIOSITY, HOW DOES THE SVD OF THE BOUNDARY ENTRIES LOOK LIKE ?
%                 b = MESH.faceDOFSall;
%
% Ub  = bsxfun(@times,UpsilonDEF(b,:)',S_upsilonDEF)' ;
% DATALOCaaa.TOL = 1e-4 ;
% [ Ub_b,S_b,V_b] = WSVDT( Ub,Mintf,DATALOCaaa) ;
%


% ORTHOGONAL COMPLEMENT (visualization purposes)


ncomp  = size(BASIS_DEF.UpsilonDEF,2)-size(PhiDEFbs,2) ;

if ncomp >1 && ~isempty(SNAPcompl)
    
    PhiDEFcomp = SprojDEF_operator(PhiDEFbs,Mdom, BASIS_DEF.UpsilonDEF) ;
    DATALOC.TOL = 1e-6 ;
    DATALOC.Mchol = MdomCHOL ;
    [ PhiDEFcomp,S,~,~] = WSVDT( PhiDEFcomp,[],DATALOC) ;
    
    PhiDEFcomp = PhiDEFcomp(:,1:ncomp)  ;
    
    DATA.NAME_BASE = 'Complementary_orth';
    PlotModesDEF_SubdomainLevel(DATA,PhiDEFcomp,MESH);
    
    
else
    disp(['There are no complementary (inelastic modes ) here, only elastic'])
    PhiDEFcomp = [] ;
end



BASIS_DEF.PhiDEFbs = PhiDEFbs ;
BASIS_DEF.PhiDEFcomp = PhiDEFcomp ;


