function [BASIS_DEF,MESH,DATA,DATAoffline,Kstiff] ...
    = DeformModesNONL_INVARIANTSC(OPERFE,MESH,DATA,...
    SNAPdisp,DATAcommon,DATAoffline,Vintf,PhiRB,Mdom,Mintf,MATPRO)
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

% Stiffness matrix

Kstiff = StiffMatrixRecover(MATPRO,OPERFE)  ;


% step 3) INVARIANT MODES
% -------------------------
[PhiDEFbs_INV,MdomCHOL,SNAPbasic_COMPLEM] = InvariantModes_basic(SNAPbasic,DATAcommon,PhiRB,Mdom,DATAoffline,MESH,DATA,Kstiff) ;

% step 4) COMPLEMENTARY BASIC MODES --> bASIC MODES


METHOD =  'NO_COMPLEMENTARY_TESTS' ; %'BASED_ON_FORCES'

switch  METHOD
    case 'NO_COMPLEMENTARY_TESTS'
        [PhiDEFbs] = ComplementaryModes_NOTESTS...
            (PhiDEFbs_INV,DATAcommon,SNAPbasic,MdomCHOL,MESH,PhiRB,Vintf,Mintf,DATA,MATPRO,OPERFE,Mdom,Kstiff) ;
    case 'BASED_ON_FORCES'
        %
        [PhiDEFbs] = ComplementaryModes_basicF...
            (PhiDEFbs_INV,DATAcommon,SNAPbasic,MdomCHOL,MESH,PhiRB,Vintf,Mintf,DATA,MATPRO,OPERFE,Mdom,SNAPbasic_COMPLEM,Kstiff) ;
    otherwise
        % MEthod before Apr-9th-2024
        [PhiDEFbs] = ComplementaryModes_basic...
            (PhiDEFbs_INV,DATAcommon,SNAPbasic,MdomCHOL,MESH,PhiRB,Vintf,Mintf,DATA,MATPRO,OPERFE,Mdom,SNAPbasic_COMPLEM,Kstiff) ;
end



%%% STEP 5. PURGING RIGID BODY DISPLACEMENTS, complementary tests. Basis
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
A{1} = MdomCHOL*PhiDEFbs ;
for iproj = 1:length(SNAPcompl)
    A{iproj+1} =  MdomCHOL*SNAPcompl{iproj} ;
end
TOL = DATAoffline.TOLSVD_complementary_modes*ones(size(SNAPcompl)) ;
RELTOL = [0,TOL] ;
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




