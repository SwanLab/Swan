function [PhiDEFbs,MintfINV_chol,DATA] = ComplementaryModes_NOTESTS(PhiDEFbs_all,DATAcommon,SNAPbasic,MdomCHOL,...
    MESH,PhiRB,Vintf,Mintf,DATA,MATPRO,OPERFE,Mdom,Kstiff,DATAoffline)
% Determination of complementary basic modes (complementary to invariant modes (PhiDEFbs_all), which are determined in
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/DOMAINdecompHROM/GENERAL_BUBBLE/InvariantModes_basic.m)
% The actual output is PhiDEFbs, which are the set of basic, elastic modes.
% SNAPbasic --> All basic snapshots (elastic)
% Vintf ---> Fictitious interface modes
% PhiRB --> Rigid body modes, domain
% Modification of ComplementaryModes_basicF.m
% No complementary tests are required here
% COMMENTS BY CHATGPT  (6-MAY-2025)


% =========================================================================
% FUNCTION: ComplementaryModes_NOTESTS
% =========================================================================
% PURPOSE:
% Constructs the final set of basic deformation modes (`PhiDEFbs`) 
% by identifying and separating:
%   - Invariant modes: those aligned with reactive forces at the interface
%   - Complementary modes: additional modes orthogonal to invariants
%
% This version does **not** require complementary test simulations.
%
% CONTEXT:
% - Used in domain decomposition-based reduced order modeling (FE-HROM)
% - Builds on deformation and interface modes for static condensation
% - Adaptation of `ComplementaryModes_basicF.m` but without test problems
%
% INPUTS:
% - PhiDEFbs_all      : Full set of candidate deformation modes
% - DATAcommon        : Common data structure (includes element type info)
% - SNAPbasic         : Basic displacement snapshots (e.g., from offline stage)
% - MdomCHOL          : Cholesky decomposition of domain mass matrix (optional)
% - MESH              : FE mesh info (face DOFs, connectivity, etc.)
% - PhiRB             : Rigid body modes
% - Vintf             : Fictitious interface deformation modes
% - Mintf             : Interface mass matrix
% - DATA              : General simulation data
% - MATPRO            : Material properties
% - OPERFE            : Operators for FE assembly
% - Mdom              : Mass matrix on the domain
% - Kstiff            : Global stiffness matrix
% - DATAoffline       : Configuration flags (e.g., SVD settings, tolerances)
%
% OUTPUTS:
% - PhiDEFbs          : Final basic elastic deformation modes (invariant + complementary)
% - MintfINV_chol     : Cholesky factor of inverse interface mass matrix (for future use)
% - DATA              : Updated with indices and references to selected modes
%
% KEY STEPS:
% 1. Determine how many basic deformation modes are needed (based on element type).
% 2. Identify invariant modes by computing their alignment with interface reactions.
% 3. Orthogonalize both the invariant modes and the interface deformation modes.
% 4. Compute complementary modes if invariant modes are insufficient.
% 5. Ensure orthogonality of complementary modes with respect to rigid body and invariant modes.
% 6. Assemble the final deformation basis matrix `PhiDEFbs` and verify independence.
%
% VALIDATION:
% - Principal angles between mode sets are used to assess linear independence.
% - Cosines of angles near 0 may indicate near-dependence and potential ill-posedness.
% - Visual inspection of modes is recommended via `PlotModesDEF_SubdomainLevel` and `PlotModesGID_Interf`.
%
% FILE REFERENCES:
% - InvariantModes_basic.m
% - Test case: /TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/02_BEAMQ8_IprofE.mlx
% =========================================================================


% 
% JAHO, 10-Apr-2024, Balmes 185, Barcelona
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/02_BEAMQ8_IprofE.mlx
% -------------------------------
if nargin == 0
    load('tmp.mat')
end

%
DATAcommon = DefaultField(DATAcommon,'INFO_BASIC_DEFORMATIONAL_MODES',[]) ;


if ~isfield(DATAcommon.INFO_BASIC_DEFORMATIONAL_MODES,'NUMBER_OF_DEFORMATIONAL_MODES')
    % disp('You have not specified THE TOTAL NUMBER OF basic DEFORMATIONAL MODES FOR THIS TYPE OF ELEMENT')
    %disp(['Retrieving this number from the type of EIF element '])
    
    switch DATAcommon.TypeFunctionDisplacementInterfaces
        
        case '1D_LINEAR'
            nDEF_basic = 1;
            
        case 'HEXAHEDRA_QUADRATIC'
            nDEF_basic = 72 ; % 3*26-6
            
            
        case 'HEXAHEDRA_LINEAR'
            nDEF_basic = 18 ; % 3*8-6
            
        case 'QUADRILATERAL_QUADRATIC'
            nDEF_basic = 13 ; % 2*8-3
            
        case 'QUADRILATERAL_LINEAR'
            nDEF_basic = 5 ; % 2*4-3
            
        case 'QUADRILATERAL_LINEAR_UNCOUPLED'
            nDEF_basic = 13 ; % 2*8-3
        otherwise
            
            nDEF_basic = size(Vintf,2)-size(PhiRB,2);
            
            
    end
    
    disp(['Number of intef. def modes = ',num2str(nDEF_basic)])
    
else
    nDEF_basic = DATAcommon.INFO_BASIC_DEFORMATIONAL_MODES.NUMBER_OF_DEFORMATIONAL_MODES;
    
    %  See for instance /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/LARGE_STRAINS/InputDataFunctions/Q8_for_beamBEHAVIOR_48.m
end
%
%n_inv = size(PhiDEFbs_all,2) ;

%if nDEF_basic == n_inv
%   disp(['No need for determining complementary deformational modes for elastic behavior'])
%  PhiDEFbs = PhiDEFbs_all ;
%else

% disp('--------------------------------------------------------------------------------')

% STEP 1) WE BEGIN BY DETERMINING THE DEFORMATIONAL PART OF THE FICTITIOUS
% INTERFACE MODES
[VintfDEF,COORbnd,CNbREN,b,Mintf_chol] = DefPartInterfaceModes(MESH,PhiRB,Vintf,Mintf,DATA,nDEF_basic,DATAoffline) ;

% -------------------------------------------------------------------------------------------------------------
%  STEP 2) DETERMINE THE nINV modes which are more aligned to the reactive modes caused by the invariant modes
% --------------------------------------------------------------------------------------------------------------
b = MESH.faceDOFSall;
% Reactive forces at interface nodes
X_lambda = Kstiff(b,:)*PhiDEFbs_all ;  %

% DATA.NAME_BASE = 'prueba2';
% PlotModesDEF_SubdomainLevel(DATA,PhiDEFbs_all,MESH);
%
%
%  PlotModesGID_Interf(COORbnd,CNbREN,MESH,-X_lambda,DATA,'prueba3') ;




if DATAoffline.USE_CHOLESKY_DECOMPOSITION == 1
    
    MintfINV = inv(Mintf) ;
    [ PsiSE_all,S_lambda,V_lambda,MintfINV_chol]= WSVDT(X_lambda,MintfINV) ; % (orthogonaliz.)
    % ------------------
    % Deformational fict. modes (orthogoComplementaryModes_NOTESTSnaliz.)
    DATAkk.Mchol = Mintf_chol ;
    [VintfDEF,~,~,Mintf_chol ]= WSVDT(VintfDEF,[],DATAkk) ; % These are the deformational fictitious interface modes (orhogonalized)
    
else
    MintfINV = [] ;
    [ PsiSE_all,S_lambda,V_lambda,MintfINV_chol]= WSVDT(X_lambda,MintfINV) ; % (orthogonaliz.)
    % ------------------
    % Deformational fict. modes (orthogoComplementaryModes_NOTESTSnaliz.)
    %  DATAkk.Mchol = Mintf_chol ;
    [VintfDEF,~,~,Mintf_chol ]= WSVDT(VintfDEF,[]) ; % These are the deformational fictitious interface modes (orhogonalized)
    
end

if size(VintfDEF,2) ~=nDEF_basic
    
    disp(['Number of deformational intef. modes = ',num2str(size(VintfDEF,2))])
    disp(['Number of "basic" deformational modes for this type of EIF element = ',num2str(nDEF_basic)])
    
    disp(['Chech how the intf. modes are constructed; if the element is curved, introduce the inverse mapping option'])
    disp(['Set [].ORDER_INVERSE_ELEMENT_inverse_mapping to a higher order' ])
    
    error('')
    
end

% PRINCIPAL ANGLES AND PRINCIPAL VECTORS
[U_v,S_v,V_v] = SVDT(PsiSE_all'*VintfDEF) ;


disp(['Cosine angles  formed by INVARIANT fict. def. interface modes and  invariant reactive forces'])
S_v

% MODIFICATION INTRODUCED IN
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/03_PLATE_Q8.mlx
DATAoffline = DefaultField(DATAoffline,'TOL_cosine_ANGLES_WORK_DONE_INTERFACES',1e-4) ;
IIIA = find(S_v < DATAoffline.TOL_cosine_ANGLES_WORK_DONE_INTERFACES) ;
% See also
% /home/joaquin/Desktop/CURRENT_TASKS/PAPERS_2020_onwards/12_EIFEM_EXTENS/EIFEM_bub.pdf
if ~isempty(IIIA)
    n_inv = IIIA(1)-1 ;
    disp('-----------------------------------------------------------------------------------------------------------------')
    disp(['There are ',num2str(length(IIIA)),' reactive modes that do not contribute to the work done by the fict. interfaces'])
    disp(['(because the cosine of their principal angles  is below ',num2str(DATAoffline.TOL_cosine_ANGLES_WORK_DONE_INTERFACES),'  '])
    disp(['Accordingly, we shall only keep ',num2str(n_inv),' invariant modes'])
    disp('It should be noted that this will necessary affect the accuracy in reproducing each of the training test')
    disp('-----------------------------')
    
  
    
    
else
    % PhiDEFbs_inv = PhiDEFbs_all ;
    n_inv = length(S_v) ;
    disp(['Number of invariant modes = ',num2str(n_inv)])
end



% REACTION MODES MORE ALIGNED WITH THE INTERFACE MODES
PsiSE_inv = PsiSE_all*U_v(:,1:n_inv) ;
% DEFORMATIONAL MODES ASSOCIATED TO THE ABOVE MODES
%  \PhiDEF^{inv}  =   \PhiDEF^{all}  \V_{\lambda} \diag{\S_{\lambda}}^{-1}  \U_v(:,1:n_{inv})

% include_SingularVAlues =0 ;
%
% if include_SingularVAlues == 0
%PhiDEFbs_inv = PhiDEFbs_all*V_lambda*U_v(:,1:n_inv) ;
%else
% AFTER 10-mARCH-2025, sEE /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/07_WING.mlx
PhiDEFbs_inv = PhiDEFbs_all*V_lambda ;
PhiDEFbs_inv = bsxfun(@times,PhiDEFbs_inv',1./S_lambda)' ;
PhiDEFbs_inv = PhiDEFbs_inv*U_v(:,1:n_inv) ;
PLOT_MORE_ALIGNED_MODES= 1;
if PLOT_MORE_ALIGNED_MODES == 1
    DATA.NAME_BASE = 'AlignedDEFmodes';
    PlotModesDEF_SubdomainLevel(DATA,PhiDEFbs_inv,MESH);
    
end


% ORTHOGONALIZATION
DATAlocc.Mchol = MdomCHOL ;
[ PhiDEFbs_inv,SSS_inv,VVV_inv,~] = WSVDT( PhiDEFbs_inv,[],DATAlocc) ;
%end





%NameLoc = 'prueba_reactions';
%PlotModesGID_Interf(COORbnd,CNbREN,MESH,PsiSE_inv,DATA,NameLoc) ;
% REACT = Kstiff(b,:)*PhiDEFbs_inv  ;
% PlotModesGID_Interf(COORbnd,CNbREN,MESH,REACT,DATA,NameLoc) ;

% Invariant component of the deformational fict. modes aligned with PsiSE_inv
VintDEF_INV = VintfDEF*V_v(:,1:n_inv) ;
if PLOT_MORE_ALIGNED_MODES == 1
    NameLoc= 'AlignedDEFmodes_interf';
    PlotModesGID_Interf(COORbnd,CNbREN,MESH,VintDEF_INV,DATA,NameLoc) ;
end


%
% PRUEBAS  = 1;
% if PRUEBAS == 1
%     disp('BORRAR ESTO')
%     DATAkk.Mchol = MintfINV_chol ;
%     [X_lambda_1,S_lambda_1,V_lambda_1]= WSVDT(X_lambda(:,1),[],DATAkk) ; % These are the deformational fictitious interface modes (orhogonalized)
%
%     PsiSE_inv_1 = PsiSE_inv(:,1)  ;
%      VintDEF_INV_1 = VintDEF_INV(:,1) ;
%     [UUU,SSS,VVV] = SVDT(X_lambda_1'*VintDEF_INV_1) ;
%
%     [UUU2,SSS2,VVV2] = SVDT(PsiSE_inv_1'*VintDEF_INV_1) ;
%
%     AUX1 = PhiDEFbs_all*V_lambda ;
%     AUX1 = bsxfun(@times,AUX1',1./S_lambda)' ;
%     PhiDEFbs_inv_new = AUX1*U_v(:,1:n_inv) ;
%
%      DATA.NAME_BASE = 'prueba2';
%      PlotModesDEF_SubdomainLevel(DATA,PhiDEFbs_inv_new,MESH);
%     NameLoc = 'prueba_valign';
%
%   PlotModesGID_Interf(COORbnd,CNbREN,MESH,VintDEF_INV,DATA,NameLoc) ;
%
% end



if  n_inv == nDEF_basic
    
    disp(['... No need to determine complementary modes '])
    
    PhiDEFbs = PhiDEFbs_inv;
else
    
    disp('----------------------------------------------------------------------------------------------')
    disp(['Computing complementary basic modes (aside from those regarded as "invariant")    '])
    nDEF_compb = nDEF_basic-n_inv;
    disp(['Number of such complementary deformational modes =',num2str(nDEF_compb)])
    
    
    
    
    % ------------------------------------------------
    % STEP 3) COMPLEMENTARY FICT. MODES
    % ---------------------------------------------------
    VintDEF_compl = VintfDEF  - VintDEF_INV*(VintDEF_INV'*Mintf*VintfDEF) ;
    
    if DATAoffline.USE_CHOLESKY_DECOMPOSITION == 1
        DATAlocW.Mchol = Mintf_chol ;
        DATAlocW.TOL = 1e-6;
        [VintDEF_compl,SSpsv,VVpsv] = WSVDT(VintDEF_compl,[],DATAlocW) ;
        
    else
        %  DATAlocW.Mchol = Mintf_chol ;
        DATAlocW.TOL = 1e-6;
        [VintDEF_compl,SSpsv,VVpsv] = WSVDT(VintDEF_compl,[],DATAlocW) ;
        
    end
    
    disp(['Angles (in degrees) formed by COMPLEMENTARY fict. def. interface modes and  invariant reactive forces'])
    
    [UUpv,SSpv,VVpv] = SVDT(PsiSE_inv'*VintDEF_compl) ;
    
    acosd(SSpv)
    
    
    % Plotting both sets of modes
    
    disp(['Deformed shapes of the  deformational fict. interface modes: first ',num2str(size(VintDEF_INV,2)) ' most aligned with invariant modes ',[]])
    disp(['Remaining  ',num2str(size(VintDEF_compl,2)) ' most orthogonal to invariants ',[]]) ;
    
    NameLoc =     'FictDEFModAlignINVARIANTS+Complem' ;
    
    PlotModesGID_Interf(COORbnd,CNbREN,MESH,[VintDEF_INV,VintDEF_compl],DATA,NameLoc) ;
    
    
    
    %*******************************************************************************
    % STEP 4) COMPLEMENTARY MODES --> DEFORMATIONAL BASIS MATRIX
    %------------------------------------------------------------------------
    disp('Searching for the complementary modes...')
    disp('---------------------------------------------------')
    %     SNAPbasic_COMPLEM = SprojDEF_operator(PhiDEFbs_all,Mdom,SNAPbasic_COMPLEM) ;
    %     DATAlocc.TOL = 1e-10;
    %     DATAlocc.Mchol = MdomCHOL ;
    %     [ PhiDEFcompl,SbsC,VbsC,~] = WSVDT( SNAPbasic_COMPLEM,[],DATAlocc) ;
    %     nCOMP_avail = size(PhiDEFcompl,2) ;
    %     disp(['Number of available deformational modes complementary to the invariants = ',num2str(nCOMP_avail)]) ;
    %     if nDEF_compb > nCOMP_avail
    %         error('The number of available complementary modes is smaller than the number of necessary complementary modes ')
    %     end
    
    PhiDEFcompl = zeros(size(PhiDEFbs_all,1),size(VintDEF_compl,2)) ;
    
    PhiDEFcompl(b,:) = VintDEF_compl ;
    s = setdiff(1:size(PhiDEFbs_all,1),b) ; % Interior DOFs
    
    PhiDEFcompl(s,:) = Kstiff(s,s)\(-Kstiff(s,b)*PhiDEFcompl(b,:)) ;
    
    
    disp(['Deformed shapes PhiDEFcompl before orthogonalization '])
    
    DATA.NAME_BASE = 'PhiDEFcompl';
    PlotModesDEF_SubdomainLevel(DATA,PhiDEFcompl,MESH);
    
    
    
    disp(['Deformed shapes PhiDEFcompl after orthogonalization '])
    % PhiDEFcompl = SprojDEF_operator(PhiDEFbs_inv,Mdom,PhiDEFcompl) ;
    
    PhiDEFcompl = SprojDEF_operator([PhiRB,PhiDEFbs_inv],Mdom,PhiDEFcompl) ;% Change introduced
    % on April 24th 2025, see /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/10_AUXETIC_2D.mlx
    
    DATAlocc.TOL = 1e-10;
    
    if DATAoffline.USE_CHOLESKY_DECOMPOSITION == 1
        
        DATAlocc.Mchol = MdomCHOL ;
        [ PhiDEFcompl,SbsC,VbsC,~] = WSVDT( PhiDEFcompl,[],DATAlocc) ;
        
    else
        DATAlocc.Mchol = [] ;
        [ PhiDEFcompl,SbsC,VbsC,~] = WSVDT( PhiDEFcompl,[],DATAlocc) ;
    end
    
    
    nCOMP_avail = size(PhiDEFcompl,2) ;
    disp(['Number of available deformational modes complementary to the invariants = ',num2str(nCOMP_avail)]) ;
    if nDEF_compb > nCOMP_avail
        error('The number of available complementary modes is smaller than the number of necessary complementary modes ')
    end
    
    PhiDEFbs = [PhiDEFbs_inv,PhiDEFcompl] ;
    DATA.INDEXES_LINEAR_MODES_DEFORMATION.INVARIANT = [1:size(PhiDEFbs_inv,2)] ;   % 28-Feb-2025
    DATA.INDEXES_LINEAR_MODES_DEFORMATION.COMPLEMENTARY = [size(PhiDEFbs_inv,2)+1:size(PhiDEFbs,2)] ;
    
    
    
    
    
    %     % ---------------------------------------------------------------------------
    %     % STEP 5) FORCES INDUCED BY THE COMPLEMENTARY MODES (Basis matrix)
    %     % ---------------------------------------------------------------------------
    %     DATAlocK.Mchol = MintfINV_chol ;
    %     DATAlocK.TOL = 1e-10 ;
    %     PsiDEFcompl = Kstiff(b,:)*PhiDEFcompl ;
    %     [PsiDEFcompl, SSp,VVp]= WSVDT(PsiDEFcompl,[],DATAlocK) ; % (orthogonaliz.)
    %
    %     % *********************************************
    %     % STEP 6) Subspace of   PsiDEFcompl most aligned with VintDEF_compl
    %     %*******************************************************************************
    %
    %     [UUxxx,SSxx,VVxxx] = SVDT(PsiDEFcompl'*VintDEF_compl) ;
    %
    %     disp([ 'Cosine Angles formed by COMPLEMENTARY fict. def. interface modes and  COMPLEMENTARY reactive forces'])
    %     (SSxx)
    %
    %
    %     [UUxxxS,SSxxS,VVxxxS] = SVDT(PsiDEFcompl'*VintDEF_INV) ;
    %
    %     disp(['Cosine Angles   formed by INVARIANT fict. def. interface modes and  COMPLEMENTARY reactive forces'])
    %     (SSxxS)
    
    
    
    
end


DATA.NAME_BASE = 'PhiDEF_inv_and_compl';
PlotModesDEF_SubdomainLevel(DATA,[PhiRB,PhiDEFbs],MESH);

disp('-----------------------------------------------------------------')
disp(['REACTION FORCES CAUSED BY basic deformational modes  REACT =  K*PhiDEFbs'])
disp(['-----------------------------------------------------------------------------------'])
DATA.NAME_BASE = '';
NameLoc =     'Reac_PhiDEFbs' ;
%     NAME_MODES_FOLDER = [cd,filesep,'MODES',filesep];
%     NAME_MODES_DISP = [NAME_MODES_FOLDER,DATA.NAME_BASE,NameLoc] ;
%     NameFileMesh = [NAME_MODES_DISP,'.msh'] ;
%     NameFile_res = [NAME_MODES_DISP,'.res'] ;
%
% %    disp(['Deformed shapes of the  deformational fict. interface modes: first ',num2str(size(VintDEF_INV,2)) ' most aligned with invariant modes ',[]])
%  %   disp(['Remaining  ',num2str(size(VintDEF_compl,2)) ' most orthogonal to invariants ',[]]) ;
%
%     DATALOCaaa = [];
%     GidPostProcessModesDOML(COORbnd,CNbREN,MESH.TypeElementB,[Kstiff(b,:)*PhiDEFbs],DATA.MESH.posgp,NameFileMesh,NameFile_res,[],DATALOCaaa) ;
%
%      NameLoc =     'FictDEFModAlignINVARIANTS+Complem' ;
%
PlotModesGID_Interf(COORbnd,CNbREN,MESH,[Kstiff(b,:)*PhiDEFbs],DATA,NameLoc) ;


%end
