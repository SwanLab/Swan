function [PhiDEFbs] = ComplementaryModes_basicF(PhiDEFbs_INV,DATAcommon,SNAPbasic,MdomCHOL,...
    MESH,PhiRB,Vintf,Mintf,DATA,MATPRO,OPERFE,Mdom,SNAPbasic_COMPLEM,Kstiff)
% Determination of complementary basic modes (complementary to invariant modes (PhiDEFbs_INV), which are determined in
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/DOMAINdecompHROM/GENERAL_BUBBLE/InvariantModes_basic.m)
% The actual output is PhiDEFbs, which are the set of basic, elastic modes.
% SNAPbasic --> All basic snapshots (elastic)
% Vintf ---> Fictitious interface modes
% PhiRB --> Rigid body modes, domain
% Modification of ComplementaryModes_basic.m (this is based on forces rather than on displacements)
% JAHO, 9-Apr-2024, CAmpus Nord UPC, Barcelona
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/02_BEAMQ8_IprofE.mlx
% -------------------------------
if nargin == 0
    load('tmp1.mat')
end


if ~isfield(DATAcommon.INFO_BASIC_DEFORMATIONAL_MODES,'NUMBER_OF_DEFORMATIONAL_MODES')
    error('You have not specified THE TOTAL NUMBER OF DEFORMATIONAL MODES FOR THIS TYPE OF ELEMENT')
    %  See for instance /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/LARGE_STRAINS/InputDataFunctions/Q8_for_beamBEHAVIOR_48.m
end

nDEF_basic = DATAcommon.INFO_BASIC_DEFORMATIONAL_MODES.NUMBER_OF_DEFORMATIONAL_MODES;
nINVmodes = size(PhiDEFbs_INV,2) ;

if nDEF_basic == nINVmodes
    disp(['No need for determining complementary deformational modes for elastic behavior'])
    PhiDEFbs = PhiDEFbs_INV ;
else
    disp('----------------------------------------------------------------------------------------------')
    disp(['Computing complementary basic modes (aside from those regarded as "invariant")    '])
    nDEF_compb = nDEF_basic-nINVmodes;
    disp(['Number of such complementary deformational modes =',num2str(nDEF_compb)])
    disp('--------------------------------------------------------------------------------')
    
    % STEP 1) WE BEGIN BY DETERMINING THE DEFORMATIONAL PART OF THE FICTITIOUS
    % INTERFACE MODES
    [VintfDEF,COORbnd,CNbREN,b] = DefPartInterfaceModes(MESH,PhiRB,Vintf,Mintf,DATA) ;
    
    % -------------------------------------------------------------------------------------------------------------
    %  STEP 2) DETERMINE THE nINV modes which are more aligned to the reactive modes caused by the invariant modes
    % --------------------------------------------------------------------------------------------------------------
    b = MESH.faceDOFSall;
    % Reactive forces at interface nodes
    PsiSE_INV = Kstiff(b,:)*PhiDEFbs_INV ;  %   
    
    MintfINV = inv(Mintf) ; 
    [ Basis_PsiSE_INV,SSp,VVp,MintfINV_chol]= WSVDT(PsiSE_INV,MintfINV) ; % (orthogonaliz.)
    % ------------------
    % Deformational fict. modes (orthogonaliz.)
    [VintDEF_orthogonalized,~,~,Mintf_chol ]= WSVDT(VintfDEF,Mintf) ; % These are the deformational fictitious interface modes (orhogonalized)
    % PRINCIPAL ANGLES AND PRINCIPAL VECTORS
    [UUpv,SSpv,VVpv] = SVDT(Basis_PsiSE_INV'*VintDEF_orthogonalized) ;
    
    VintDEF_INV = VintDEF_orthogonalized*VVpv ;
     disp(['Angles (in degrees) formed by INVARIANT fict. def. interface modes and  invariant reactive forces'])
    acosd(SSpv)
    % ------------------------------------------------
    % STEP 3) COMPLEMENTARY FICT. MODES     
    % ---------------------------------------------------
    VintDEF_compl = VintDEF_orthogonalized - VintDEF_INV*(VintDEF_INV'*Mintf*VintDEF_orthogonalized) ;
     
    DATAlocW.Mchol = Mintf_chol ; 
    DATAlocW.TOL = 1e-6; 
    [VintDEF_compl,SSpsv,VVpsv] = WSVDT(VintDEF_compl,[],DATAlocW) ;
    
    disp(['Angles (in degrees) formed by COMPLEMENTARY fict. def. interface modes and  invariant reactive forces'])
    
        [UUpv,SSpv,VVpv] = SVDT(Basis_PsiSE_INV'*VintDEF_compl) ;
        
        acosd(SSpv)

    
    % Plotting both sets of modes 
    
    NameLoc =     'FictDEFModAlignINVARIANTS+Complem' ;
    NAME_MODES_FOLDER = [cd,filesep,'MODES',filesep];
    NAME_MODES_DISP = [NAME_MODES_FOLDER,DATA.NAME_BASE,NameLoc] ;
    NameFileMesh = [NAME_MODES_DISP,'.msh'] ;
    NameFile_res = [NAME_MODES_DISP,'.res'] ;
    
    disp(['Deformed shapes of the  deformational fict. interface modes: first ',num2str(size(VintDEF_INV,2)) ' most aligned with invariant modes ',[]])
    disp(['Remaining  ',num2str(size(VintDEF_compl,2)) ' most orthogonal to invariants ',[]]) ; 
    
    DATALOCaaa = [];
    GidPostProcessModesDOML(COORbnd,CNbREN,MESH.TypeElementB,[VintDEF_INV,VintDEF_compl],DATA.MESH.posgp,NameFileMesh,NameFile_res,[],DATALOCaaa) ;
    %*******************************************************************************
    % STEP 4) COMPLEMENTARY MODES --> DEFORMATIONAL BASIS MATRIX
    %------------------------------------------------------------------------
    disp('Searching for the complementary modes...')
    disp('---------------------------------------------------')
    SNAPbasic_COMPLEM = SprojDEF_operator(PhiDEFbs_INV,Mdom,SNAPbasic_COMPLEM) ;
    DATAlocc.TOL = 1e-10;
    DATAlocc.Mchol = MdomCHOL ;
    [ PhiDEFcompl,SbsC,VbsC,~] = WSVDT( SNAPbasic_COMPLEM,[],DATAlocc) ;
    nCOMP_avail = size(PhiDEFcompl,2) ;
    disp(['Number of available deformational modes complementary to the invariants = ',num2str(nCOMP_avail)]) ;     
    if nDEF_compb > nCOMP_avail
        error('The number of available complementary modes is smaller than the number of necessary complementary modes ')
    end
    
    % ---------------------------------------------------------------------------
    % STEP 5) FORCES INDUCED BY THE COMPLEMENTARY MODES (Basis matrix)
    % ---------------------------------------------------------------------------
    DATAlocK.Mchol = MintfINV_chol ; 
    DATAlocK.TOL = 1e-10 ; 
    PsiDEFcompl = Kstiff(b,:)*PhiDEFcompl ; 
     [PsiDEFcompl, SSp,VVp]= WSVDT(PsiDEFcompl,[],DATAlocK) ; % (orthogonaliz.)
     
    % *********************************************
    % STEP 6) Subspace of   PsiDEFcompl most aligned with VintDEF_compl
    %*******************************************************************************
    
        [UUxxx,SSxx,VVxxx] = SVDT(PsiDEFcompl'*VintDEF_compl) ;
        
      disp(['Angles (in degrees) formed by COMPLEMENTARY fict. def. interface modes and  COMPLEMENTARY reactive forces'])
    acosd(SSxx)
    
    
     [UUxxxS,SSxxS,VVxxxS] = SVDT(PsiDEFcompl'*VintDEF_INV) ;
        
      disp(['Angles (in degrees) formed by INVARIANT fict. def. interface modes and  COMPLEMENTARY reactive forces'])
    acosd(SSxxS)
    
    error('Option not pursued any further....')
    
    
%     
%     
%     %indCOMPLbasic = setdiff(1:size(SNAPbasic,2),DATAcommon.INFO_BASIC_DEFORMATIONAL_MODES.INDICES_for_invariant_TRAINING_TESTS) ;
%     % Purging rigid-body displacements
%     SNAPbasicCOMPL = SprojDEF_operator(PhiRB,Mdom, SNAPbasic_COMPLEM) ;
%     DATAlocc.TOL = 1e-10;
%     DATAlocc.Mchol = MdomCHOL ;
%     [ PhiDEFcomple,Sbs,Vbs,~] = WSVDT( SNAPbasicCOMPL,[],DATAlocc) ;
%     
%     % COMPUTATION OF COMPLENTARY, BASIC DEFORMATIONAL MODES
%     % --------------------------------------------------------
%     
%     
%     % Recovering stiffness matrix
%     
%     
%     NORM_USED_SCALAR_PRODUCT ='Euclidean' ;  'Kbb'   ;  'Mintf' ;
%     
%     Kstiff =[] ;
%     
%     switch  NORM_USED_SCALAR_PRODUCT
%         case 'Kbb'
%             
%             
%             Kbb = Kstiff(b,b) ;
%             error('To be implemented....')
%             
%         case 'Mintf'
%             
%         case 'Euclidean'
%             PhiDEFbs = CompBasicDefMODES_euclid(VintfDEF,PhiDEFbs_INV,b,DATA,COORbnd,CNbREN,MESH,SNAPbasic_COMPLEM,Mdom,MdomCHOL)  ;
%             
%             
%         otherwise
%             
%             error('Option not implemented yet')
%             
%             
%     end
    
    
end

