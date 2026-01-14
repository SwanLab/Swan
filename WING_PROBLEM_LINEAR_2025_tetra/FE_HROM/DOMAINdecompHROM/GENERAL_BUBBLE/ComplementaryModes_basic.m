function [PhiDEFbs] = ComplementaryModes_basic(PhiDEFbs_INV,DATAcommon,SNAPbasic,MdomCHOL,...
    MESH,PhiRB,Vintf,Mintf,DATA,MATPRO,OPERFE,Mdom,SNAPbasic_COMPLEM,Kstiff)
% Determination of complementary basic modes (complementary to invariant modes (PhiDEFbs_INV), which are determined in
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/DOMAINdecompHROM/GENERAL_BUBBLE/InvariantModes_basic.m)
% The actual output is PhiDEFbs, which are the set of basic, elastic modes.
% SNAPbasic --> All basic snapshots (elastic)
% Vintf ---> Fictitious interface modes
% PhiRB --> Rigid body modes, domain
%
% JAHO, 7-Apr-2024, Molinos Marfagones, Cartagena
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
    disp(['Computing complementary basic tests (aside from those regarded as "invariant")    '])
    nDEF_compb = nDEF_basic-nINVmodes;
    disp(['Number of such complementary deformational modes =',num2str(nDEF_compb)])
    disp('--------------------------------------------------------------------------------')
    
    %indCOMPLbasic = setdiff(1:size(SNAPbasic,2),DATAcommon.INFO_BASIC_DEFORMATIONAL_MODES.INDICES_for_invariant_TRAINING_TESTS) ;
    % Purging rigid-body displacements
    SNAPbasicCOMPL = SprojDEF_operator(PhiRB,Mdom, SNAPbasic_COMPLEM) ;
    DATAlocc.TOL = 1e-10;
    DATAlocc.Mchol = MdomCHOL ;
    [ PhiDEFcomple,Sbs,Vbs,~] = WSVDT( SNAPbasicCOMPL,[],DATAlocc) ;
    
    % COMPUTATION OF COMPLENTARY, BASIC DEFORMATIONAL MODES
    % --------------------------------------------------------
    % Fictititious INTERFACE MODES, deformational part (VintfDEF)
    
    [VintfDEF,COORbnd,CNbREN,b] = DefPartInterfaceModes(MESH,PhiRB,Vintf,Mintf,DATA) ;
    
    
    % Recovering stiffness matrix
    
    
    NORM_USED_SCALAR_PRODUCT ='Euclidean' ;  'Kbb'   ;  'Mintf' ;
    
    Kstiff =[] ;
    
    switch  NORM_USED_SCALAR_PRODUCT
        case 'Kbb'
            
            
            Kbb = Kstiff(b,b) ;
            error('To be implemented....')
            
        case 'Mintf'
            
        case 'Euclidean'
            PhiDEFbs = CompBasicDefMODES_euclid(VintfDEF,PhiDEFbs_INV,b,DATA,COORbnd,CNbREN,MESH,SNAPbasic_COMPLEM,Mdom,MdomCHOL)  ;
            
            
        otherwise
            
            error('Option not implemented yet')
            
            
    end
    
    
end

