function [PhiDEFbs_INV,MdomCHOL,SNAPbasic_COMPLEM] = InvariantModes_basic(SNAPbasic,DATAcommon,PhiRB,Mdom,DATAoffline,MESH,DATA,Kstiff)
% Determination of invariant modes, for elastic deformational modes
% JAHO, 7-Apr-2024, Molinos Marfagones, Cartagena 
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/02_BEAMQ8_IprofE.mlx
% -------------------------------
if nargin == 0
    load('tmp1.mat')
end
b = MESH.faceDOFSall;

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

if ~isfield(DATAcommon.INFO_BASIC_DEFORMATIONAL_MODES,'INDICES_all_prescribed_TRAINING_TESTS')

        error('You have not specified which are the indices of the tests in which we have to search for the complementary modes to those which are invariantes')

end


% Purging rigid-body displacements
SNAPinv = SprojDEF_operator(PhiRB,Mdom,SNAPbasic(:,DATAcommon.INFO_BASIC_DEFORMATIONAL_MODES.INDICES_for_invariant_TRAINING_TESTS)) ;

SNAPbasic_COMPLEM = SprojDEF_operator(PhiRB,Mdom,SNAPbasic(:,DATAcommon.INFO_BASIC_DEFORMATIONAL_MODES.INDICES_all_prescribed_TRAINING_TESTS)) ;
 

%DATAoffline.TOLSVD_INVARIANT_BASIC_MODES = 1e-4 ;   % Determine invariant modes
DATAoffline = DefaultField(DATAoffline,'TOLSVD_INVARIANT_BASIC_MODES',1e-4) ;


DATAlocc.TOL = DATAoffline.TOLSVD_INVARIANT_BASIC_MODES;
[ PhiDEFbs,Sbs,Vbs,MdomCHOL] = WSVDT( SNAPinv,Mdom,DATAlocc) ;

disp(['Singular Values Deformational Modes (in the search for invariant modes)'])

Sbs/Sbs(1)

RATIO_SV = Sbs(1:end-1)./Sbs(2:end) ;

if DATAoffline.TOLSVD_INVARIANT_BASIC_MODES >0 
    THinv = 1e20 ; 
else
    THinv = 100 ; 
end

DATAoffline = DefaultField(DATAoffline,'THRESHOLD_RATIO_CONSECUTIVE_SINGLE_VALUES_INVARIANT_DEF_MODES',THinv) ;


III =   find(RATIO_SV>DATAoffline.THRESHOLD_RATIO_CONSECUTIVE_SINGLE_VALUES_INVARIANT_DEF_MODES) ;

if isempty(III)
    nINVmodes = length(Sbs) ; 
else
nINVmodes = III(1) ;
end
disp('***********************************************************************************************************')
disp('Searching for invariant modes')
disp(['Total number of  deformational modes (for the first set of training tests) = ',num2str(length(Sbs))]) ;
disp(['Number of  deformational modes that may be deemed INVARIANT = ',num2str(nINVmodes)]) ;
if ~isempty(III)
disp(['Jump in singular values (ratio consecutive singular values) =',num2str(RATIO_SV(nINVmodes))])
disp('***********************************************************************************************************')
else
    disp(['tolerance SVD = ',num2str(DATAoffline.TOLSVD_INVARIANT_BASIC_MODES)])
end

% USE_CRITERIO_FORCES  =0 ; 
% 
% if  USE_CRITERIO_FORCES == 1
% disp(['Criterion based on forces'])
% [PhiDEFbs_b,SSdd,VVdd ]= SVDT(PhiDEFbs(b,1:nINVmodes)) ; 
% Forces = Kstiff*PhiDEFbs  ; 
% [PsiDEFb,SSbf,VVbf ]= SVDT(Forces(b,1:nINVmodes))  ; 
% 
% [UUc,SSc,VVc] = SVDT(PsiDEFb'*PhiDEFbs_b) ; 
% 
% 
% end



PhiDEFbs_INV = PhiDEFbs(:,1:nINVmodes) ;

DATA.NAME_BASE = 'Invariants';
PlotModesDEF_SubdomainLevel(DATA,PhiDEFbs_INV,MESH);

disp(['*****************************************************************************************'])



[UUUdef,SSS,VVV] = SVDT(PhiDEFbs_INV(b,:)) ; 
disp(['Singular values   PhiDEFbs_INV(b,:)'])
SSS/SSS(1)

PsiSE_INV = Kstiff*PhiDEFbs_INV ; 

[PsiSE_INV,SSSse,VVVse] = SVDT(PsiSE_INV(b,:)) ; 

disp(['Singular values  self-equilibrated modes arising from PsiSE_INV'])
SSSse/SSSse(1)

disp(['Cosine angles formed by invariant deformational modes and their conjugate interface forces'])

[Uconj,Sconj,Vconj ] = SVDT(PsiSE_INV'*UUUdef) ; 
Sconj

