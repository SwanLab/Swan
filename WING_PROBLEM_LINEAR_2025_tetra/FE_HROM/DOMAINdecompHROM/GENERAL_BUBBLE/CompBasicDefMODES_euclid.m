function PhiDEFbs = CompBasicDefMODES_euclid(VintfDEF,PhiDEFbs_INV,b,DATA,COORbnd,CNbREN,MESH,SNAPbasic_COMPLEM,Mdom,MdomCHOL)
% This function finds the basic, complementary modes to those regarded as
% "invariant"
% VintDEF = Deformational component of the interface, polynomial modes
% PhiDEFbs_INV = Invariant modes
% b = Boundary DOFS
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/02_BEAMQ8_IprofE.mlx
% JAHO, 7-Apr-2024, Sunday, Molinos Marfagones, Cartagena
% -------------------------------------------------------------------------------------------------------------------
if nargin == 0
    load('tmp.mat')
end


% Most aligned of PhiDEFbs_INV(b,:) with VintfDEF



%% S1)
VintDEF_orthogonalized = SVDT(VintfDEF) ; % These are the deformational fictitious interface modes (orhogonalized)

%%% S2)
[PhiDEFbs_INV_b_orthogonalized,SSbinv,VVbinv] = SVDT(PhiDEFbs_INV(b,:)) ;   % These are the invariant modes at the interface, orthogonalized
% Notice that
% PhiDEFbs_INV_b_orthogonalized*diag(SSbinv)*VVbinb' = PhiDEFbs_INV(b,:)
% Thus PhiDEFbs_INV_b_orthogonalized =  PhiDEFbs_INV(b,:)*VVbinb*diag(SSbinv)^{-1}
% The above is for the interface nodes. For all the nodes, we have
% that
%  PhiDEFbs_INV_b_orthogonalized_ALL = PhiDEFbs_INV(b,:)*VVbinb*diag(SSbinv)^{-1}
% S3)
coeff = bsxfun(@times,VVbinv',1./SSbinv)'  ;
PhiDEFbs_INV_b_orthogonalized_ALL = PhiDEFbs_INV*coeff ;  % This is the extension of PhiDEFbs_INV_b_orthogonalized for all DOFs


%%%% S4...Principal angles
[UU,SS,VV] = SVDT(VintDEF_orthogonalized'*PhiDEFbs_INV_b_orthogonalized) ;
Vdef_aligned_inv = VintDEF_orthogonalized*UU ;  % These are the nINV fictitious modes most aligned with invariant modes
PhiDEFbs_INV_aligned_interface = PhiDEFbs_INV_b_orthogonalized_ALL*VV ;   % These are the  invariant modes, reorderd so as
% to match  Vdef_aligned_inv

disp(['Cosines angled formed by invariant def. modes and interface modes'])
SS

NameLoc =     'DispIntf_DEF_invariant' ;
NAME_MODES_FOLDER = [cd,filesep,'MODES',filesep];
NAME_MODES_DISP = [NAME_MODES_FOLDER,DATA.NAME_BASE,NameLoc] ;
NameFileMesh = [NAME_MODES_DISP,'.msh'] ;
NameFile_res = [NAME_MODES_DISP,'.res'] ;

disp(['Deformed shapes of deformational interface modes aligned with invariant modes (GID) --> '])
DATALOCaaa = [];

GidPostProcessModesDOML(COORbnd,CNbREN,MESH.TypeElementB,[Vdef_aligned_inv],DATA.MESH.posgp,NameFileMesh,NameFile_res,[],DATALOCaaa) ;


disp(['Deformational invariant modes corresponding to the preceding interface modes (GID) -->'])





DATA.NAME_BASE = 'Invariants_aligned_intf';
PlotModesDEF_SubdomainLevel(DATA,PhiDEFbs_INV_aligned_interface,MESH);


disp('----------------------------------------------------------------------------')

% BASIS VECTORS WHICH ARE ORTHOGONAL TO Vdef_aligned_inv
Vdef_orthog_inv=  VintfDEF - Vdef_aligned_inv*Vdef_aligned_inv'*VintfDEF ;
[Vdef_orthog_inv,SS,VV] = SVDT(Vdef_orthog_inv) ;

disp(['Deformed shapes of deformational interface modes aligned with orthogonal complement of  invariant modes (GID) --> '])
DATALOCaaa = [];

NameLoc =     'DispIntf_DEF_orth_compl_invariant' ;
NAME_MODES_FOLDER = [cd,filesep,'MODES',filesep];
NAME_MODES_DISP = [NAME_MODES_FOLDER,DATA.NAME_BASE,NameLoc] ;
NameFileMesh = [NAME_MODES_DISP,'.msh'] ;
NameFile_res = [NAME_MODES_DISP,'.res'] ;

GidPostProcessModesDOML(COORbnd,CNbREN,MESH.TypeElementB,[Vdef_orthog_inv],DATA.MESH.posgp,NameFileMesh,NameFile_res,[],DATALOCaaa) ;

% LET US COPE  NOW WITH THE "REAL" PART OF THE INTERFACE DISPLACEMENTS
% SNAPbasic_COMPLEM, PhiDEFbs_INV

% ORTHGONAL COMPLEMENT OF PhiDEFbs_INV 

SNAPbasic_COMPLEM = SprojDEF_operator(PhiDEFbs_INV,Mdom,SNAPbasic_COMPLEM) ;


%%% S2)
 TOL = 1e-6; 
DATAlocA.ISREL = 1; 
[PhiCOMPLEb,SSbinv,VVbinv] = SVDT(SNAPbasic_COMPLEM(b,:),TOL,DATAlocA ) ;   % These are the invariant modes at the interface, orthogonalized
coeff = bsxfun(@times,VVbinv',1./SSbinv)'  ;
PhiCOMPLE_all = SNAPbasic_COMPLEM*coeff ;  % This is the extension of PhiDEFbs_INV_b_orthogonalized for all DOFs

[UU,SS,VV] = SVDT(Vdef_orthog_inv'*PhiCOMPLEb) ;
Vdef_orthog_inv_aligned = Vdef_orthog_inv*UU ;  % These are the nINV fictitious modes most aligned with invariant modes
PhiDEFbs_COMP = PhiCOMPLE_all*VV ;   % These are the  invariant modes, reorderd so as

 disp(['Cosines angled formed by real complementary interface modes and  fictititious interface modes (compl. to invariants)'])
SS


disp(['Deformed shapes of complementary deformational   modes  --> '])
DATALOCaaa = [];




DATA.NAME_BASE = 'DEFMODES_complentary_basic';
PlotModesDEF_SubdomainLevel(DATA,PhiDEFbs_COMP,MESH);

DATALOCaaaa.Mchol = MdomCHOL; 
PhiDEFbs_COMP = WSVDT(PhiDEFbs_COMP,[],DATALOCaaaa) ; 

PhiDEFbs = [PhiDEFbs_INV,PhiDEFbs_COMP] ; 




disp(['Deformed shapes of fictititious complementary interface     modes (aligned with their domain counterpart)  --> '])
DATALOCaaa = [];

NameLoc =     'DispIntf_DEF_orth_compl_invariant' ;
NAME_MODES_FOLDER = [cd,filesep,'MODES',filesep];
NAME_MODES_DISP = [NAME_MODES_FOLDER,DATA.NAME_BASE,NameLoc] ;
NameFileMesh = [NAME_MODES_DISP,'.msh'] ;
NameFile_res = [NAME_MODES_DISP,'.res'] ;

GidPostProcessModesDOML(COORbnd,CNbREN,MESH.TypeElementB,[Vdef_orthog_inv_aligned],DATA.MESH.posgp,NameFileMesh,NameFile_res,[],DATALOCaaa) ;

 [UUU,SSS,VVV] = SVDT(PhiDEFbs(b,:)) ; 
disp(['Singular values   PhiDEFbs(b,:)'])
SSS/SSS(1)
