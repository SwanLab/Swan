function [Z,Wred,Bred,BstD,nstrain,Hred,Cglo,wSTs,BasisS] =  ...
    ReducedPointsHomogenizationFLUC(nameWORKSPACE,DATAINfint,DATA_CUBATURE,DATA_GENGAUSS,DATAIN)

if nargin == 0
    load('tmp1.mat')
elseif nargin == 3
    DATA_GENGAUSS = [] ; DATAIN= [] ; 
elseif nargin == 4
    DATAIN = []; 
end

%%% LOAD DATA
% --------------------------------
load(nameWORKSPACE,'COOR','DOFr','DOFm','DISPLACEMENTS','Gb','strainINPglo','K','STRESSES','wSTs','Bst','Cglo',...
    'MaterialType','CN','Nst','CONNECTb','TypeElement','DISPmacro','posgp')

% Default Values
DATA_GENGAUSS = DefaultField(DATA_GENGAUSS,'TOL_LOC_InternalForces_SVD',1e-6) ;
DATAINfint = DefaultField(DATAINfint,'TOL_LOC_InternalForces',DATA_GENGAUSS.TOL_LOC_InternalForces_SVD) ;

DATA_GENGAUSS = DefaultField(DATA_GENGAUSS,'TOL_ECM',1e-6) ;
DATA_CUBATURE = DefaultField(DATA_CUBATURE,'TOL',DATA_GENGAUSS.TOL_ECM) ;



[nnode,ndim]= size(COOR) ;
nDOFS = nnode*ndim ;
% 1. Unconstrained DOFs
DOFf = (1:nnode*ndim)' ;
DOFf([DOFr ;DOFm]) = [] ;
DOFl = [DOFm ; DOFf] ;
% Matrix A, such that d= A d_L + U
A = sparse(nDOFS,length(DOFl));
A(DOFr,1:length(DOFm)) = Gb;
A(DOFm,1:length(DOFm)) = speye( length(DOFm));
A(DOFf,length(DOFm)+1:end) = speye( length(DOFf));


% 2. Matrix D  (relating DISPmacro with strainINPglo)
D = DISPmacro/strainINPglo ;

% 2. SVD displacements (DOFl), fluctutattions 
% ---------------------
FLUCT = DISPLACEMENTS-DISPmacro ; 


TOL_SVD =1e-6 ;
DATAsvd.RELATIVE_SVD  = 1;
[BasisU,Sl,Vl] = RSVDT(FLUCT(DOFl,:),TOL_SVD,0,0,DATAsvd) ;

BasisUplot = A*BasisU ; 
DATAprint  = [] ; DOFl_LOC  =[] ; 
DATAprint.NameFile_base = [DATAINfint.CURRENT_FOLDER,filesep,'GIDPOST',filesep,DATAINfint.NAME_DATA_INPUT_LOC,'_modesfluc'] ; 
NameFileMesh=[] ; 
DATAprint.MaterialType = MaterialType ; 
DATAprint = GidPostProcessModes_loc(COOR,CN,TypeElement,BasisUplot,posgp,NameFileMesh,DATAprint,DOFl_LOC);


% REduced "A" matrix
Ared = A*BasisU ;
% Reduced stiffness

CHECK_consistency =1;
if CHECK_consistency == 1
    
    Kred = Ared'*(K*Ared) ;
    
    % Input contribution
    fred = Ared'*K*D ;
    % Checking that everything is consistent
    q_train = -Kred\(fred*strainINPglo) ;
    d_reconst = Ared*q_train + D*strainINPglo;
    % HOMOGENIZATION
    % STRAINS
    error_DISPL = DISPLACEMENTS-d_reconst;
    error_DISPL = norm(error_DISPL,'fro')/norm(DISPLACEMENTS,'fro')*100  ; % CORRECT !!!
end


[BasisS,SingVal_stress,Vs] = RSVDT(STRESSES,TOL_SVD,0,0,DATAsvd) ;


DATAINfint = DefaultField(DATAINfint,'BRED_USING_THE_ENTIRE_BASISU',0) ; 
Bred = Bst*Ared ;
SingVal_stress = SingVal_stress/SingVal_stress(1) ; 
 

[BasisF, SingVal_F,nstrain,SNAPforceS,BasisS,VrightVal_F,DATAOUTfint]= BasisFfromStress(BasisS,SingVal_stress,...
    Bred, wSTs,DATAINfint) ;

% % Candidates points
 INDSEL = RestrictedDomainForECMpoints(DATA_GENGAUSS,COOR,Nst,CN,...
     CONNECTb) ;

DATA_CUBATURE.IND_POINTS_CANDIDATES = INDSEL ;
% Set of points and weights
[Z,Wred]= EmpiricalCubatureMethod(BasisF,SingVal_F,wSTs,DATA_CUBATURE) ;


% Printing reduced set of elements
%--------------------------------------
DATA_REFMESH.CN = CN ; DATA_REFMESH.COOR = COOR ;       DATA_REFMESH.MaterialType  = MaterialType ; 
  DATA_REFMESH.TypeElement  = TypeElement ; 
HYPERREDUCED_VARIABLES.WdomRED = Wred; 
DATAIN.NAME_WS_MODES = nameWORKSPACE ; 
ngaus = size(CN,2) ; 
HYPERREDUCED_VARIABLES.setElements = large2smallREP(Z,ngaus) ;

DATAIN = PrintingGid_ECMpoints(DATAIN,DATA_REFMESH,HYPERREDUCED_VARIABLES) ;
% -------------------------------------

 



BstD = Bst*D ;


Kred = Ared'*(K*Ared) ;

% Input contribution
fred = Ared'*K*D ;

Hred =  Kred\fred;






