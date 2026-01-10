function [Z,Wred,Bred,BstD,nstrain,Hred,Cglo,wSTs,BasisS] =  ...
    ReducedPointsHomogenization(nameWORKSPACE,DATAINfint,DATA_CUBATURE,DATA_GENGAUSS,DATAIN)

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
    'MaterialType','CN','Nst','CONNECTb','TypeElement')

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


% 2. Matrix D
Umatrix = DISPLACEMENTS -A*DISPLACEMENTS(DOFl,:) ;
D = Umatrix/strainINPglo ;

% 2. SVD displacements (DOFl)
% ---------------------
TOL_SVD =0 ;
DATAsvd.RELATIVE_SVD  = 1;
[BasisU,Sl,Vl] = RSVDT(DISPLACEMENTS(DOFl,:),TOL_SVD,0,0,DATAsvd) ;

% REduced "A" matrix
Ared = A*BasisU ;
% Reduced stiffness

CHECK_consistency =0;
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

if DATAINfint.BRED_USING_THE_ENTIRE_BASISU ==1
    error('This option was abandoned in favor of  .... == 0 ')
%     TOL_SVD =0 ;
%     DATAsvd.RELATIVE_SVD  = 1;
%     [BasisUall,Sl,Vl] = RSVDT(DISPLACEMENTS,TOL_SVD,0,0,DATAsvd) ;
    BredINP = [Bred,Bst*D];  
else
    BredINP = Bred ;
end

[BasisF, SingVal_F,nstrain,SNAPforceS,BasisS,VrightVal_F,DATAOUTfint]= BasisFfromStress(BasisS,SingVal_stress,...
    BredINP, wSTs,DATAINfint) ;

% % Candidates points
% INDSEL = RestrictedDomainForECMpoints(DATA_GENGAUSS,COOR,Nst,CN,...
%     CONNECTb) ;

%DATA_CUBATURE.IND_POINTS_CANDIDATES = INDSEL ;
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


% 
% 
%  Kred = Ared'*(K*Ared) ;
%     
%     % Input contribution
%     fred = Ared'*K*D ;
%     % Checking that everything is consistent
%     q_train = -Kred\(fred*strainINPglo) ;
  %   d_reconst = Ared*q_train + D*strainINPglo;D*strainINPglo;
%     % HOMOGENIZATION
%  

BD = Bst*D ; 

%%%%%%%
Zindexes = small2large(Z,nstrain);
CgloRED =  Cglo(Zindexes,Zindexes) ; 
BredAST =  Bred(Zindexes,:) ; 
BDast=  BD(Zindexes,:) ; 
WFE_red= wSTs(Z)  ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q_train = -Kred\(fred*strainINPglo) ;

STRAINS_reducedpoints =  BredAST*q_train +  BDast*strainINPglo; ; 

STRESSES_reducedpoints = CgloRED*STRAINS_reducedpoints  ; 


for istrain = 1:nstrain 
    COMP = istrain:nstrain:size(STRESSES_reducedpoints,1) ; 
    for  itraining = 1:size(STRESSES_reducedpoints,2)
    STRESSES_reducedpoints(COMP,itraining) = STRESSES_reducedpoints(COMP,itraining)./WFE_red ;
    end
end


  %ITRAINING = 4 ; 

  


STRESSES_FE = STRESSES(Zindexes,:) ;  


errorREDU = norm(STRESSES_reducedpoints-STRESSES_FE)/norm(STRESSES_FE) ; 




     coeff = (BasisS(Zindexes,:)'*BasisS(Zindexes,:))\BasisS(Zindexes,:)' ;
   ReconsStresses = BasisS*coeff ;

   
   
   STRESSES_RED_ALL = ReconsStresses*STRESSES_reducedpoints  ; 
   
   
   errorALL = norm(STRESSES_RED_ALL-STRESSES)/norm(STRESSES) ; 



