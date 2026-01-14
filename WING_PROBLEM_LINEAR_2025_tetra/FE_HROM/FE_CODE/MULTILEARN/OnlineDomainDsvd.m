clc
clear all
% Structural analysis using domain-wise reduced-order data
% Prototype for one single domain
% See Domain_Decom_SVD.m
% Online stage. For the offline stage see DomainPartitionModal.m
%---------------------------------------------------------------
% INPUTS
%-----------------------------------------------------------
if exist('INPUT_PERIODIC')==0
    addpath('FE_CODE') ;
end
INPUT_DATAFILE  ='DATA_BEAMsolid' ; 'DATA_BEAMslice';  
DATA.NAMEWS = ['DATAWS/',INPUT_DATAFILE,'_OFFLINE','.mat'] ;
outputFE = ['DATAWS/',INPUT_DATAFILE,'_WS.mat'] ;

nDEF = 5; % Number of deformation modes
nREACT =nDEF ;
DATA.TypeUnitCell = 'HEXAG_2D_SQUARE';
% END INPUTS
% ----------------------------------------------------------


% Retrieval of reference mesh, coordinates, modes...
load(outputFE,'Bst','wSTs','Cglo','Fpnt','CNb','Tnod','TypeElementB','CONNECTb','K'); Kfe = K ; 
load(DATA.NAMEWS,'U','S','V','CNref','NODESref','COOR','TypeElement','posgp','BasisRdef','ElementsREF')

%,'Cglo','Bst','wSTs'

% Change renumbering
COORori = COOR ;
COOR = COOR(NODESref,:) ;
CN = CNref;
for inode = 1:length(NODESref)
    III = find(NODESref(inode) == CNref) ;
    CN(III) = inode ;
end
% % ----------------------------------
% Basis matrix of deformation modes
% --------------------------------
BasisUdef = U(:,1:nDEF) ;

% -----------------------------
% Basis matrix rigid body modes
% -----------------------------
% Reference point
[NODESfaces,NODEREF] =  PointPlanesRBODY(COOR,CN,DATA) ;
nnode = size(COOR,1) ;  ndim = size(COOR,2) ;
CoordinatesChange = repmat(COOR(NODEREF,:),nnode,1);  %
COORrel = COOR - CoordinatesChange ;  % Coordinates with respect to the reference point
BasisUrb =  zeros(ndim*nnode,3) ;
BasisUrb(1:ndim:end,1) = 1;
BasisUrb(2:ndim:end,2) = 1;
BasisUrb(1:ndim:end,3) = -COORrel(:,2);
BasisUrb(2:ndim:end,3) = COORrel(:,1);
%-----------------------------------------------
% Basis MAtrix for reaction forces
% ----------------------------------
% Body forces
f1 = NODESfaces{1} ;
f2 = NODESfaces{3} ;
f1 = small2large(f1,ndim) ;
f2 = small2large(f2,ndim) ;
f = [f1; f2] ;
BasisRrb = zeros(size(BasisUrb)) ;
BasisRrb(f,:) = BasisUrb(f,:) ;
% Deformation (self-equilibrium)
BasisRdef = BasisRdef(:,1:nREACT) ;

% Stiffness matrix reference domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of internal forces
ngaus = size(posgp,2) ;
ListGauss = small2large(ElementsREF,ngaus) ;
%Wdom =  wSTs(ListGauss);
if ndim == 2
    nstrain = 3 ;
end
%wDIAG = CompWeightDiag(Wdom,nstrain)  ;
ListGaussDOFS =  small2large(ListGauss,nstrain) ;
ListDOFs =  small2large(NODESref,ndim) ;
%%%
celasDOM = Cglo(ListGaussDOFS,ListGaussDOFS) ;
BstDOM = Bst(ListGaussDOFS,ListDOFs)*BasisUdef ;  % Celas already include the weights


Kdom = (BstDOM)'*(celasDOM*BstDOM) ;

%%%%%%%%%%%%%%%%5

FpntDOM = zeros(size(Fpnt)) ;
FpntDOM(ListDOFs) = Fpnt(ListDOFs)  ; % Point loads

CNbDOM = cell(size(CNb)) ;
TnodDOM = cell(size(CNbDOM)) ;
for idim = 1:ndim
    if ~isempty(CNb{idim})
        [CNbDOM{idim} ListEbndDOM]= ElemBnd(CNb{idim},NODESref) ;
        TnodDOM{idim} = Tnod{idim}(ListEbndDOM,:) ;
    end
end
ftracDOM = FtracCOMPvect(COORori,CNbDOM,TypeElementB,FpntDOM,TnodDOM,CONNECTb) ;
ftracDOM = ftracDOM(ListDOFs) ;

fextDOM = ftracDOM ;



%%%%%%%%%%%%%%%5

%%%%% SYSTEM MATRIX
nRB = 3;
K = {} ;
K{1,1} = BasisUdef(f1,:)'*BasisUdef(f1,:) ;
K{1,2} = Kdom ;
K{1,3} = zeros(size(Kdom,1),nREACT) ;
K{2,1} = Kdom ;
K{2,2}=  zeros(nDEF,nDEF) ;
K{2,3} = -BasisUdef'*BasisRdef ;
K{3,1} = K{1,3}' ;
K{3,2} = K{2,3}' ;
K{3,3} = BasisRdef(f2,:)'*BasisRdef(f2,:) ;

rRB = -(BasisUrb'*BasisRrb)\(BasisUrb'*fextDOM) ;
reactRB = BasisRrb*rRB ;

F{1,1} = zeros(size(K{1,1},1),1) ;
F{2,1} = BasisUdef'*fextDOM +  BasisUdef'*reactRB ;
F{3,1} = BasisRdef(f2,:)'*( -reactRB(f2)  ) ;

K = cell2mat(K) ;
F = cell2mat(F) ;

x = K\F ;

qDEF = x(1:nDEF) ;
LagrangeMult = x(nDEF+1:nDEF+nDEF) ;
rDEF = x(2*nDEF+1:end) ;

% Therefore
d = BasisUdef*qDEF;


%%
reactDEF = BasisRdef*rDEF ;
reactDOM = reactDEF+reactRB ;

% qDEFalt = Kdom\(BasisUdef'*(fextDOM +reactDOM))
% d = BasisUdef*qDEFalt;

% Reaction forces

% REsultants


stressDOM = celasDOM*BstDOM*qDEF ;

stressGLO = stressDOM;
ncomp = ngaus*nstrain*size(CN,1)  ;
if ndim == 3
    indGID = [1 2 3 6 4 5] ;
else
    indGID = [1 2 3 ] ;
end



for iglo = 1:nstrain
    stressGLO(iglo:nstrain:ncomp)  =  stressDOM(indGID(iglo):nstrain:ncomp) ;
end


 

GidPostProcess(COOR,CN,TypeElement,d,[], ...
    [],  reactDOM,[''],posgp,'ONLINE',[],DATA);






