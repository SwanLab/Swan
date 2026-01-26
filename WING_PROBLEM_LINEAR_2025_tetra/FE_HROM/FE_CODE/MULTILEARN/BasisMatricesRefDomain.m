function [COORrel,BasisUrb,BasisUdef,BasisRdef,BasisRrb,f1,f2,fextDOM,KdomRED,f1NOD,f2NOD] ...
    = BasisMatricesRefDomain(COOR,NODESref,CNref,U,DATA,nDEF,nREACT,BasisRdef,posgp,ElementsREF,Cglo,...
    Bst,Fpnt,CNb,Tnod,TypeElementB,CONNECTb)

%%%%%%%%%%%%%%%%%%%% REFERENCE DOMAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%
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
[NODESfaces,NODEREF,f1NOD, f2NOD] =  PointPlanesRBODY(COOR,CN,DATA) ;
nnode = size(COOR,1) ;  ndim = size(COOR,2) ;
CoordinatesChange = repmat(COOR(NODEREF,:),nnode,1);  %
COORrel = COOR - CoordinatesChange ;  % Coordinates with respect to the reference point
if  ndim==2
    BasisUrb =  zeros(ndim*nnode,3) ;
    BasisUrb(1:ndim:end,1) = 1;
    BasisUrb(2:ndim:end,2) = 1;
    BasisUrb(1:ndim:end,3) = -COORrel(:,2);
    BasisUrb(2:ndim:end,3) = COORrel(:,1);
else
    
    BasisUrb =  zeros(ndim*nnode,6) ;
    BasisUrb(1:3:end,1) = 1 ;
    BasisUrb(2:3:end,2) = 1;
    BasisUrb(3:3:end,3) = 1;
    % Rotation modes
    BasisUrb(2:3:end,4) = COORrel(:,3);
    BasisUrb(3:3:end,4) = -COORrel(:,2);
    
    BasisUrb(1:3:end,5) = -COORrel(:,3);
    BasisUrb(3:3:end,5) = COORrel(:,1);
    
    BasisUrb(1:3:end,6) = COORrel(:,2);
    BasisUrb(2:3:end,6) = -COORrel(:,1);
    
    
    %      MODES_ROT = zeros(ndim*nnode,3) ;
    %     MODES_ROT(2:3:end,1) =  COORrve{idomain}(:,3) ;
    %     MODES_ROT(3:3:end,1) =  -COORrve{idomain}(:,2) ;
    %
    %     MODES_ROT(1:3:end,2) =  -COORrve{idomain}(:,3) ;
    %     MODES_ROT(3:3:end,2) =  COORrve{idomain}(:,1) ;
    %
    %     MODES_ROT(1:3:end,3) =  COORrve{idomain}(:,2) ;
    %     MODES_ROT(2:3:end,3) =  -COORrve{idomain}(:,1) ;
    
    
end
%-----------------------------------------------
% Basis MAtrix for reaction forces
% ----------------------------------
f1 = small2large(f1NOD,ndim) ; % Degrees of freedom face F1 (local)
f2 = small2large(f2NOD,ndim) ; % Degrees of freedom face F2 (local)
f = [f1; f2] ;
BasisRrb = zeros(size(BasisUrb)) ;
BasisRrb(f,:) = BasisUrb(f,:) ;
% Deformation (self-equilibrium)
BasisRdef = BasisRdef(:,1:nREACT) ;

%---------------------------------------------
% Stiffness matrix reference domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of internal forces
ngaus = size(posgp,2) ;
ListGauss = small2large(ElementsREF,ngaus) ;
%Wdom =  wSTs(ListGauss);
if ndim == 2
    nstrain = 3 ;
elseif ndim ==3
    nstrain = 6;
    
end
%wDIAG = CompWeightDiag(Wdom,nstrain)  ;
%dbstop('87')
ListGaussDOFS =  small2large(ListGauss,nstrain) ;
ListDOFs =  small2large(NODESref,ndim) ;
%%%
celasDOM = Cglo(ListGaussDOFS,ListGaussDOFS) ;
BstDOM = Bst(ListGaussDOFS,ListDOFs)*BasisUdef ;  % Celas already include the weights
KdomRED = (BstDOM)'*(celasDOM*BstDOM) ;
%%%%%%%%%%%%%%%%5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRACTION FORCES = EXETERNAL FORCES (for a single domain )
FpntDOM = zeros(size(Fpnt)) ;
FpntDOM(ListDOFs) = Fpnt(ListDOFs)  ; % Point loads

CNbDOM = cell(size(CNb)) ;
TnodDOM = cell(size(CNbDOM)) ;
for idim = 1:ndim
    if ~isempty(CNb{idim})
        %  dbstop('106')
        [CNbDOM{idim} ListEbndDOM]= ElemBnd_noremove(CNb{idim},NODESref) ;
        TnodDOM{idim} = Tnod{idim}(ListEbndDOM,:) ;
    end
end
ftracDOM = FtracCOMPvect(COORori,CNbDOM,TypeElementB,FpntDOM,TnodDOM,CONNECTb) ;
ftracDOM = ftracDOM(ListDOFs) ;
fextDOM = ftracDOM ;