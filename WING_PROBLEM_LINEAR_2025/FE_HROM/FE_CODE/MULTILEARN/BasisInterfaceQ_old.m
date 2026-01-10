function [BasisINTtruncated,f1,f2] = ...
    BasisInterfaceQ_old(BasisUdef,SingVal_disp,DATAINM,COORref,CNref,DATA,BasisRdef,DATAOUT,Ui_Si,SingVal_reac,...
    Ui_Si_reac)

 
if nargin == 0
    load('tmp2.mat')
end


DATAINM = DefaultField(DATAINM,'NMODES_interface',[]) ;
if isempty(DATAINM.NMODES_interface)
    DATAINM = DefaultField(DATAINM,'SURPLUS_NUMBER_OF_INTERFACE_DISPLACEMENTS',0) ;
    nDEFr = size(BasisRdef{1},2) ;
    nMODESint =  nDEFr -DATAINM.SURPLUS_NUMBER_OF_INTERFACE_DISPLACEMENTS ;      %DATAINM.NMODES_interface ;
else
    nMODESint = DATAINM.NMODES_interface ;
end

ndim = size(COORref,2) ;
DATA = [] ;
[NODESbound,DATALIM] =  PointPlanesRBODY_GEN(COORref,CNref,DATA) ;
NODESfaces = NODESbound.PLANE ;
% Def. basis weighted by singular values
%DATAINM = DefaultField(DATAINM,'IntModesWeightedBySingValues',0); 
%if DATAINM.IntModesWeightedBySingValues == 1
DATAINM = DefaultField(DATAINM,'USE_REACTIONS_FOR_INTERFACEMODES',0); 
f1 = small2large(NODESfaces{1},ndim) ;
f2 = small2large(NODESfaces{3},ndim) ;
if DATAINM.USE_REACTIONS_FOR_INTERFACEMODES == 1
     BasisUdef = bsxfun(@times,BasisRdef{1}',SingVal_reac{1})';
     Ui_Si = Ui_Si_reac ; 
elseif DATAINM.USE_REACTIONS_FOR_INTERFACEMODES == 2
    BasisUdef = bsxfun(@times,BasisUdef{1}',SingVal_disp{1})';
    
else
    BasisUdef = bsxfun(@times,BasisUdef{1}',SingVal_disp{1})';
end
%else
 %    BasisUdef = BasisUdef{1} ; 
%end
%



DATAINM = DefaultField(DATAINM,'MODES_INTERFACE_RIGID_BODY',0)   ;


DATAINM = DefaultField(DATAINM,'MODES_INTERFACE_INDVIDUAL_PROJECTS',0) ; 
if DATAINM.MODES_INTERFACE_INDVIDUAL_PROJECTS ==0
   Xa = [BasisUdef(f1,:),BasisUdef(f2,:)] ;
   BasisINT = [] ;
else
    if isempty(Ui_Si)
        error(['OPTION DATAINM.MODES_INTERFACE_INDVIDUAL_PROJECTS  only valid when modes are determined project-wise'])
    end
    NDISPpr= DATAINM.NMODES_PROJECT_DISPLACEMENTS ;
    BasisINT = [] ;
    iini = 1 ;
    for ipro = 1:length(NDISPpr)
        ifin = (iini-1) + NDISPpr(ipro) ;
        BasisUdef_loc = Ui_Si(:,iini:ifin) ;
        Xa_i = [BasisUdef_loc(f1,:),BasisUdef_loc(f2,:)] ;
        [BasisINT_i S_i]= SVDT(Xa_i,0) ; 
    %    BasisINT_i = bsxfun(@times,BasisINT_i',S_i)';
        nmodesi = min(nMODESint(ipro),length(S_i)); 
        BasisINT = [BasisINT  BasisINT_i(:,1:nmodesi) ]; 
        iini = ifin+1 ;
    end
    [BasisINT S]= SVDT(BasisINT,0); 
    Xa =BasisINT ; 
end

if DATAINM.MODES_INTERFACE_RIGID_BODY == 0
    
    if isempty(BasisINT)
    [BasisINT S]= SVDT(Xa,0) ;
    
    %  nMODESint = max(1,nDEFr-DATAINM.SURPLUS_NUMBER_OF_INTERFACE_DISPLACEMENTS) ;
    BasisINTtruncated = BasisINT(:,1:nMODESint) ;
    else
        BasisINTtruncated = BasisINT ; 
    end
elseif DATAINM.MODES_INTERFACE_RIGID_BODY ~=0
    % Matrix BasisINTrb
    % -----------------------------------------
    COOR_FACE = COORref(NODESfaces{1},:) ;
    COORrefPOINT = sum(COOR_FACE,1)/size(COOR_FACE,1); % Center of gravity
    COORrel = bsxfun(@minus,COOR_FACE',COORrefPOINT')'; % Relative coordinates
    DATAINM.ORTHOGONAL_RIGID_BODY_MODES= 2   ;
    BasisINTrb = ConstructBasisRigidBody(COORrel,DATAINM) ; %
    %%%%% BasisINTdef
    % nDEFr = size(BasisRdef{1},2) ;
     
    DATAINM = DefaultField(DATAINM,'ModesINTERFACE_RigidB_select',1:size(BasisINTrb,2)) ; 
  
    BasisINTrb = BasisINTrb(:,DATAINM.ModesINTERFACE_RigidB_select) ; 
    nRB = size(BasisINTrb,2) ;
    if sum(nMODESint) < nRB
        error('When  MODES_INTERFACE_RIGID_BODY = 1,the number of reaction modes cannot be less than the number of rigid body modes')
    end
    Xa = Xa -BasisINTrb*(BasisINTrb\Xa) ; % Orthogonal complement
    [BasisINTdef S]= SVDT(Xa,0) ;
    nMODESintDEF = sum(nMODESint)-nRB ; %   max(0,nMODESintnDEFr-nRB-DATAINM.SURPLUS_NUMBER_OF_INTERFACE_DISPLACEMENTS) ;
    BasisINT =[BasisINTrb,BasisINTdef] ;
    BasisINTtruncated = [BasisINTrb, BasisINTdef(:,1:nMODESintDEF) ];
    
 
end







DATAINM = DefaultField(DATAINM,'PrintInterfaceModes',1) ;

if DATAINM.PrintInterfaceModes == 1
    
    %
    CNinterface = ElemBnd(DATAOUT.CONNECTb,NODESfaces{1});
    COORinterface = COORref ;
    
    NODESref =  DATAOUT.NODESref ;
    
    DOFl = [] ;
    DATA.NODES = NODESfaces{1}' ;
    
    
    
    posgp = [];
    NAME_MODES ='INTERFACE' ;
    disp('Printing interface modes')
    GidPostProcessModes(COORref,CNinterface,DATAOUT.TypeElementB,BasisINT,posgp,NAME_MODES,DATA,DOFl);
    disp('End')
end





% Matrix BasisINTrb
% -----------------------------------------
% f1 = NODESfaces{1} ; % Nodes involved face 1
% COOR_FACE = COORref(f1,:) ;
% COORrefPOINT = sum(COOR_FACE,1)/size(COOR_FACE,1); % Center of gravity
% COORrel = bsxfun(@minus,COOR_FACE',COORrefPOINT')'; % Relative coordinates
% BasisINTrb = ConstructBasisRigidBody(COORrel,DATAINM) ; %
%%%%% BasisINT
