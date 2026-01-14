function [BasisINTtruncated,f1,f2] = ...
    BasisInterfaceQ(BasisUdef,SingVal_disp,DATAINM,COORref,CNref,DATA,BasisRdef,DATAOUT,Ui_Si,SingVal_reac,Ui_Si_reac)


if nargin == 0
    load('tmp2.mat')
    DATAINM.USE_REACTIONS_FOR_INTERFACEMODES  = 2;
end

DATAINM = DefaultField(DATAINM,'TOL_SVD_InterfaceModes_PROJ',1e-5) ;
DATAINM = DefaultField(DATAINM,'NMODES_interface',[]) ;


if isempty(DATAINM.NMODES_interface)
    error('Specify DATAINM.NMODES_interface')
else
    nMODESint = DATAINM.NMODES_interface ;
end

ndim = size(COORref,2) ;
DATA = [] ;
[NODESbound,DATALIM] =  PointPlanesRBODY_GEN(COORref,CNref,DATA) ;
NODESfaces = NODESbound.PLANE ;
f1 = small2large(NODESfaces{1},ndim) ;
f2 = small2large(NODESfaces{3},ndim) ;
% Def. basis weighted by singular values
DATAINM = DefaultField(DATAINM,'USE_REACTIONS_FOR_INTERFACEMODES',0);

% Methods developed before 15-Dec-2017
if DATAINM.USE_REACTIONS_FOR_INTERFACEMODES <2
[BasisINTtruncated,BasisINT] = BasisInterfaceQ_methodsA(DATAINM,BasisUdef,SingVal_disp,DATA,BasisRdef,...
    Ui_Si,SingVal_reac,Ui_Si_reac,f1,f2,nMODESint,COORref,NODESfaces) ;
elseif DATAINM.USE_REACTIONS_FOR_INTERFACEMODES == 2
    error('This option has not been properly assessed')
    [BasisINTtruncated,BasisINT] = BasisInterfaceQ_NEW(DATAINM,BasisUdef,SingVal_disp,DATA,BasisRdef,...
    Ui_Si,SingVal_reac,Ui_Si_reac,f1,f2,nMODESint,COORref,NODESfaces) ;
else
    
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
