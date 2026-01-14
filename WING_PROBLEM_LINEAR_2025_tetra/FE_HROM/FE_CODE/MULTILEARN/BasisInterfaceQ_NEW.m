function [BasisINTtruncated,BasisINT] = BasisInterfaceQ_NEW(DATAINM,BasisUdef,SingVal_disp,DATA,BasisRdef,...
    Ui_Si,SingVal_reac,Ui_Si_reac,f1,f2,nMODESint,COORref,NODESfaces)

if nargin == 0
    load('tmp.mat')
end
DATAINM = DefaultField(DATAINM,'MODES_INTERFACE_RIGID_BODY',0)   ;
DATAINM = DefaultField(DATAINM,'MODES_INTERFACE_INDVIDUAL_PROJECTS',0) ;
if DATAINM.MODES_INTERFACE_INDVIDUAL_PROJECTS ==0
    error(['Option non compatible with MODES_INTERFACE_RIGID_BODY = 2'])
end
if isempty(Ui_Si)
    error(['OPTION DATAINM.MODES_INTERFACE_INDVIDUAL_PROJECTS  only valid when modes are determined project-wise'])
end

APPROACH = 1;
NDISPpr= DATAINM.NMODES_PROJECT_DISPLACEMENTS ;

if  APPROACH == 1
WeightedSingVal =1 ;  % BasisINT is orthogonal, but not orthonomal (it is weighted by singular values)
% BASIS MATRIX FOR INTERFACE DISPLACEMENTS  (all projects)
 % BasisINTdisp  is a basis matrix for the displacement of the interfaces
 % of the element. It can be weighted by the corresponding singular
 % value
BasisINTdisp =  BasisIntf_projects(NDISPpr,Ui_Si,WeightedSingVal,nMODESint,f1,f2) ; 
% However, using BasisINTdisp as interface modes is conducive to
% instabilities. We are only interested in its projection onto the space
% spanned by reaction modes. Accordingly, we calculate now such a space (of the same dimension as BasisINTdisp) 
WeightedSingVal =0;
BasisINTreac =  BasisIntf_projects(NDISPpr,Ui_Si_reac,WeightedSingVal,nMODESint,f1,f2) ; 
% Projection 
Xa = BasisINTreac*(BasisINTreac'*BasisINTdisp) ; 

elseif APPROACH == 2
   error('Does not work')
% Compute basis matrix for reaction modes at the interface (common to f1 and f2)
   WeightedSingVal = 0 ; 
   Xr =  BasisIntf_projects(NDISPpr,Ui_Si_reac,WeightedSingVal,nMODESint,f1,f2) ;
   [BasisINTreac S]= SVDT(Xr,0) ;  
   
   iini = 1 ;
Xa = [] ;

for ipro = 1:length(NDISPpr)
    ifin = (iini-1) + NDISPpr(ipro) ;
    BasisUdef_loc = Ui_Si(:,iini:ifin) ;
    Xa_i = [BasisUdef_loc(f1,:),BasisUdef_loc(f2,:)] ;
    % Projection onto BasisINTreac
    Xa_proy = BasisINTreac*(BasisINTreac'*Xa_i) ; 
    Xa = [Xa Xa_proy] ; 
%     [BasisINT_i S_i]= SVDT(Xa_proy,0) ;
%     if   WeightedSingVal == 1
%         S_i = S_i/(S_i(1)) ;
%         BasisINT_i = bsxfun(@times,BasisINT_i',S_i)';
%     end
%     nmodesi = min(nMODESint(ipro),length(S_i));
%     BasisINT = [BasisINT  BasisINT_i(:,1:nmodesi) ];

    iini = ifin+1 ;
    
end

 

   % 
     Xa = BasisINTreac ; 
    
end
 
%     
%  BasisINT_i = bsxfun(@times,BasisINT_i',S_i)'
% TOL_LOC =  DATAINM.TOL_SVD_InterfaceModes_PROJ
% [BasisINT S]= SVDT(BasisINT,TOL_LOC);
% Xa =BasisINT ;

 TOL_LOC =  DATAINM.TOL_SVD_InterfaceModes_PROJ ;
if DATAINM.MODES_INTERFACE_RIGID_BODY == 0
    % No rigid body modes
   
    [BasisINT S]= SVDT(Xa,TOL_LOC) ; 
    BasisINTtruncated = BasisINT ; 
    
elseif DATAINM.MODES_INTERFACE_RIGID_BODY ~=0
    % Matrix rigid body modes
    % -----------------------------------------
    DATAINM.ORTHOGONAL_RIGID_BODY_MODES= 2   ;    
    BasisINTrb =  BasisRigidBodyInterface(COORref,NODESfaces,DATAINM)  ;
    DATAINM = DefaultField(DATAINM,'ModesINTERFACE_RigidB_select',1:size(BasisINTrb,2)) ;    
    BasisINTrb = BasisINTrb(:,DATAINM.ModesINTERFACE_RigidB_select) ;
    nRB = size(BasisINTrb,2) ;
    if sum(nMODESint) < nRB
        error('When  MODES_INTERFACE_RIGID_BODY = 1,the number of reaction modes cannot be less than the number of rigid body modes')
    end
    Xa = Xa -BasisINTrb*(BasisINTrb\Xa) ; % Orthogonal complement
   
    
    [BasisINTdef S]= SVDT(Xa,TOL_LOC) ;
    nTOT = min(length(S)+nRB,sum(nMODESint)) ;
    nMODESintDEF = nTOT-nRB ; %   max(0,nMODESintnDEFr-nRB-DATAINM.SURPLUS_NUMBER_OF_INTERFACE_DISPLACEMENTS) ;
    BasisINT =[BasisINTrb,BasisINTdef] ;
    BasisINTtruncated = [BasisINTrb, BasisINTdef(:,1:nMODESintDEF) ];
    
    
end




