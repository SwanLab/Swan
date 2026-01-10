function [BasisINTtruncated,BasisINT] = BasisInterfaceQ_methodsA(DATAINM,BasisUdef,SingVal_disp,DATA,BasisRdef,...
    Ui_Si,SingVal_reac,Ui_Si_reac,f1,f2,nMODESint,COORref,NODESfaces)


if DATAINM.USE_REACTIONS_FOR_INTERFACEMODES == 1
    BasisUdef = bsxfun(@times,BasisRdef{1}',SingVal_reac{1})';
    Ui_Si = Ui_Si_reac ;
elseif DATAINM.USE_REACTIONS_FOR_INTERFACEMODES == 0
    BasisUdef = bsxfun(@times,BasisUdef{1}',SingVal_disp{1})';
    
   % BasisUdef = BasisUdef(:,1:2) ; 
 
    
else
    error('Option not implemented')
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
    MULTIPLY_singular_VALUES_BASIS_MATRIX_INTERF =0; 
    for ipro = 1:length(NDISPpr)
        ifin = (iini-1) + NDISPpr(ipro) ;
        ifin = min(ifin,size(Ui_Si,2)) ; 
        BasisUdef_loc = Ui_Si(:,iini:ifin) ;
        Xa_i = [BasisUdef_loc(f1,:),BasisUdef_loc(f2,:)] ;
        [BasisINT_i S_i]= SVDT(Xa_i,0) ;
         if   MULTIPLY_singular_VALUES_BASIS_MATRIX_INTERF == 1
             S_i = S_i/(S_i(1)) ; 
             BasisINT_i = bsxfun(@times,BasisINT_i',S_i)';
         end

        %    BasisINT_i = bsxfun(@times,BasisINT_i',S_i)';
        nmodesi = min(nMODESint(ipro),length(S_i));
        BasisINT = [BasisINT  BasisINT_i(:,1:nmodesi) ];
        iini = ifin+1 ;
    end
    
    TOL_LOC =  DATAINM.TOL_SVD_InterfaceModes_PROJ ;
    
    [BasisINT S]= SVDT(BasisINT,TOL_LOC);
    Xa =BasisINT ;
end

if DATAINM.MODES_INTERFACE_RIGID_BODY == 0
    
    if isempty(BasisINT)
        [BasisINT S]= SVDT(Xa,0) ;
        
        %  nMODESint = max(1,nDEFr-DATAINM.SURPLUS_NUMBER_OF_INTERFACE_DISPLACEMENTS) ;
%         DATAINM = DefaultField(DATAINM,'NMODES_interface_CHOOSE',1:nMODESint) ; 
%         NMODESCHOOSE = DATAINM.NMODES_interface_CHOOSE ; 
        BasisINTtruncated = BasisINT(:,1:nMODESint) ;
    else
        BasisINTtruncated = BasisINT ;
    end
elseif DATAINM.MODES_INTERFACE_RIGID_BODY ~=0
    % Matrix BasisINTrb
    % -----------------------------------------
    DATAINM.ORTHOGONAL_RIGID_BODY_MODES= 2   ;
    
    BasisINTrb =  BasisRigidBodyInterface(COORref,NODESfaces,DATAINM)  ;
    
    %%%%% BasisINTdef
    % nDEFr = size(BasisRdef{1},2) ;
    
    DATAINM = DefaultField(DATAINM,'ModesINTERFACE_RigidB_select',1:size(BasisINTrb,2)) ;
    
    BasisINTrb = BasisINTrb(:,DATAINM.ModesINTERFACE_RigidB_select) ;
    nRB = size(BasisINTrb,2) ;
    if sum(nMODESint) < nRB
        error('When  MODES_INTERFACE_RIGID_BODY = 1,the number of reaction modes cannot be less than the number of rigid body modes')
    end
    Xa = Xa -BasisINTrb*(BasisINTrb\Xa) ; % Orthogonal complement
    TOL_LOC =  DATAINM.TOL_SVD_InterfaceModes_PROJ ;
    
    [BasisINTdef S]= SVDT(Xa,TOL_LOC) ;
    nTOT = min(length(S)+nRB,sum(nMODESint)) ;
    nMODESintDEF = nTOT-nRB ; %   max(0,nMODESintnDEFr-nRB-DATAINM.SURPLUS_NUMBER_OF_INTERFACE_DISPLACEMENTS) ;
    BasisINT =[BasisINTrb,BasisINTdef] ;
    BasisINTtruncated = [BasisINTrb, BasisINTdef(:,1:nMODESintDEF) ];
    
    
end




