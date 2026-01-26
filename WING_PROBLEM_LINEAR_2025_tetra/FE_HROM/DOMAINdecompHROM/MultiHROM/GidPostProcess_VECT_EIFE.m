function strainCOARSE_history = GidPostProcess_VECT_EIFE(MESH,DATA,icluster,DATAinpGID,PROPMAT,DATALOC)

if nargin==0
    load('tmp.mat')
    
end
COOR = MESH.COOR;
CN = MESH.CN ;
TypeElement = MESH.TypeElement;
MaterialType = MESH.MaterialType;

DATALOC.DATAcommon = DefaultField(DATALOC.DATAcommon,'FACTOR_MULTIPLYING_DEFORMATIONAL_RB_amplitudes',[1,1]) ; % Just For visualization purposes
DATA.FACTOR_MULTIPLYING_DEFORMATIONAL_RB_amplitudes = DATALOC.DATAcommon.FACTOR_MULTIPLYING_DEFORMATIONAL_RB_amplitudes ;
NameFile_msh = DATA.PRINT.NAME_FILE_MSH{icluster} ;
NameFile_res= DATA.PRINT.NAME_FILE_RES{icluster} ;

% % For domain decomposition purpooses (not useful now )
% DATA = DefaultField(DATA,'MATERIAL_ORIGINAL',[]) ;
% DATA.MATERIAL =  DATA.MATERIAL_ORIGINAL ;
% DATA = DefaultField(DATA,'PostProcessWithNOSLICES',1) ;
% if DATA.PostProcessWithNOSLICES == 1 && ~isempty(DATA.MATERIAL) && (~isempty(DATA.MakeMeshByRepetition.nDOM) && DATA.MakeMeshByRepetition.nDOM(1)>1 )
%     nmat = length(DATA.MATERIAL.PLY) ;
%     ndom = prod(DATA.MakeMeshByRepetition.nDOM); %(1) ;
%     TOTnmat = nmat*ndom ;
%     MAT_TypeORIG = (1:TOTnmat) ;
%     MAT_TypeNOSLICES = repmat((1:nmat)',ndom,1) ;
%     NewMaterialType = zeros(size(MaterialType));
%
%     for imat = 1:length(MAT_TypeORIG)
%         INDLOC =find(MaterialType==imat) ;
%         NewMaterialType(INDLOC) = MAT_TypeNOSLICES(imat) ;
%     end
%
%     MaterialType =  NewMaterialType ;
% end

%
% if ~iscell(CN)
%     DATA.NAMEMESH = {'MESH'} ;
%     CN = {CN} ;
%     MaterialType = {MaterialType} ;
%     TypeElement = {TypeElement} ;
%     GAUSSV_PROP = {DATA.PRINT.GAUSSV_PROP} ;
%     %  GAUSS_SNAP = DATA.STORE.NAME_MATFILE_STORE{icluster};
%     posgp = {DATA.MESH.posgp} ;
% else
%     error('Option not implemented...9-Dec-2020')
%     GAUSSV_PROP = DATA.PRINT.GAUSSV_PROP ;
%     %  GAUSS_SNAP = DATA.PRINT.GAUSS_SNAP ;
%     posgp  = DATA.MESH.posgp;
% end
if exist(DATA.STORE.NAME_MATFILE_STORE{icluster})
    if DATA.STORE.COMPRESS_WITH_SVD == 1
        load( DATA.STORE.NAME_MATFILE_STORE{icluster},'SNAP_cluster','STEP_LOC') ;
    else
        load( DATA.STORE.NAME_MATFILE_STORE{icluster},'SNAP','STEP_LOC') ;
    end
else
    return
end


[STEPS_PRINT,iaaa,ibbb ]= intersect(STEP_LOC,DATA.PRINT.NSTEPS) ;

if DATA.STORE.COMPRESS_WITH_SVD == 1
    fff= fieldnames(SNAP_cluster) ;
    SNAP = [];
    for  iii = 1:length(fff)
        if ~isempty( SNAP_cluster.(fff{iii}).V)
            V = SNAP_cluster.(fff{iii}).V(iaaa,:) ;
            S = SNAP_cluster.(fff{iii}).S  ;
            SV = bsxfun(@times,V',S)';
            SNAP.(fff{iii}) = SNAP_cluster.(fff{iii}).U*SV' ;
            
        else
            SNAP.(fff{iii}) = [];
        end
        
    end
else
    fff= fieldnames(SNAP) ;
    for  iii = 1:length(fff)
        SNAP.(fff{iii}) =   SNAP.(fff{iii})(:,iaaa) ;
    end
    
end
TIME_PRINT = DATA.STEPS(STEPS_PRINT) ;
DATA.TIME_PRINT = TIME_PRINT ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = SNAP.DISP;%(:,STEPS_PRINT) ; % JAHO, 19-Nov-2023


if size(CN,2) == 26
    % Hexahedra element in which we have eliminated the interior point
    % We need now to restore such a point in order to be processesd by GID
    COORintPOINT = zeros(size(CN,1),size(COOR,2)) ;
    nnodeE = 26 ;
    ndim = 3;
    for  inodeE = 1:nnodeE
        COORintPOINT  = COORintPOINT + COOR(CN(:,inodeE),:) ;
    end
    COORintPOINT = COORintPOINT/nnodeE ;
    % DISPLACEMENTS
    dMIDPOINT = zeros(size(COORintPOINT,1)*ndim,size(d,2)) ;
    for itime = 1:size(d,2)
        dLOC = reshape(d(:,itime),ndim,[] )' ;
        dMIDPOINT_loc = zeros(size(CN,1),ndim) ;
        for  inodeE = 1:nnodeE
            dMIDPOINT_loc   = dMIDPOINT_loc + dLOC(CN(:,inodeE),:) ;
        end
        dMIDPOINT_loc = dMIDPOINT_loc'/nnodeE;
        dMIDPOINT(:,itime) =  dMIDPOINT_loc(:) ;  ;
    end
    
    nnodeOLD  = size(COOR,1) ;
    OTHERinputs.COORprintCOARSE = [COOR;COORintPOINT] ;
    CNintpoint = nnodeOLD + (1:size(CN,1))';
    OTHERinputs.CNprintCOARSE = [CN,CNintpoint] ;
    OTHERinputs.dPRINTcoarse = [d; dMIDPOINT] ;
    
else
    
    OTHERinputs.COORprintCOARSE = [] ; OTHERinputs.CNprintCOARSE = [] ; OTHERinputs.dPRINTcoarse = [] ;
end

stressGLO_REF = SNAP.PK2STRESS; %(:,STEPS_PRINT);  % These are the stresses in the unrotated configuration
DATA.NameFile_msh = NameFile_msh;
DATA.NameFile_res = NameFile_res;
if ~isempty(SNAP.VONMISES_CAUCHY_STRESS)
    VON_mises_COARSE= zeros(size(CN,1),size(SNAP.VONMISES_CAUCHY_STRESS,2));
    for itime = 1:size(SNAP.VONMISES_CAUCHY_STRESS,2)
        VONloc  = reshape(SNAP.VONMISES_CAUCHY_STRESS(:,itime),DATA.MESH.ngaus_STRESS,[]) ;
        VON_mises_COARSE(:,itime) = max(VONloc,[],1)  ;
    end
  %  VON_mises_COARSE = VON_mises_COARSE; (:,STEPS_PRINT);
else
    VON_mises_COARSE = [] ;
end

SNAP = DefaultField(SNAP,'InternalVarStrain',[]) ;
if ~isempty(SNAP.InternalVarStrain)
    InternalVarStrain= zeros(size(CN,1),size(SNAP.InternalVarStrain,2));
    for itime = 1:size(SNAP.InternalVarStrain,2)
        VONloc  = reshape(SNAP.InternalVarStrain(:,itime),DATA.MESH.ngaus_STRESS,[]) ;
        InternalVarStrain(:,itime) = max(VONloc,[],1)  ;
    end
  %  InternalVarStrain = InternalVarStrain(:,STEPS_PRINT) ;
else
    InternalVarStrain = [] ;
end





DATALOC = DefaultField(DATALOC,'FineScaleDomainsToShow',[])  ;
DATA.FineScaleDomainsToShow = DATALOC.FineScaleDomainsToShow  ;
DATALOC = DefaultField(DATALOC,'ReconstructionStressesUsingCoarseScaleDISP',1) ;
DATA.ReconstructionStressesUsingCoarseScaleDISP = DATALOC.ReconstructionStressesUsingCoarseScaleDISP  ;

%OTHERinputs.DATAinpGID = DATAinpGID;

strainCOARSE_history = PostProcess_EIFE_vect(COOR,CN, PROPMAT,MaterialType,...
    TypeElement,DATA,MESH.TRANSF_COORD,d,MESH.Vrot,stressGLO_REF,VON_mises_COARSE,InternalVarStrain,OTHERinputs);





%
%
% NAME_INPUT_DATA = 'MESH';
% IND_ELEM_MESHES = GidMesh2DFE_multi(NameFile_msh,COOR,CN,NAME_INPUT_DATA,MaterialType,TypeElement,DATA.NAMEMESH);
%
%
%
% % Writing results file
% % Nodal vectors
% % -------------------
% NODESV_PROP = DATA.PRINT.NODESV_PROP ;
% NODAL_VECTOR = [] ;
% NAMES_VAR = fieldnames(NODESV_PROP) ;
% ivar_vect = 0 ;
% ivar_scalar = 0 ;
% NODAL_SCALAR = [];
% for iname = 1:length(NAMES_VAR)
%     locname = NAMES_VAR{iname} ;
%     LOCFF= NODESV_PROP.(locname) ;
%     LOCFF = DefaultField(LOCFF,'PRINT',0) ;
%     if LOCFF.PRINT == 1 && ~isempty(SNAP.(locname))
%
%         % Components
%         COMP = NODESV_PROP.(locname).COMP ;
%         LEGEND = NODESV_PROP.(locname).LEGEND ;
%         TYPE = NODESV_PROP.(locname).TYPE ;
%
%         % Extract data
%         % ---------------
%
%
%         switch TYPE
%             case 'Vector'
%                 ivar_vect = ivar_vect +1 ;
%                 NODAL_VECTOR(ivar_vect).NAME = LEGEND ;
%                 NODAL_VECTOR(ivar_vect).COMP = COMP;
%                 NODAL_VECTOR(ivar_vect).VAR = SNAP.(locname);
%                 NODAL_VECTOR(ivar_vect).NODES = 1:size(COOR,1) ;
%                 NODAL_VECTOR(ivar_vect).NAME_SNAP = locname ;
%             case 'Scalar'
%                 ivar_scalar = ivar_scalar +1 ;
%                 NODAL_SCALAR(ivar_scalar).NAME = LEGEND ;
%                 NODAL_SCALAR(ivar_scalar).COMP = COMP;
%                 NODAL_SCALAR(ivar_scalar).VAR =  SNAP.(locname);
%                 NODAL_SCALAR(ivar_scalar).NODES = 1:size(COOR,1) ;
%                 NODAL_SCALAR(ivar_scalar).NAME_SNAP = locname ;
%
%         end
%
%     end
%
% end
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MESH = [] ;
%
% for imesh = 1:length(CN)
%     MESH(imesh).TypeElement = TypeElement{imesh} ;
%     MESH(imesh).posgp = posgp{imesh} ;
%     MESH(imesh).NAMEMESH = DATA.NAMEMESH{imesh};
%
%     % -------------------
%     NAMES_VAR = fieldnames(GAUSSV_PROP{imesh}) ;
%     ivar_matrix = 0 ;
%     ivar_scalar = 0 ;
%     ivar_vector = 0 ;
%     for iname = 1:length(NAMES_VAR)
%         locname = NAMES_VAR{iname} ;
%
%         GAUSSV_PROP{imesh}.(locname) = DefaultField(GAUSSV_PROP{imesh}.(locname),'PRINT',0) ;
%
%         if GAUSSV_PROP{imesh}.(locname).PRINT == 1
%
%             % Components
%             LEGEND = GAUSSV_PROP{imesh}.(locname).LEGEND ;
%             TYPE = GAUSSV_PROP{imesh}.(locname).TYPE ;
%             var= SNAP.(locname) ;
%
%             % RECONSTRUCTION
%             DATAinpGID = DefaultField(DATAinpGID,'OPERreconstr',[]);
%             DATAinpGID.OPERreconstr = DefaultField(DATAinpGID.OPERreconstr,locname,[]) ;
%             if ~isempty(DATAinpGID.OPERreconstr.(locname))
%                 var = DATAinpGID.OPERreconstr.(locname).coeff*var;
%                 var = DATAinpGID.OPERreconstr.(locname).BASIS*var ;
%             end
%
%
%             switch TYPE
%                 case 'Matrix'
%                     COMP = GAUSSV_PROP{imesh}.(locname).COMP ;
%
%                     ivar_matrix = ivar_matrix +1 ;
%
%                     INDICES=  GAUSSV_PROP{imesh}.(locname).CONVERSION_GID ;
%                     %   [AAAA,INDICESinv ]= sort(INDICES) ;
%
%                     nstrainACTUAL = length(INDICES) ;
%                     DATA.MESH.nstrainACTUAL = nstrainACTUAL ;
%                     if DATA.MESH.nstrain ~= nstrainACTUAL
%                         % This happens in plane-stress !!!
%                         nelemGAUSS=  size(var,1)/DATA.MESH.nstrain ;
%                         nnew = nelemGAUSS*length(INDICES);
%                         varNEW  = zeros(nnew,size(var,2)) ;
%
%                         for istress = 1:DATA.MESH.nstrain
%                             %   if istress <=DATA.MESH.nstrain
%                             indOLD = istress:DATA.MESH.nstrain:size(var,1) ;
%                             istressNEW =  find(INDICES ==istress) ;
%                             indNEW = istressNEW:nstrainACTUAL:size(varNEW,1) ;
%                             varNEW(indNEW,:) = var(indOLD,:) ;
%                             %   end
%                         end
%                     else
%                         varNEW  = zeros(size(var)) ;
%
%                         for istress = 1:length(INDICES)
%                             indOLD = istress:DATA.MESH.nstrain:size(var,1) ;
%                             indNEW = INDICES(istress):DATA.MESH.nstrain:size(var,1) ;
%                             varNEW(indNEW,:) = var(indOLD,:) ;
%                         end
%
%                     end
%
%
%                     % Conversion to GID indidces
%
%
%                     MESH(imesh).GAUSS_MATRIX(ivar_matrix).NAME = LEGEND ;
%                     MESH(imesh).GAUSS_MATRIX(ivar_matrix).COMP = COMP;
%                     MESH(imesh).GAUSS_MATRIX(ivar_matrix).VAR = varNEW ;
%                     MESH(imesh).GAUSS_MATRIX(ivar_matrix).ELEMENTS = IND_ELEM_MESHES{imesh} ;
%                 case 'Scalar'
%                     ivar_scalar = ivar_scalar +1 ;
%                     MESH(imesh).GAUSS_SCALAR(ivar_scalar).NAME = LEGEND ;
%                     MESH(imesh).GAUSS_SCALAR(ivar_scalar).VAR =var ;
%                     MESH(imesh).GAUSS_SCALAR(ivar_scalar).ELEMENTS = IND_ELEM_MESHES{imesh} ; ;
%
%                 case 'Vector'
%                     COMP = GAUSSV_PROP{imesh}.(locname).COMP ;
%                     ivar_vector = ivar_vector +1 ;
%                     MESH(imesh).GAUSS_VECTOR(ivar_vector).NAME = LEGEND ;
%                     MESH(imesh).GAUSS_VECTOR(ivar_vector).COMP = COMP;
%                     MESH(imesh).GAUSS_VECTOR(ivar_vector).VAR = var ;
%                     MESH(imesh).GAUSS_VECTOR(ivar_vector).ELEMENTS = IND_ELEM_MESHES{imesh} ; ;
%
%             end
%
%         end
%
%     end
%
%
% end
%
% DATA.TIME_DISCRETIZATION = TIME_PRINT ;
% GidResults2DFE_multi_LARGE(NameFile_res,DATA.MESH.ndim,NODAL_VECTOR,NODAL_SCALAR,MESH,DATA,DATAinpGID);
%
% cddd = cd ;
% NAMEFILEOPEN =  [NameFile_res] ;
%
%
% %  if DATA.COMPRESS_GID == 1
% %     disp('Compressing GID FILE ...')
% %     TEXT = ['gid -PostResultsToBinary  ',NAMEFILEOPEN, ' ',NAMEFILEOPEN] ;
% %     unix(TEXT) ;
% %     delete(NameFile_msh)
% %     disp('Done')
% % end
%
%
% disp('open GID FILE:')
% %clipboard('copy',NAMEFILEOPEN)
% disp(NAMEFILEOPEN)
%
% DATA = DefaultField(DATA,'OPEN_GID',0) ;
% if DATA.OPEN_GID ==1
%     TTTT = ['gidpost ',NAMEFILEOPEN] ;
%     unix(TTTT);
% end