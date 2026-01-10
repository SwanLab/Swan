function  strainCOARSE_history = GidPostProcess_VECT_EIFEbub(MESH,DATA,icluster,DATAinpGID,PROPMAT,DATALOC)

if nargin==0
    load('tmp.mat')
elseif nargin == 3
    DATAinpGID = [] ;
    
end
COOR = MESH.COOR;
MESH = DefaultField(MESH,'AUXILIAR',[]); 
if isempty(MESH.AUXILIAR)
    CN = MESH.CN ;
    TypeElement = MESH.TypeElement;
    MaterialType_coarse_GID = MESH.MaterialType;
    
else
    % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/08_GENERAL_MODES.mlx
    CN = MESH.AUXILIAR.CN  ;
    TypeElement = MESH.AUXILIAR.TypeElement  ;
    MaterialType_coarse_GID = ones(size(CN,1),1) ;
end
MaterialType = MESH.MaterialType;


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
d = SNAP.DISP;%(:,STEPS_PRINT) ;

OTHERinputs.COORprintCOARSE = [] ; OTHERinputs.CNprintCOARSE = [] ; OTHERinputs.dPRINTcoarse = [] ;

if isfield(MESH,'SMOOTH_from_UNCOUP_TO_SUPPORT')
    OTHERinputs.dPRINTcoarse = MESH.SMOOTH_from_UNCOUP_TO_SUPPORT*d ; 
    % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/06_UNCOUPLED.mlx
    OTHERinputs.COORprintCOARSE =  MESH.COOR_SUPPORT ; 
    OTHERinputs.CNprintCOARSE = MESH.CN_SUPPORT ; 
 
end



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
  
end

if isfield(SNAP,'PK2STRESS')
stressGLO_REF = SNAP.PK2STRESS; %(:,STEPS_PRINT);  % These are the stresses in the unrotated configuration
else
    stressGLO_REF = [] ; 
end
DATA.NameFile_msh = NameFile_msh;
DATA.NameFile_res = NameFile_res;
SNAP = DefaultField(SNAP,'VONMISES_CAUCHY_STRESS',[]) ; 

if ~isempty(SNAP.VONMISES_CAUCHY_STRESS) && isempty(MESH.AUXILIAR)
    VON_mises_COARSE= zeros(size(CN,1),size(SNAP.VONMISES_CAUCHY_STRESS,2));
    for itime = 1:size(SNAP.VONMISES_CAUCHY_STRESS,2)
        VONloc  = reshape(SNAP.VONMISES_CAUCHY_STRESS(:,itime),DATA.MESH.ngaus_STRESS,[]) ;
        VON_mises_COARSE(:,itime) = max(VONloc,[],1)  ;
    end
    VON_mises_COARSE = VON_mises_COARSE; %(:,STEPS_PRINT);
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
    InternalVarStrain = InternalVarStrain(:,STEPS_PRINT) ;
else
    InternalVarStrain = [] ;
end





DATALOC = DefaultField(DATALOC,'FineScaleDomainsToShow',[])  ;
DATA.FineScaleDomainsToShow = DATALOC.FineScaleDomainsToShow  ;
DATALOC = DefaultField(DATALOC,'ReconstructionStressesUsingCoarseScaleDISP',1) ;
DATA.ReconstructionStressesUsingCoarseScaleDISP = DATALOC.ReconstructionStressesUsingCoarseScaleDISP  ;
OTHERinputs.MaterialType_coarse_GID = MaterialType_coarse_GID; 
strainCOARSE_history = PostProcess_EIFE_vectBUB(COOR,CN, PROPMAT,MaterialType,...
    TypeElement,DATA,MESH.TRANSF_COORD,d,MESH.Vrot,stressGLO_REF,VON_mises_COARSE,InternalVarStrain,OTHERinputs);


