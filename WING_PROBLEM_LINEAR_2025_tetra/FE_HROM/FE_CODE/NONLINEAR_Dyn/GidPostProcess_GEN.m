function GidPostProcess_GEN(COOR,CN,TypeElement,DATA,NAME_BASE_GIDfiles,NameFileMesh,MaterialType,...
    NODES_SNAP,GAUSS_SNAP,OPERfe,NODESV_PROP,GAUSSV_PROP);
% Post-processing of results using GID
%dbstop('5')
if nargin==0
    load('tmp1.mat')
end

% Name of the mesh file
%if ~isempty(NameFileMesh)
%[dummy1 NameFileMeshHERE dummy2]= fileparts(NameFileMesh) ;
%else
NameFileMeshHERE = '' ;
%end

DATA = DefaultField(DATA,'POST_LabelFileGid','') ;
[NAME_INPUT_DATA_base ,NAME_INPUT_DATA]= fileparts(NAME_BASE_GIDfiles) ;

NameFile_msh = [NAME_INPUT_DATA_base,filesep,'GIDPOST',filesep,NameFileMeshHERE,NAME_INPUT_DATA,DATA.POST_LabelFileGid ,'.msh'] ;
% Name of the results file
NameFile_res= [NAME_INPUT_DATA_base,filesep,'GIDPOST',filesep,NameFileMeshHERE,NAME_INPUT_DATA,DATA.POST_LabelFileGid,'.res'] ;

if ~exist([NAME_INPUT_DATA_base,filesep,'GIDPOST'],'dir')
    mkdir([NAME_INPUT_DATA_base,filesep,'GIDPOST'])
end

%error('Fix this')
DATA.MATERIAL = DATA.MATERIAL_ORIGINAL ; 
DATA = DefaultField(DATA,'PostProcessWithNOSLICES',1) ; 
if DATA.PostProcessWithNOSLICES == 1 && ~isempty(DATA.MATERIAL) && (~isempty(DATA.MakeMeshByRepetition.nDOM) && DATA.MakeMeshByRepetition.nDOM(1)>1 )
    nmat = length(DATA.MATERIAL.PLY) ;
    ndom = prod(DATA.MakeMeshByRepetition.nDOM); %(1) ;
    TOTnmat = nmat*ndom ;
    MAT_TypeORIG = (1:TOTnmat) ;
    MAT_TypeNOSLICES = repmat((1:nmat)',ndom,1) ;
    NewMaterialType = zeros(size(MaterialType));
    
    for imat = 1:length(MAT_TypeORIG)
        INDLOC =find(MaterialType==imat) ;
        NewMaterialType(INDLOC) = MAT_TypeNOSLICES(imat) ;
    end
    
    MaterialType =  NewMaterialType ;
end


if ~iscell(CN)
    DATA.NAMEMESH = {'MESH'} ;
    CN = {CN} ;
    MaterialType = {MaterialType} ;
    TypeElement = {TypeElement} ;
    GAUSSV_PROP = {GAUSSV_PROP} ;
    GAUSS_SNAP = {GAUSS_SNAP} ;
    posgp = {DATA.posgp} ;
else
    posgp  = DATA.posgp;
end
IND_ELEM_MESHES = GidMesh2DFE_multi(NameFile_msh,COOR,CN,NAME_INPUT_DATA,MaterialType,TypeElement,DATA.NAMEMESH);


% Writing results file
% Nodal vectors
% -------------------
NODAL_VECTOR = [] ;
NAMES_VAR = fieldnames(NODESV_PROP) ;
ivar_vect = 0 ;
ivar_scalar = 0 ;
NODAL_SCALAR = []; 
for iname = 1:length(NAMES_VAR)
    locname = NAMES_VAR{iname} ;
    if NODESV_PROP.(locname).PRINT == 1
        
        % Components
        COMP = NODESV_PROP.(locname).COMP ;
        LEGEND = NODESV_PROP.(locname).LEGEND ;
        TYPE = NODESV_PROP.(locname).TYPE ;
        
        switch TYPE
            case 'Vector'
                ivar_vect = ivar_vect +1 ;
                NODAL_VECTOR(ivar_vect).NAME = LEGEND ;
                NODAL_VECTOR(ivar_vect).COMP = COMP;
                NODAL_VECTOR(ivar_vect).VAR = NODES_SNAP.(locname) ;
                NODAL_VECTOR(ivar_vect).NODES = 1:size(COOR,1) ;
            case 'Scalar'
                ivar_scalar = ivar_scalar +1 ;
                NODAL_SCALAR(ivar_scalar).NAME = LEGEND ;
                NODAL_SCALAR(ivar_scalar).COMP = COMP;
                NODAL_SCALAR(ivar_scalar).VAR = NODES_SNAP.(locname) ;
                NODAL_SCALAR(ivar_scalar).NODES = 1:size(COOR,1) ;
                
        end
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MESH = [] ;

for imesh = 1:length(CN)
    MESH(imesh).TypeElement = TypeElement{imesh} ;
    MESH(imesh).posgp = posgp{imesh} ;
    MESH(imesh).NAMEMESH = DATA.NAMEMESH{imesh};
    
    % -------------------
    NAMES_VAR = fieldnames(GAUSSV_PROP{imesh}) ;
    ivar_matrix = 0 ;
    ivar_scalar = 0 ;
    for iname = 1:length(NAMES_VAR)
        locname = NAMES_VAR{iname} ;
        if GAUSSV_PROP{imesh}.(locname).PRINT == 1
            
            % Components
            LEGEND = GAUSSV_PROP{imesh}.(locname).LEGEND ;
            TYPE = GAUSSV_PROP{imesh}.(locname).TYPE ;
            
            switch TYPE
                case 'Matrix'
                    COMP = GAUSSV_PROP{imesh}.(locname).COMP ;
                    
                    ivar_matrix = ivar_matrix +1 ;
                    
                    INDICES=  GAUSSV_PROP{imesh}.(locname).CONVERSION_GID ;
                    var= GAUSS_SNAP{imesh}.(locname) ;
                    varNEW  = zeros(size(var)) ;
                    % Conversion to GID indidces
                    for istress = 1:length(INDICES)
                        indOLD = istress:OPERfe.nstrain:size(var,1) ;
                        indNEW = INDICES(istress):OPERfe.nstrain:size(var,1) ;
                        varNEW(indNEW,:) = var(indOLD,:) ;
                    end
                    
                    MESH(imesh).GAUSS_MATRIX(ivar_matrix).NAME = LEGEND ;
                    MESH(imesh).GAUSS_MATRIX(ivar_matrix).COMP = COMP;
                    MESH(imesh).GAUSS_MATRIX(ivar_matrix).VAR = varNEW ;
                    MESH(imesh).GAUSS_MATRIX(ivar_matrix).ELEMENTS = IND_ELEM_MESHES{imesh} ;
                case 'Scalar'
                    ivar_scalar = ivar_scalar +1 ;
                    MESH(imesh).GAUSS_SCALAR(ivar_scalar).NAME = LEGEND ;
                    MESH(imesh).GAUSS_SCALAR(ivar_scalar).VAR = GAUSS_SNAP{imesh}.(locname) ;
                    MESH(imesh).GAUSS_SCALAR(ivar_scalar).ELEMENTS = IND_ELEM_MESHES{imesh} ; ;
                    
            end
            
        end
        
    end
    
    
end


GidResults2DFE_multi_NEW(NameFile_res,OPERfe.ndim,NODAL_VECTOR,NODAL_SCALAR,MESH,DATA,OPERfe);

cddd = cd ;
NAMEFILEOPEN =  [NameFile_res] ;


%  if DATA.COMPRESS_GID == 1
%     disp('Compressing GID FILE ...')
%     TEXT = ['gid -PostResultsToBinary  ',NAMEFILEOPEN, ' ',NAMEFILEOPEN] ;
%     unix(TEXT) ;
%     delete(NameFile_msh)
%     disp('Done')
% end


disp('open GID FILE:')
clipboard('copy',NAMEFILEOPEN)
disp(NAMEFILEOPEN)

DATA = DefaultField(DATA,'OPEN_GID',0) ;
if DATA.OPEN_GID ==1
    TTTT = ['gidpost ',NAMEFILEOPEN] ;
    unix(TTTT);
end