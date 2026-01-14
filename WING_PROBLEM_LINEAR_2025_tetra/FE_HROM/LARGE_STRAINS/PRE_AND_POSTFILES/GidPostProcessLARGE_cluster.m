function GidPostProcessLARGE_cluster(MESH,DATA,icluster,DATAinpGID,CLUSTER_SEQUENCE)

if nargin==0
    load('tmp.mat')
elseif nargin == 3
    DATAinpGID = [] ;
    
end
COOR = MESH.COOR;
CN = MESH.CN ;
TypeElement = MESH.TypeElement;
MaterialType = MESH.MaterialType;


NameFile_msh = DATA.PRINT.NAME_FILE_MSH{icluster} ;
NameFile_res= DATA.PRINT.NAME_FILE_RES{icluster} ;

% For domain decomposition purpooses (not useful now )
DATA = DefaultField(DATA,'MATERIAL_ORIGINAL',[]) ;
DATA.MATERIAL =  DATA.MATERIAL_ORIGINAL ;
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
    GAUSSV_PROP = {DATA.PRINT.GAUSSV_PROP} ;
    %  GAUSS_SNAP = DATA.STORE.NAME_MATFILE_STORE{icluster};
    posgp = {DATA.MESH.posgp} ;
else
    error('Option not implemented...9-Dec-2020')
    GAUSSV_PROP = DATA.PRINT.GAUSSV_PROP ;
    %  GAUSS_SNAP = DATA.PRINT.GAUSS_SNAP ;
    posgp  = DATA.MESH.posgp;
end
if exist(DATA.STORE.NAME_MATFILE_STORE{icluster})
    load( DATA.STORE.NAME_MATFILE_STORE{icluster},'SNAP_cluster','STEP_LOC') ;
else
    return
end


NAME_INPUT_DATA = 'MESH';
IND_ELEM_MESHES = GidMesh2DFE_multi(NameFile_msh,COOR,CN,NAME_INPUT_DATA,MaterialType,TypeElement,DATA.NAMEMESH);


[STEPS_PRINT,iaaa,ibbb ]= intersect(STEP_LOC,DATA.PRINT.NSTEPS) ;
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
TIME_PRINT = DATA.STEPS(STEPS_PRINT) ;


% Writing results file
% Nodal vectors
% -------------------
NODESV_PROP = DATA.PRINT.NODESV_PROP ;
NODAL_VECTOR = [] ;
NAMES_VAR = fieldnames(NODESV_PROP) ;
ivar_vect = 0 ;
ivar_scalar = 0 ;
NODAL_SCALAR = [];
for iname = 1:length(NAMES_VAR)
    locname = NAMES_VAR{iname} ;
    LOCFF= NODESV_PROP.(locname) ;
    LOCFF = DefaultField(LOCFF,'PRINT',0) ;
    if LOCFF.PRINT == 1 && ~isempty(SNAP.(locname))
        
        % Components
        COMP = NODESV_PROP.(locname).COMP ;
        LEGEND = NODESV_PROP.(locname).LEGEND ;
        TYPE = NODESV_PROP.(locname).TYPE ;
        
        % Extract data
        % ---------------
        
        
        switch TYPE
            case 'Vector'
                ivar_vect = ivar_vect +1 ;
                NODAL_VECTOR(ivar_vect).NAME = LEGEND ;
                NODAL_VECTOR(ivar_vect).COMP = COMP;
                NODAL_VECTOR(ivar_vect).VAR = SNAP.(locname);
                NODAL_VECTOR(ivar_vect).NODES = 1:size(COOR,1) ;
                NODAL_VECTOR(ivar_vect).NAME_SNAP = locname ;
            case 'Scalar'
                ivar_scalar = ivar_scalar +1 ;
                NODAL_SCALAR(ivar_scalar).NAME = LEGEND ;
                NODAL_SCALAR(ivar_scalar).COMP = COMP;
                NODAL_SCALAR(ivar_scalar).VAR =  SNAP.(locname);
                NODAL_SCALAR(ivar_scalar).NODES = 1:size(COOR,1) ;
                NODAL_SCALAR(ivar_scalar).NAME_SNAP = locname ;
                
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
    ivar_vector = 0 ;
    for iname = 1:length(NAMES_VAR)
        locname = NAMES_VAR{iname} ;
        
        
        
        if GAUSSV_PROP{imesh}.(locname).PRINT == 1
            
            % Components
            LEGEND = GAUSSV_PROP{imesh}.(locname).LEGEND ;
            TYPE = GAUSSV_PROP{imesh}.(locname).TYPE ;
            var= SNAP.(locname) ;
            
            % RECONSTRUCTION
            %  DATAinpGID = DefaultField(DATAinpGID,'OPERreconstr',[]);
            %  DATAinpGID.OPERreconstr = DefaultField(DATAinpGID.OPERreconstr,locname,[]) ;
            OPERreconst = DATAinpGID.OPERreconstr.(locname) ; 
            
            if isstruct(OPERreconst.BASIS)
                 nrows = size(DATAinpGID.OPERreconstr.(locname).BASIS.BASIS,2) ;
            else
                nrows = size(DATAinpGID.OPERreconstr.(locname).BASIS{1},1) ;
            end
            
            
            varNEW = zeros(nrows,size(var,2)) ;
            for  iprintt = 1:length(STEPS_PRINT)
                itimeLOC = STEPS_PRINT(iprintt) ;
                iCL =  CLUSTER_SEQUENCE(itimeLOC) ;
                %           if ~isempty(DATAinpGID.OPERreconstr.(locname))
                nrowsCOEFF = size(OPERreconst.coeff{iCL},2) ;
                varLOC = OPERreconst.coeff{iCL}*var(1:nrowsCOEFF,iprintt);
                if isstruct(OPERreconst.BASIS)
                    % varLOC = OPERreconst.BASIS.BASIS*OPERreconst.BASIS.coeff{iCL}*varLOC ;
                     varLOC =OPERreconst.BASIS.coeff{iCL}*varLOC ;
                else
                varLOC = OPERreconst.BASIS{iCL}*varLOC ;
                end
                varNEW(:,iprintt)  = varLOC ;
                %          end
            end
            
            % 
            
            if isstruct(OPERreconst.BASIS)
                var =   OPERreconst.BASIS.BASIS*varNEW ;
            else
                
                var = varNEW ;
            end
            
            
            switch TYPE
                case 'Matrix'
                    COMP = GAUSSV_PROP{imesh}.(locname).COMP ;
                    
                    ivar_matrix = ivar_matrix +1 ;
                    
                    INDICES=  GAUSSV_PROP{imesh}.(locname).CONVERSION_GID ;
                    %   [AAAA,INDICESinv ]= sort(INDICES) ;
                    
                    nstrainACTUAL = length(INDICES) ;
                    DATA.MESH.nstrainACTUAL = nstrainACTUAL ;
                    if DATA.MESH.nstrain ~= nstrainACTUAL
                        % This happens in plane-stress !!!
                        nelemGAUSS=  size(var,1)/DATA.MESH.nstrain ;
                        nnew = nelemGAUSS*length(INDICES);
                        varNEW  = zeros(nnew,size(var,2)) ;
                        
                        for istress = 1:DATA.MESH.nstrain
                            %   if istress <=DATA.MESH.nstrain
                            indOLD = istress:DATA.MESH.nstrain:size(var,1) ;
                            istressNEW =  find(INDICES ==istress) ;
                            indNEW = istressNEW:nstrainACTUAL:size(varNEW,1) ;
                            varNEW(indNEW,:) = var(indOLD,:) ;
                            %   end
                        end
                    else
                        varNEW  = zeros(size(var)) ;
                        
                        for istress = 1:length(INDICES)
                            indOLD = istress:DATA.MESH.nstrain:size(var,1) ;
                            indNEW = INDICES(istress):DATA.MESH.nstrain:size(var,1) ;
                            varNEW(indNEW,:) = var(indOLD,:) ;
                        end
                        
                    end
                    
                    
                    % Conversion to GID indidces
                    
                    
                    MESH(imesh).GAUSS_MATRIX(ivar_matrix).NAME = LEGEND ;
                    MESH(imesh).GAUSS_MATRIX(ivar_matrix).COMP = COMP;
                    MESH(imesh).GAUSS_MATRIX(ivar_matrix).VAR = varNEW ;
                    MESH(imesh).GAUSS_MATRIX(ivar_matrix).ELEMENTS = IND_ELEM_MESHES{imesh} ;
                case 'Scalar'
                    ivar_scalar = ivar_scalar +1 ;
                    MESH(imesh).GAUSS_SCALAR(ivar_scalar).NAME = LEGEND ;
                    MESH(imesh).GAUSS_SCALAR(ivar_scalar).VAR =var ;
                    MESH(imesh).GAUSS_SCALAR(ivar_scalar).ELEMENTS = IND_ELEM_MESHES{imesh} ; ;
                    
                case 'Vector'
                    COMP = GAUSSV_PROP{imesh}.(locname).COMP ;
                    ivar_vector = ivar_vector +1 ;
                    MESH(imesh).GAUSS_VECTOR(ivar_vector).NAME = LEGEND ;
                    MESH(imesh).GAUSS_VECTOR(ivar_vector).COMP = COMP;
                    MESH(imesh).GAUSS_VECTOR(ivar_vector).VAR = var ;
                    MESH(imesh).GAUSS_VECTOR(ivar_vector).ELEMENTS = IND_ELEM_MESHES{imesh} ; ;
                    
            end
            
        end
        
    end
    
    
end

DATA.TIME_DISCRETIZATION = TIME_PRINT ;
GidResults2DFE_multi_LARGE_cluster(NameFile_res,DATA.MESH.ndim,NODAL_VECTOR,NODAL_SCALAR,MESH,DATA,DATAinpGID,CLUSTER_SEQUENCE);

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
%clipboard('copy',NAMEFILEOPEN)
disp(NAMEFILEOPEN)

DATA = DefaultField(DATA,'OPEN_GID',0) ;
if DATA.OPEN_GID ==1
    TTTT = ['gidpost ',NAMEFILEOPEN] ;
    unix(TTTT);
end