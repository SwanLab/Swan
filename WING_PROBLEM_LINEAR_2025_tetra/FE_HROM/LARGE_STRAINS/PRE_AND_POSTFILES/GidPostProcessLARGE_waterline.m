function GidPostProcessLARGE_waterline(MESH,DATA,icluster,DATAinpGID)

if nargin==0
    load('tmp.mat')
    DATAinpGID.IS_HROM_SIMULATION = 1; 
elseif nargin == 3
    DATAinpGID = [] ;
    
end

% COORDINATES 
% -----------
nnodeSOLID = size(MESH.COOR,1) ; % Number of Nodes of the solid mesh 
COOR = [MESH.COOR; DATA.WATERLINE_MESH.COOR] ;  % Coordinate matrix 
CNwaterline = DATA.WATERLINE_MESH.CN + nnodeSOLID ;  % Connectivities free surface (waterline)
nmat = max(MESH.MaterialType) ;
MaterialTypeRED=  ones(size(MESH.MaterialType)) + nmat ;  
MaterialTypeWS = ones(size(MESH.HYDRO.CNb,1),1)*(nmat+1) ;
MaterialTypeWL = ones(size(CNwaterline,1),1)*(nmat+2) ;

%  GAUSS_SNAP = DATA.STORE.NAME_MATFILE_STORE{icluster};


DATAinpGID = DefaultField(DATAinpGID,'IS_HROM_SIMULATION',0) ; 
DATAinpGID.NOSELECTION = [] ; 

if DATAinpGID.IS_HROM_SIMULATION == 0
    % No reduced meshes
    ECM_data_STRING = {'NOSELECTION','NOSELECTION','NOSELECTION'} ;
    CNglo = {MESH.CN,MESH.HYDRO.CNb,CNwaterline} ;
    DATA.NAMEMESH = {'INTERIOR','WETSURFACES','WATERLINE'} ;
    
    TypeElement = {MESH.TypeElement,MESH.TypeElementB,DATA.WATERLINE_MESH.TypeElement} ;
    MaterialTypeGLO = {MESH.MaterialType,MaterialTypeWS,MaterialTypeWL} ;
    posgp = {DATA.MESH.posgp,[],[]} ;
    GAUSSV_PROP = {DATA.PRINT.GAUSSV_PROP,[],[]} ;


    
else
    ECM_data_STRING = {'NOSELECTION','ECMdata','ECMdata_press','NOSELECTION'} ;
    CNglo = {MESH.CN,MESH.CN,MESH.HYDRO.CNb,CNwaterline} ;
    DATA.NAMEMESH = {'INTERIOR','INTERIOR_RED','WETSURFACES','WATERLINE'} ;
    TypeElement = {MESH.TypeElement,MESH.TypeElement,MESH.TypeElementB,DATA.WATERLINE_MESH.TypeElement} ;
    MaterialTypeGLO = {MESH.MaterialType,MaterialTypeRED,MaterialTypeWS,MaterialTypeWL} ;
    posgp = {DATA.MESH.posgp,DATA.MESH.posgp,[],[]} ;
     GAUSSV_PROP = {DATA.PRINT.GAUSSV_PROP,DATA.PRINT.GAUSSV_PROP,[],[]} ;
end


setElements = cell(size(ECM_data_STRING)) ; setGaussPoints = setElements ; setIndexSTRESSES = setElements ;
setIndexVECTOR = setElements ;
CN = setElements ; % {CN,MESH.HYDRO.CNb(setElements{imesh},:),CNwaterline} ;
MaterialType = MaterialTypeGLO ; 

for imesh = 1:length(ECM_data_STRING)
    DATAinpGID = DefaultField(DATAinpGID,ECM_data_STRING{imesh},[]) ;
    if ~isempty(DATAinpGID.(ECM_data_STRING{imesh}))
        % ONLY SHOWING HROM MESH
        % ---------------------------------
        setElements{imesh} = DATAinpGID.(ECM_data_STRING{imesh}).setElements ;
    else
        % if  imesh == 1
        setElements{imesh} = 1:size(CNglo{imesh},1) ;
        % elseif imesh == 2
        %     setElements{imesh} = 1:size(MESH.HYDRO.CNb,1) ;
        %  end
    end
    setGaussPoints{imesh} = small2large(setElements{imesh},DATA.MESH.ngaus) ;
    setIndexSTRESSES{imesh}= small2large(setGaussPoints{imesh},DATA.MESH.nstrain) ;
    setIndexVECTOR{imesh}= small2large(setGaussPoints{imesh},DATA.MESH.ndim) ;
    CN{imesh} = CNglo{imesh}(setElements{imesh},:) ; 
    MaterialType{imesh} = MaterialTypeGLO{imesh}(setElements{imesh},:) ; 
end
  


 

 

NameFile_msh = DATA.PRINT.NAME_FILE_MSH{icluster} ;
NameFile_res= DATA.PRINT.NAME_FILE_RES{icluster} ;

 




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
                NODAL_VECTOR(ivar_vect).NODES = 1:size(MESH.COOR,1) ;
                NODAL_VECTOR(ivar_vect).NAME_SNAP = locname ;
            case 'Scalar'
                ivar_scalar = ivar_scalar +1 ;
                NODAL_SCALAR(ivar_scalar).NAME = LEGEND ;
                NODAL_SCALAR(ivar_scalar).COMP = COMP;
                NODAL_SCALAR(ivar_scalar).VAR =  SNAP.(locname);
                NODAL_SCALAR(ivar_scalar).NODES = 1:size(MESH.COOR,1) ;
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
    if  ~isempty(GAUSSV_PROP{imesh})
        
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
                DATAinpGID = DefaultField(DATAinpGID,'OPERreconstr',[]);
                DATAinpGID.OPERreconstr = DefaultField(DATAinpGID.OPERreconstr,locname,[]) ;
                if ~isempty(DATAinpGID.OPERreconstr.(locname))
                    var = DATAinpGID.OPERreconstr.(locname).coeff*var;
                    var = DATAinpGID.OPERreconstr.(locname).BASIS*var ;
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
                        MESH(imesh).GAUSS_MATRIX(ivar_matrix).VAR = varNEW(setIndexSTRESSES{imesh},:) ;
                        MESH(imesh).GAUSS_MATRIX(ivar_matrix).ELEMENTS = IND_ELEM_MESHES{imesh} ;
                    case 'Scalar'
                        ivar_scalar = ivar_scalar +1 ;
                        MESH(imesh).GAUSS_SCALAR(ivar_scalar).NAME = LEGEND ;
                        MESH(imesh).GAUSS_SCALAR(ivar_scalar).VAR =var(setGaussPoints{imesh},:) ;
                        MESH(imesh).GAUSS_SCALAR(ivar_scalar).ELEMENTS = IND_ELEM_MESHES{imesh} ; 
                        
%                         setGaussPoints = small2large(setElements,DATA.MESH.ngaus) ;
% setIndexSTRESSES= small2large(setGaussPoints,DATA.MESH.nstrain) ;
% setIndexVECTOR= small2large(setGaussPoints,DATA.MESH.ndim) ;

                        
                    case 'Vector'
                        COMP = GAUSSV_PROP{imesh}.(locname).COMP ;
                        ivar_vector = ivar_vector +1 ;
                        MESH(imesh).GAUSS_VECTOR(ivar_vector).NAME = LEGEND ;
                        MESH(imesh).GAUSS_VECTOR(ivar_vector).COMP = COMP;
                        MESH(imesh).GAUSS_VECTOR(ivar_vector).VAR = var(setIndexVECTOR{imesh},:) ;
                        MESH(imesh).GAUSS_VECTOR(ivar_vector).ELEMENTS = IND_ELEM_MESHES{imesh} ; ;
                        
                end
                
            end
            
        end
    end
    
    
end

DATA.TIME_DISCRETIZATION = TIME_PRINT ;
GidResults2DFE_multi_LARGE(NameFile_res,DATA.MESH.ndim,NODAL_VECTOR,NODAL_SCALAR,MESH,DATA,DATAinpGID);

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