function [NameFile_msh,NameFile_res,COORprint,CNprint,NAME_INPUT_DATA,...
    TypeElementPRINT,MaterialTypegloPRINT,NAMEMESH,DATAIN,DATAmeshR]=...
    Construct3DstructureROMrve(MESH2D,MESH3D,DATAIN,MESH_LATERAL,FOLDER,DATA_REFMESH)
% Plot 3D domains along with 1D skeleton
% JAHO, 22 Jan-2018
if nargin == 0
    load('tmp1.mat')
end

% --------------------------------------
% 1D mesh ------------------------------
% --------------------------------------
iacumNODES =  size(MESH2D.COORall,1) ; % Connectivities index (cumulative, for printing global 3D mesh)
COORprint = [MESH2D.COORall] ;  % Coordinates global GID mesh
CNprint =  {MESH2D.CNall};  % Connectivities global GID mesh
TypeElementPRINT=  {MESH2D.TypeElement} ;  % Type of element global GID mesh
MaterialTypegloPRINT = {MESH2D.MaterialType} ; % {ones(size(MESH2D.CN,1),1)} ; % Type of element global GID mesh
NAMEMESH =   {'2D'} ;
iMESH = 1;  % Number of meshes constructed so far
%imatACUM = 1 ;  % Number of GID materials defined so far
imatACUM = length(unique(MESH2D.MaterialType)) ;

% For lateral meshes
CNprint_LAT = {} ;
TypeElementPRINT_LAT = {} ;
MaterialTypegloPRINT_LAT = {} ;
NAMEMESH_LAT = {} ;
iMESH_LAT = 0 ;


 %save(DATAIN.nameGRAPHS,'DATAmeshR')  ; % Used for reconstructing 3D meshes with GID

 
 

for iENTITY  = 1: length(MESH2D.PROP) %% Loop over entities (either beams or joints)
    itypeGLO = MESH2D.PROP(iENTITY).INDEX_ENTITY ;  % INDEX of the 3D mesh associated to "iENTITY"
    switch MESH2D.PROP(iENTITY).TYPE ;  % Type of entity
        case {'RVE'}
            % It is an RVE
            for iLOC = 1:length(itypeGLO)  % Loop over types of RVES (mixed beams)
                iMESH = iMESH +1 ;  % Index of mesh to print
                itypeSLICE = itypeGLO(iLOC) ; % Type of slice
                [MESHcluster,DATAmeshR] =MeshCluster3D_RVE(MESH2D,MESH3D,iENTITY,itypeSLICE,DATAIN,DATA_REFMESH) ;  % Gathering COORDINATES and CONNECTIVITIES
               
                COORprint = [COORprint; MESHcluster.COOR] ; % Coordinates to print
                CNprint{iMESH} = MESHcluster.CN  + iacumNODES;  % Connectivities to print. We sum the number of accumulated nodes
                TypeElementPRINT{iMESH} =  MESH3D.RVES(itypeSLICE).DATA3D.TypeElement ;
                MaterialTypegloPRINT{iMESH} = MESHcluster.MaterialType + imatACUM ;
                
                
                NAME_LOC_MSH = [MESH3D.RVES(itypeSLICE).NAME] ;
                [~,NAME_LOC_MSH] = fileparts(NAME_LOC_MSH) ;
                NAMEMESH{iMESH} = NAME_LOC_MSH ;
                if ~isempty(MESHcluster.CNlateral)
                    iMESH_LAT = iMESH_LAT +1 ;
                    CNprint_LAT{iMESH_LAT} = MESHcluster.CNlateral  + iacumNODES;
                    if ~isempty(MESH3D.RVES(itypeSLICE).DATA3D.TypeElementB_coarse )
                        TypeElementPRINT_LAT{iMESH_LAT} =  MESH3D.RVES(itypeSLICE).DATA3D.TypeElementB_coarse ;
                    else
                        TypeElementPRINT_LAT{iMESH_LAT} =  MESH3D.RVES(itypeSLICE).DATA3D.TypeElementB  ;
                    end
                    MaterialTypegloPRINT_LAT{iMESH_LAT} = MESHcluster.MaterialTypeglo_lat + imatACUM ;
                    NAMEMESH_LAT{iMESH_LAT} = [NAME_LOC_MSH,'_lat'] ;
                    imatACUM = imatACUM + length(unique(MESHcluster.MaterialType )) +1;
                    
                else
                    imatACUM = imatACUM + length(unique(MESHcluster.MaterialType )) ; % Updating number of materials
                    
                end
                iacumNODES = iacumNODES + size(MESHcluster.COOR,1) ;  % Updating number of nodes
            end
%         case 'JOINT'
%             
%             iMESH = iMESH +1 ;  % Index of mesh to print
%             itypeJOINT = itypeGLO ; % Type of joint
%             [MESHcluster] =MeshCluster3Djoints(MESH2D,MESH3D,iENTITY,itypeJOINT,[]) ;  % Gathering COORDINATES and CONNECTIVITIES
%             COORprint = [COORprint; MESHcluster.COOR] ; % Coordinates to print
%             CNprint{iMESH} = MESHcluster.CN  + iacumNODES;  % Connectivities to print
%             TypeElementPRINT{iMESH} =  MESH3D.JOINTS(itypeSLICE).DATA3D.TypeElement ;
%             MaterialTypegloPRINT{iMESH} = MESHcluster.MaterialType + imatACUM ;
%             
%             NAME_LOC_MSH = [MESH3D.JOINTS(itypeJOINT).DATA3D.NAME_SLICE] ;
%             [~,NAME_LOC_MSH] = fileparts(NAME_LOC_MSH) ;
%             NAMEMESH{iMESH} = NAME_LOC_MSH ;
%             
%             iacumNODES = iacumNODES + size(MESHcluster.COOR,1) ;
%             imatACUM = imatACUM + length(unique(MESHcluster.MaterialType )) ;
        otherwise
            error('Option not implemented')
            
    end
    
end

DATAIN = DefaultField(DATAIN,'PlotUnified3Dmeshes',1) ;

if DATAIN.DISP3D_lateral == 1
    DATAIN.PlotUnified3Dmeshes = 0 ;
end



if   DATAIN.PlotUnified3Dmeshes  == 1
    %     COORprint,CNprint,NAME_INPUT_DATA,...
    %     TypeElementPRINT,MaterialTypegloPRINT,NAMEMESH
    CNprint_glo{1} = CNprint{1} ;
    CNprint_glo{2}= cell2mat(CNprint(2:end)');
    CNprint =  CNprint_glo ;
    
    TypeElementPRINT_glo{1}=TypeElementPRINT{1} ;
    TypeElementPRINT_glo{2}= TypeElementPRINT{2};
    TypeElementPRINT = TypeElementPRINT_glo ;
    
    MaterialTypegloPRINT_glo{1} = MaterialTypegloPRINT{1} ;
    MaterialTypegloPRINT_glo{2}= cell2mat(MaterialTypegloPRINT(2:end)');
    MaterialTypegloPRINT = MaterialTypegloPRINT_glo;
    
    
    NAMEMESH_glo{1} ='1D' ;
    NAMEMESH_glo{2} ='3D' ;
    NAMEMESH = NAMEMESH_glo ;
else
    if  ~isempty(CNprint_LAT)
        
        for imeshLAT = 1:length(CNprint_LAT)
            CNprint{end+1} = CNprint_LAT{imeshLAT} ;
            TypeElementPRINT{end+1} = TypeElementPRINT_LAT{imeshLAT} ;
            MaterialTypegloPRINT{end+1} = MaterialTypegloPRINT_LAT{imeshLAT} ;
            NAMEMESH{end+1} = NAMEMESH_LAT{imeshLAT} ;
            
            
        end
        
    end
end




NAME_INPUT_DATA = [FOLDER,filesep,DATAIN.NAME_INPUT_DATA] ;





%[NameFile_msh,NameFile_res] = NameFilePostProcessGid(DATAIN,NAME_INPUT_DATA) ;
DATAIN = DefaultField(DATAIN,'NameFile_msh',[])  ; 
if isempty(DATAIN.NameFile_msh)
    [NameFile_msh,NameFile_res] = NameFilePostProcessGid(DATAIN,NAME_INPUT_DATA) ;
else
    NameFile_msh = DATAIN.NameFile_msh ; 
    NameFile_res = DATAIN.NameFile_res ; 
    
end



%