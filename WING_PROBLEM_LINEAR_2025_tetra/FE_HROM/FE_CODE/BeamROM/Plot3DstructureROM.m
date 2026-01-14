function [NameFile_msh,NameFile_res,COORprint,CNprint,NAME_INPUT_DATA,...
      TypeElementPRINT,MaterialTypegloPRINT,NAMEMESH]= Plot3DstructureROM(MESH1D,MESH3D,DATAIN)
% Plot 3D domains along with 1D skeleton
% JAHO, 22 Jan-2018
if nargin == 0
    load('tmp1.mat')
end

iacumCN =  size(MESH1D.COOR,1) ; % Connectivities index (cumulative, for printing global 3D mesh)
COORprint = [MESH1D.COOR] ;  % Coordinates global GID mesh
CNprint =  {MESH1D.CN};  % Connectivities global GID mesh
TypeElementPRINT=  {'Linear'} ;  % Type of element global GID mesh
MaterialTypegloPRINT = {ones(size(MESH1D.CN,1),1)} ; % Type of element global GID mesh
NAMEMESH =   {'1D'} ;
itypeACUM = 1;
imatACUM = 1 ;




for iENTITY  = 1: length(MESH1D.PROP) %% Loop over entities (either beams or joints)
    itypeGLO = MESH1D.PROP(iENTITY).INDEX_ENTITY ;  % INDEX of the 3D mesh associated to "iENTITY"
    switch MESH1D.PROP(iENTITY).TYPE ;  % Type of entity
        case {'BEAM'}
            % It is a beam
            for iLOC = 1:length(itypeGLO)  % Loop over types of slices (mixed beams)
                itypeACUM = itypeACUM +1 ;  % Index of mesh to print
                itypeSLICE = itypeGLO(iLOC) ; % Type of slice
                [MESHcluster] =MeshCluster3D(MESH1D,MESH3D,iENTITY,itypeSLICE,[]) ;  % Gathering COORDINATES and CONNECTIVITIES
                COORprint = [COORprint; MESHcluster.COOR] ; % Coordinates to print
                CNprint{itypeACUM} = MESHcluster.CN  + iacumCN;  % Connectivities to print
                TypeElementPRINT{itypeACUM} =  MESH3D.SLICES(itypeSLICE).DATA3D.TypeElement ;
                MaterialTypegloPRINT{itypeACUM} = MESHcluster.MaterialType + imatACUM ;
                NAMEMESH{itypeACUM} = [MESH3D.SLICES(itypeSLICE).DATA3D.NAME_SLICE] ;
                iacumCN = iacumCN + size(MESHcluster.COOR,1) ;
                imatACUM = imatACUM + length(unique(MESHcluster.MaterialType )) ;
            end
        case 'JOINT'
            
            if all(MESH1D.SUCCESSFUL_MATCH{iENTITY})
                itypeACUM = itypeACUM +1 ;  % Index of mesh to print
                itypeJOINT = itypeGLO ; % Type of slice
                [MESHcluster] =MeshCluster3Djoints(MESH1D,MESH3D,iENTITY,itypeJOINT,[]) ;  % Gathering COORDINATES and CONNECTIVITIES
                COORprint = [COORprint; MESHcluster.COOR] ; % Coordinates to print
                CNprint{itypeACUM} = MESHcluster.CN  + iacumCN;  % Connectivities to print
                TypeElementPRINT{itypeACUM} =  MESH3D.SLICES(itypeSLICE).DATA3D.TypeElement ;
                MaterialTypegloPRINT{itypeACUM} = MESHcluster.MaterialType + imatACUM ;
                NAMEMESH{itypeACUM} = [MESH3D.JOINTS(itypeJOINT).DATA3D.NAME_SLICE] ;
                iacumCN = iacumCN + size(MESHcluster.COOR,1) ;
                imatACUM = imatACUM + length(unique(MESHcluster.MaterialType )) ;
            end
            
    end
    
end




NAME_INPUT_DATA = [] ;
%GID PRINTING
% [dummy1 NAMEMESH_STRUCTURE    dummy2]= fileparts(MESH1D.NAME) ;
% NameFile_msh = ['GIDPOST',filesep,NAMEMESH_STRUCTURE,'_',DATAIN.NAME_project,'.msh'] ;


[NameFile_msh,NameFile_res] = NameFilePostProcessGid(DATAIN,DATAIN.NAME_INPUT_DATA) ;

% 
% IND_ELEM_MESHES = GidMesh2DFE_multi(NameFile_msh,COORprint,CNprint,NAME_INPUT_DATA,MaterialTypegloPRINT,TypeElementPRINT,NAMEMESH);
% 
% 
% 
% 
% disp(['Open GID file  : ']) ;
% disp([NameFile_msh]) ;