function Plot3DstructureRVE(MESH2D,MESH3D,DATAIN)
% Plot 3D domains along with 2D skeleton
% JAHO, 22 Jan-2018/20-July-2018
if nargin == 0
    load('tmp.mat')
end

iacumCN =  size(MESH2D.COORvert,1) ; % Connectivities index (cumulative, for printing global 3D mesh)
COORprint = [MESH2D.COORvert] ;  % Coordinates global GID mesh
CNprint =  {MESH2D.CNvert};  % Connectivities global GID mesh
TypeElementPRINT=  {MESH2D.TypeElement} ;  % Type of element global GID mesh
MaterialTypegloPRINT = {ones(size(MESH2D.CN,1),1)} ; % Type of element global GID mesh
NAMEMESH =   {'2D'} ;
itypeACUM = 1;
imatACUM = 1 ;
for iENTITY  = 1: length(MESH2D.PROP) %% Loop over entities (either beams or joints)
    itypeGLO = MESH2D.PROP(iENTITY).INDEX_ENTITY ;  % INDEX of the 3D mesh associated to "iENTITY"
    switch MESH2D.PROP(iENTITY).TYPE ;  % Type of entity
        case {'RVE'}
            % It is a beam
 %           for iLOC = 1:length(itypeGLO)  % Loop over types of slices (mixed beams)
                itypeACUM = itypeACUM +1 ;  % Index of mesh to print
                itypeRVE = itypeGLO ; %itypeGLO(iLOC) ; % Type of slice
                [MESHcluster] =MeshCluster3D_RVE(MESH2D,MESH3D,iENTITY,itypeRVE,[]) ;  % Gathering COORDINATES and CONNECTIVITIES
                COORprint = [COORprint; MESHcluster.COOR] ; % Coordinates to print
                CNprint{itypeACUM} = MESHcluster.CN  + iacumCN;  % Connectivities to print
                TypeElementPRINT{itypeACUM} =  MESH3D.RVES(itypeRVE).DATA3D.TypeElement ;
                MaterialTypegloPRINT{itypeACUM} = MESHcluster.MaterialType + imatACUM ;
                NAMEMESH{itypeACUM} = [MESH3D.RVES(itypeRVE).NAME] ;
                iacumCN = iacumCN + size(MESHcluster.COOR,1) ;
                imatACUM = imatACUM + length(unique(MESHcluster.MaterialType )) ;
  %          end
        otherwise
            error('Option not implemented')
            
    end
    
end




NAME_INPUT_DATA = [] ;
%GID PRINTING
[dummy1 NAMEMESH_STRUCTURE    dummy2]= fileparts(MESH2D.NAME) ;
NameFile_msh = ['GIDPOST',filesep,NAMEMESH_STRUCTURE,'_',DATAIN.NAME_project,'.msh'] ;

IND_ELEM_MESHES = GidMesh2DFE_multi(NameFile_msh,COORprint,CNprint,NAME_INPUT_DATA,MaterialTypegloPRINT,TypeElementPRINT,NAMEMESH);



disp(['Open GID file  : ']) ;
disp([ cd,filesep,NameFile_msh]) ;