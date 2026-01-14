function MESHnew = MergeMeshUncoupledEIFEM_2D_1mat(MESH,INDT,DATAcommon)
if nargin == 0
    load('tmp2.mat')
end
%--------------------------------------------------------------------------------------
% INPUT
% %------------------------------------------------------------------------------------
% EIFEM, TREATMENT OF UNIT CELLS WITH INDEPENDENT DISPLACEMENTS
%
MESHnew = MESH ;
disp('Uncoupled meshes, EIFEM, case 2D')
indMAT = cell2mat(INDT) ;
nMESH = length(indMAT);
disp(['Material index in GID identifying   elements only coupled in the 1st direction) = ',num2str(indMAT(1))]) ;
disp(['Material index in GID identifying   elements only coupled in the 2nd direction) = ',num2str(indMAT(2))]) ;

ELEMS = cell(1,nMESH) ;
for imat = 1:length(indMAT)
    ELEMS{imat} = find(MESH.MaterialType == indMAT(imat)) ;
end

% PAIRING ELEMENTS
% CENTROIDS
if length(ELEMS{1}) ~= length(ELEMS{2})
    error('The number of elements of both meshes should be the same ')
end

CENTROIDS  =cell(1,nMESH) ;
ndim  =size(MESH.COOR,2)  ;
for imat = 1:length(indMAT)
    CNloc = MESH.CN(ELEMS{imat},:) ;
    COOR_cent  = zeros(size(CNloc,1),ndim)  ;
    for inodeLOC = 1:size(CNloc,2)
        COOR_cent =  COOR_cent + MESH.COOR(CNloc(:,inodeLOC),:) ;
    end
    CENTROIDS{imat} = COOR_cent/size(CNloc,2) ;
end

[Ind1,DIST] = knnsearch(CENTROIDS{1},CENTROIDS{2}) ;

if norm(DIST) >1e-10
    error('Centroids do not coincide')
end

% NEW CONNECTIVITY MATRIX
imat = 1;
ELEMS{imat} =   ELEMS{imat}(Ind1) ;
MESHnew.CN = [MESH.CN(ELEMS{1},:),MESH.CN(ELEMS{2},:)] ;
% New material type
MESHnew.MaterialType = MESH.MaterialType(ELEMS{1}) ;

% SUPPORT MESH (FOR POST-PROCESSING PURPOSES)
% Find points with common location
% NODE_OLD_1     NODE_OLD_1 NODE_OLD_2 NODE_OLD_3 ....NODE_OLD_

[III,DDDD] = knnsearch(MESH.COOR,MESH.COOR) ;


[IndicesNEW,BBB,Indices_OLD] = unique(III) ;
% Indices_OLD =  nnodes_OLD  x  1
% Indices_OLD(inodeUNC) = inodeSUPPORT --> inodeUNC = Index node uncoupled
% mesh,  inodeSUPPORT --- INDEX NODE SUPPORT MESH
MESHnew.COOR_SUPPORT = MESH.COOR(IndicesNEW,:) ;
nnodeE = size(MESH.CN,2) ;
MESHnew.CN_SUPPORT = Indices_OLD(MESHnew.CN) ;
MESHnew.CN_SUPPORT  = reshape(MESHnew.CN_SUPPORT,size(MESHnew.CN,1),size(MESHnew.CN,2)) ; 
MESHnew.CN_SUPPORT =  MESHnew.CN_SUPPORT(:,1:nnodeE) ;
MESHnew.IndicesNODES_uncoupledMESH_TO_supportMESH = Indices_OLD ;
%
%     OPERATOR SUCH THAT   d_standard = P*d_uncoupled
nnode_UNCOUP = size(MESHnew.COOR,1) ; 
MESHnew.SMOOTH_from_UNCOUP_TO_SUPPORT = Smooth_Uncoup_to_Support_EIFEM(MESHnew.COOR_SUPPORT,IndicesNEW,nnode_UNCOUP) ; 

% SMOOTH_from_UNCOUP_TO_SUPPORT = Smooth_Uncoup_to_Support_EIFEM_fast(MESHnew.COOR_SUPPORT,IndicesNEW,nnode_UNCOUP)

PLOT_GID = 0;

if PLOT_GID == 1
    
    NameFile = 'prueba.msh' ;
    NAMEPROJ = 'dummy';
    
    GidMesh2DFE(NameFile,MESHnew.COOR_SUPPORT,MESHnew.CN_SUPPORT,NAMEPROJ,MESHnew.MaterialType,MESHnew.TypeElement)
    
    
    
end



