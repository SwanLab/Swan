function MESH = GeometryMesh(DATA)

if nargin == 0
    load('tmp1.mat')
end

% See DOCS, KW:GeometryMesh
% JAHO, 22-Nov-2020. Geometric parameters from the FE mesh DATA.NameFileMeshDATA

% JAHO, 27-JAN-2023. CHANGE IN GID'S VERSION (15 ONWARDS). ELEMENTS OF DIFFERENT SURFACES/VOLUMES
% WRITTEN IN DIFFERENT MODULES

DATA =DefaultField(DATA,'GID_VERSION_WITH_MESHES_DIFFERENT_SURFACES_SEPARATED',1) ; 
DATA = DefaultField(DATA,'READ_MATERIAL_COLUMN',1) ; 

if DATA.READ_MATERIAL_COLUMN == 1
    if DATA.GID_VERSION_WITH_MESHES_DIFFERENT_SURFACES_SEPARATED == 0
        [MESH]= ReadMeshFileStr(DATA.NameFileMeshDATA,'READ_MATERIAL_COLUMN',1)  ;
    else
        [MESH]= ReadMeshFileStr_MULT(DATA.NameFileMeshDATA,'READ_MATERIAL_COLUMN',1)  ;
    end
    
else
    % 23-March-2024
     [MESH]= ReadMeshFileStr(DATA.NameFileMeshDATA,'READ_MATERIAL_COLUMN',0)  ;
end

DATA = DefaultField(DATA,'MERGE_ELEMENTS_MESHES_UNCOUPLED_EIFEM',0) ;
if DATA.MERGE_ELEMENTS_MESHES_UNCOUPLED_EIFEM == 1
    % 23-Apr-2024, 
    % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/06_UNCOUPLED.mlx
     MESH = MergeMeshUncoupledEIFEM(MESH) ; 
end


DATA = DefaultField(DATA,'CHANGE_COORDINATES_MESH',[]) ; 

if ~isempty(DATA.CHANGE_COORDINATES_MESH)
    Xref = DATA.CHANGE_COORDINATES_MESH.REFERENCE_POINT ; 
    R    = DATA.CHANGE_COORDINATES_MESH.ROTATION_MATRIX  ; 
    T = DATA.CHANGE_COORDINATES_MESH.TRANSLATION  ;   
    
    % ROTATION 
    % --------
    X = MESH.COOR ; 
    for idim=1:size(X,2)
        X(:,idim) = MESH.COOR(:,idim)-Xref(idim) ;         
    end
    X = (R*X')'; 
    % TRANSLATION 
  
    for idim=1:size(X,2)
        MESH.COOR(:,idim) =X(:,idim) +  T(idim) ;         
    end 
    
end



nCOOR = size(MESH.COOR,1) ;
nCN = length(unique(MESH.CN));
if nCOOR ~= nCN
 
           error('Delete extra-points in the MESH. Some nodes are not used in the connectivity matrix (or you forgot to assign materials)')
 
end

% READING .dat files (faces and lines). Generated with GID's problemtype PROBLEMTYPES_GID/SLICE_JOINTS.gid
[MESH.NODES_FACES,MESH.NODES_LINES,MESH.NODES_DOM] = NodesFacesLinesGID(DATA.NameFileMeshDATA(1:end-4)) ;
%%%%%%
ndim = size(MESH.COOR,2);   %CH2D
if ndim==2
    MESH.NODES_DOM = MESH.NODES_FACES ; 
    MESH.NODES_FACES = MESH.NODES_LINES ;
 
end 


DATA = DefaultField(DATA,'Convert27NodedHexahInto26NodedElement',0) ; % = 1 M 

if DATA.Convert27NodedHexahInto26NodedElement == 1 && size(MESH.CN,2) == 27
    MESH = RemoveInteriorPointHexa27(MESH) ; 
end





% Connectivities faces 
 
Indexes_faces_bnd_element = cell(1,length(MESH.NODES_FACES)) ;
NORMALv = cell(1,length(MESH.NODES_FACES)) ; 
TANGENTv = cell(1,length(MESH.NODES_FACES)) ; ;
for iface = 1:length(MESH.NODES_FACES)
    [dummy, setBelemLOC]= ElemBnd(MESH.CNb,MESH.NODES_FACES{iface}); % elements face "iface"
    % Connectivities faces f1 and f2
    CONNECTb_faces_iface = MESH.CNb(setBelemLOC,:) ;
    % The above is the connectivity matrix for the nodes of face "iface"   
    Indexes_faces_bnd_element{iface} = setBelemLOC ; 
    
    
    
   % KW:NORMALSlocal
    [NORMALv{iface},TANGENTv{iface}] = NormalsBoundaryLocal(MESH.COOR,CONNECTb_faces_iface )  ;
    
     
end

%MESH.CONNECTb_faces = CONNECTb_faces; 
MESH.Indexes_faces_bnd_element = Indexes_faces_bnd_element; 
MESH.NormalBoundaryElementsFace = NORMALv ; 
MESH.TangentBoundaryElementsFace = TANGENTv ; 

% Moved inside on 22-March-2021
%***********************************

DATA = DefaultField(DATA,'RenumberElementsForEficiency',1) ;  
if  DATA.RenumberElementsForEficiency == 1
    disp('Renumering elements')
    [~,IndicesRenumberingElements]  = sort(MESH.CN(:,1)) ;
    MESH.CN = MESH.CN(IndicesRenumberingElements,:) ;
  %  MATPRO.celasglo = MATPRO.celasglo(:,:,IndicesRenumberingElements) ;
  %  MATPRO.dens = MATPRO.dens(IndicesRenumberingElements) ;
    MESH.MaterialType = MESH.MaterialType(IndicesRenumberingElements) ;
else
    IndicesRenumberingElements = [] ;
end
MESH.IndicesRenumberingElements =IndicesRenumberingElements ; 

