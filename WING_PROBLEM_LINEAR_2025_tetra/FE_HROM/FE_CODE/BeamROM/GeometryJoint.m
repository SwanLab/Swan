function DATA3D = GeometryJoint(JOINT,SLICES,DATAIN)
% See GeometrySlice
% Determining mesh properties of JOINT structure array
% See BeamROM.pdf
% JAHO-28-Dic-2017
%

if nargin == 0
    load('tmp.mat')
end
READ_MATERIAL_COLUMN = 1;
[DATA3D]= ReadMeshFileStr(JOINT.NAME,'READ_MATERIAL_COLUMN',READ_MATERIAL_COLUMN)  ;
if isempty(DATA3D.CNb)
    error('You forgot to generate the boundary mesh (or to assign materials )  !!! ')
end
[dummy1 NameWS    dummy2]= fileparts(JOINT.NAME) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% How to identify transmission faces  ?  This information is to be partially provided
% by variable JOINT.SKELETON_1D
JOINT = DefaultField(JOINT,'SKELETON_1D',[]) ;
if isempty(JOINT.SKELETON_1D)
    error('Variable JOINT.SKELETON_1D is empty. You must create the corresponding 1D mesh.')
end
if exist(JOINT.SKELETON_1D,'file') == 0
    error('Variable JOINT.SKELETON_1D does not exist. You must create the corresponding 1D mesh.')
end
% Reading 1D skeleton
% ------------------------
ndim  = 3 ;
[MESH1D]=ReadMeshFileStr(JOINT.SKELETON_1D,'READ_MATERIAL_COLUMN',1)  ;
if size(MESH1D.COOR,2) < ndim
    MESH1D.COOR = [MESH1D.COOR, zeros(size(MESH1D.COOR,1),1)] ;
end
% Detect 1D boundary nodes of 1D skeleton:  OUT --> BND.nodes, and BND.slices
[BND ] = DetectBoundaryNodesJoints(MESH1D.CN,MESH1D.MaterialType,JOINT.INDEX) ;
nfaces = length(BND.nodes) ; % Number of faces
% To identify the nodes of each tranmission face, we invoke function
% CoordinatesFaceNodes, that returns the coordinates of the slice in
% contact with the joint for the given 1D node
TEXT = {} ;
DATA3D.FACES = [];
for inode = 1:length(BND.nodes)
    node = BND.nodes(inode) ;
    DATA3D.FACES(inode).node1D_local = node ;
    slice = BND.slices(inode) ;
    [COORslice  CNbF1 COORface NODES_face ]= CoordinatesFaceNodes(MESH1D,node,slice,SLICES) ;
    distance_points = norm(COORslice(:,CNbF1(1,1))-COORslice(:,CNbF1(1,2))) ;  % Distance between nodes
    TOL = distance_points*1e-3 ;
    % Search neighboring nodes in MESH3D.CNb
    %     if
    %     % [TEXT] = BatchGID_FaceNodes_points(COORslice',CNbF1,TEXT) ;
    %     end
    BndNodes = unique(DATA3D.CNb(:)) ;
    [IDX DISTANCES]= knnsearch(DATA3D.COOR(BndNodes,:),COORface') ;
    if any(find(DISTANCES > TOL))
        error('Non-conforming meshes. Check tolerance TOL')
    end
    DATA3D.FACES(inode).NODES_FACES_3D = BndNodes(IDX) ;
    CENTROID= sum(DATA3D.COOR(BndNodes(IDX),:)',2)/length(IDX) ;  % Centroid face2
    DATA3D.FACES(inode).CENTROID = CENTROID ;
end
% COMPUTE_NORMALS = 0 ;
% if COMPUTE_NORMALS==1
%     NORMALSv = zeros(size(COOR,2),size(CNb,1));
%     Xloc1 = COOR(CNb(:,1),:)';
%     Xloc2 = COOR(CNb(:,2),:)';
%     Xloc3 = COOR(CNb(:,3),:)';
%     t1 = Xloc2-Xloc1 ; % Tangential vector
%     t2 = Xloc3-Xloc1 ; % Tangential vector
%     n = cross(t1,t2) ;
%     normN = sqrt(sum(n.^2,1)) ;
%     for idim = 1:size(n,1)
%         NORMALSv(idim,:)  = n(idim,:)./normN ;
%     end
%     
% end
