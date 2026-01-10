function DATA3D = GeometryRVE(SLICE,DATAIN)
% Obtain: Coordinates, Connectivities and other geometric entities
% such as face f1 and f2 ..... for domain decomposition purposes
% of GID finite element mesh NAMEMESH.
% See BeamROM.pdf
%
% %       COOR:    --> Coordinates
%               CN:  ---> Connectivities
%      TypeElement:  ---> (Hexahedra, triangle...)
%              CNb:   ---> Connectivitiy boundary elements
%     TypeElementB:   ---> Type of element boundary
%     MaterialType:   ---> Material array
%      NODES_FACES:   ---> NODES_FACES{i}, nodes pertaining to face i (labelled with GID)
%      NODES_LINES:  ---> NODES_LINES{i}, nodes pertaining to line i (labelled with GID)
%       NAME_SLICE:

% JAHO-28-Dic-2017/20-July-2018/1-Sept-2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    load('tmp2.mat')
elseif nargin == 1
    DATAIN = [] ;
end

% READING .msh file
% ------------------
DATAIN = DefaultField(DATAIN,'READ_MATERIAL_COLUMN',1) ;

%[DATA3D]= ReadMeshFileStr(SLICE.NAME,'READ_MATERIAL_COLUMN',DATAIN.READ_MATERIAL_COLUMN)  ;
[DATA3D]= ReadMeshFileStr_MULT(SLICE.NAME,'READ_MATERIAL_COLUMN',DATAIN.READ_MATERIAL_COLUMN)  ;  % JAHO, 29-Apr-2024

% %CNrenumbered = RenumberConnectivities(CN,NODES)
nCOOR = size(DATA3D.COOR,1) ;
nCN = length(unique(DATA3D.CN));
if nCOOR ~= nCN
    error('Delete extra-points in the MESH. Some nodes are not used in the connectivity matrix')
end

% READING .dat files (faces and lines). Generated with GID's problemtype PROBLEMTYPES_GID/SLICE_JOINTS.gid
[DATA3D.NODES_FACES,DATA3D.NODES_LINES] = NodesFacesLinesGID(SLICE.NAME) ;
%%%%%%
[dummy1 DATA3D.NAME_SLICE    dummy2]= fileparts(SLICE.NAME) ;
% ------------------------------------------------------
ndim = size(DATA3D.COOR,2);   %CH2D
if ndim==2
    DATA3D.NODES_FACES = DATA3D.NODES_LINES ;
end

DATAIN = DefaultField(DATAIN,'ContinuumStructuresWithoutCorners',1) ; %  = 1;
% IMP-MAY19-A
  

%%% COMPUTING NORMALS ALL BOUNDARY ELEEMENTS
% ------------------------------------------
NORMALSv = NormalsBoundary(DATA3D.COOR,DATA3D.CNb)  ;
% ----------------------------------------------
% Pairing nodes faces
DATA3D = SliceNodesFacesIdentificationRVE2D(DATA3D,NORMALSv,DATAIN) ;

%% For plotting purposes, we need the coordinates of the midside nodes

% Coordinate of midside of face 1
xCORNER = DATA3D.CENTRf{1}(1) ;
yCORNER = DATA3D.CENTRf{2}(2) ;

% CH2D -
% -----
if length(DATA3D.CENTRf{1}) == 3
    zCORNER = DATA3D.CENTRf{1}(3) ;
    BOTTOM_LEFT_COOR = [xCORNER,yCORNER,zCORNER]';
else
    BOTTOM_LEFT_COOR = [xCORNER,yCORNER]';
end



Lx = DATA3D.CENTRf{3}(1) - DATA3D.CENTRf{1}(1) ;
Ly = DATA3D.CENTRf{4}(2) - DATA3D.CENTRf{2}(2) ;
% CH2D -
if length(DATA3D.CENTRf{1}) == 3
    
    CORNER_1_CENTROID = zeros(3,1) ;
    CORNER_1_CENTROID(1) =  DATA3D.CENTRf{1}(1) ;
    CORNER_1_CENTROID(2) =  DATA3D.CENTRf{4}(2) ;
    CORNER_1_CENTROID(3) =  DATA3D.CENTRf{4}(3) ;
    DATA3D.COOR_MIDSIDE_FACE{1} = BOTTOM_LEFT_COOR +[0,Ly/2,0]' ;
    DATA3D.COOR_MIDSIDE_FACE{2} = BOTTOM_LEFT_COOR +[Lx/2,0,0]' ;
    DATA3D.COOR_MIDSIDE_FACE{3} = BOTTOM_LEFT_COOR +[Lx,Ly/2,0]' ;
    DATA3D.COOR_MIDSIDE_FACE{4} = BOTTOM_LEFT_COOR +[Lx/2,Ly,0]' ;
    
    DATA3D.COOR_CORNER{1} = BOTTOM_LEFT_COOR + [0,Ly,0]' ;
    DATA3D.COOR_CORNER{2} = BOTTOM_LEFT_COOR  ;
    DATA3D.COOR_CORNER{3} = BOTTOM_LEFT_COOR + [Lx,0,0]' ;
    DATA3D.COOR_CORNER{4} = BOTTOM_LEFT_COOR + [Lx,Ly,0]' ;
else
    DATA3D.COOR_MIDSIDE_FACE{1} = BOTTOM_LEFT_COOR +[0,Ly/2]' ;
    DATA3D.COOR_MIDSIDE_FACE{2} = BOTTOM_LEFT_COOR +[Lx/2,0]' ;
    
    DATA3D.COOR_MIDSIDE_FACE{3} = BOTTOM_LEFT_COOR +[Lx,Ly/2]' ;
    DATA3D.COOR_MIDSIDE_FACE{4} = BOTTOM_LEFT_COOR +[Lx/2,Ly]' ;
    CORNER_1_CENTROID = zeros(2,1) ;
    CORNER_1_CENTROID(1) =  DATA3D.CENTRf{1}(1) ;
    CORNER_1_CENTROID(2) =  DATA3D.CENTRf{4}(2) ;
    
    DATA3D.COOR_CORNER{1} = BOTTOM_LEFT_COOR + [0,Ly]' ;
    DATA3D.COOR_CORNER{2} = BOTTOM_LEFT_COOR  ;
    DATA3D.COOR_CORNER{3} = BOTTOM_LEFT_COOR + [Lx,0]' ;
    DATA3D.COOR_CORNER{4} = BOTTOM_LEFT_COOR + [Lx,Ly]' ;
    
end


DATA3D.CORNER_1_CENTROID = CORNER_1_CENTROID ;

%%%% Nodes from corners and sides (for continuous domains)
% ------------------------------------------------------------
[NODES_CORNERS,NODES_SIDES,NODES_FACES] = CornerSideNodes(DATA3D) ;
DATA3D.NODES_CORNERS = NODES_CORNERS ;
DATA3D.NODES_SIDES = NODES_SIDES ;



DATA3D.NODES_FACES = NODES_FACES ;  % Re-sorted inside so that
% Face1 = Corner1-Corner2-Side1 ....



end

% Length
