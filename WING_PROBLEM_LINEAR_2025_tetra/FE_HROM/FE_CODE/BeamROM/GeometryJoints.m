function DATA3D = GeometryJoints(MESH1D,MESH3D,itypeJOINT,DATAIN)
% Obtain: Coordinates, Connectivities and other geometric entities
% such as face f1 and f2 for domain decomposition purposes
% of GID finite element mesh NAMEMESH.
% See BeamROM.pdf
% JAHO-5-FEb-2018
%

if nargin == 0
    load('tmp0.mat')
end

JOINT = MESH3D.JOINTS(itypeJOINT) ; % Name of the file containing the mesh of the JOINT

% READING .msh file
% ------------------
READ_MATERIAL_COLUMN = 1;
[DATA3D]= ReadMeshFileStr(JOINT.NAME,'READ_MATERIAL_COLUMN',READ_MATERIAL_COLUMN)  ;

% READING .dat files (faces and lines). Generated with GID's problemtype PROBLEMTYPES_GID/SLICE_JOINTS.gid
[DATA3D.NODES_FACES,DATA3D.NODES_LINES] = NodesFacesLinesGID(JOINT.NAME) ;
%%%%%%
[dummy1 DATA3D.NAME_SLICE    dummy2]= fileparts(JOINT.NAME) ;
% ------------------------------------------------------
%%% COMPUTING NORMALS ALL BOUNDARY ELEEMENTS
% ------------------------------------------
%NORMALSv = NormalsBoundary(DATA3D.COOR,DATA3D.CNb)  ;
% ----------------------------------------------
% Identificataion nodes faces corresponding to the neighboring slices
% -------------------------------------------------------------------

% Centroids NODES_FACES (JOINT)  -Pseudo-centroids 
CENTROIDS = cell(size(DATA3D.NODES_FACES)) ;
for iface = 1:length(DATA3D.NODES_FACES)
    CENTROIDS{iface} = sum(DATA3D.COOR(DATA3D.NODES_FACES{iface},:),1)'/length(DATA3D.NODES_FACES{iface}) ;
end
CENTROIDS = cell2mat(CENTROIDS) ;


% NODES FACES (SLICES)
WRITE_BATCH_FILE = 0 ;
CONNECTslices = CreateNewJoint(MESH1D,MESH3D,itypeJOINT,WRITE_BATCH_FILE) ;

% Identification (pairing nodes)
% -------------------------------
DATAIN = DefaultField(DATAIN,'ToleranceIdentificationJointSliceAbsolute',[]) ; 
if isempty(DATAIN.ToleranceIdentificationJointSliceAbsolute)
TOL = 1e-1 ; % Local tolerance
CNb1 = DATA3D.CNb(1,1:2) ;
COOR1  = DATA3D.COOR(CNb1,:) ;
dCOOR = norm(COOR1(1,:)-COOR1(2,:)) ;
TOL = TOL*dCOOR ; % Absolute tolerance
else
   TOL =  DATAIN.ToleranceIdentificationJointSliceAbsolute ; 
end
%
ROTATION_SLICES = cell(1,length(CONNECTslices.NODES) ) ;
% The reference points are the CENTROIDS --- (face1,face2, face3....)
% determined  by the user via GID's interface
for iface =1:length(CONNECTslices.NODES) % Loop over nodes determined by the 1D mesh
    COORslice = CONNECTslices.COOR{iface} ; % Coordinate of the iface-th contacting face
    
    if ~isempty(COORslice)
        xG = sum(COORslice,2)/size(COORslice,2) ;  % Centroid (approximated). This should be modified !!!!
        
   
        
        
        % Identification
        ndim= 3;
        DIST  = zeros(size(CENTROIDS)) ;
        for idim = 1:size(DIST,1)
            DIST(idim,:) = CENTROIDS(idim,:)-xG(idim) ;
        end
        nDIST = sqrt(sum(DIST.^2,1)) ;
        ifaceJOINT = find(nDIST<=TOL)  ;
        if isempty(ifaceJOINT)
            disp(['Minimum distance = ',num2str(min(nDIST))])
            disp(['This distance is greater than the prescribed tolerance =',num2str(TOL)])
            disp(['Reset variable DATAIN.ToleranceIdentificationJointSliceAbsolute'])
            disp(['to value: ',num2str(min(nDIST))])
            error('Error in identifying connecting faces of joints (no coincidence of centroids). Perhaps you have to provide a skeleton file for the joint (JOINT.MESH1D_SKELETON.STRUCTURE.NAME) ---or slice faces are not correctly defined in GID')
        end
        ROTATION_SLICES{ifaceJOINT} = CONNECTslices.ROTATION{iface} ;
        COOR_joint = DATA3D.COOR(DATA3D.NODES_FACES{ifaceJOINT},:) ;
        
        if size(COOR_joint,1) ~= size(COORslice',1)
            error('Non-conforming slice/joint')
        end
        
        [IDX DISTANCES]= knnsearch(COOR_joint,COORslice') ;
        
        
        
        if any(find(DISTANCES > TOL))
            error('Non-conforming meshes. Check tolerance TOL')
        end
        DATA3D.NODES_FACES{ifaceJOINT} = DATA3D.NODES_FACES{ifaceJOINT}(IDX) ;
   
    end
    
end


DATA3D.TOLERANCE_DISTANCE_NODES = TOL ;
DATA3D.CENTROIDS = CENTROIDS ;
DATA3D.CONNECTslices = CONNECTslices.indexSLICE_GLO ;   % Slices connected to the joints
%DATA3D.ROTATIONS_CONNECTING_SLICES = ROTATION_SLICES ;



