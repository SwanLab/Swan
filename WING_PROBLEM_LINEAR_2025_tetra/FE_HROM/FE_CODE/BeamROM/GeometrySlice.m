function DATA3D = GeometrySlice(SLICE,DATAIN)
% Obtain: Coordinates, Connectivities and other geometric entities
% such as face f1 and f2 for domain decomposition purposes
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
%          CENTRf1:  ---> Coordinates centroid f1
%          CENTRf2:  ---> Coordinates centroid f2
%     ROTATIONf1f2:  ---> Rotation matrix face 2 with respect face 1, (curved slices)
%
%
%
% JAHO-28-Dic-2017
%

if nargin == 0
    load('tmp.mat')
elseif nargin == 1
    DATAIN = [] ;
end

% READING .msh file
% ------------------
DATAIN = DefaultField(DATAIN,'READ_MATERIAL_COLUMN',1) ;
[DATA3D]= ReadMeshFileStr(SLICE.NAME,'READ_MATERIAL_COLUMN',DATAIN.READ_MATERIAL_COLUMN)  ;

% READING .dat files (faces and lines). Generated with GID's problemtype PROBLEMTYPES_GID/SLICE_JOINTS.gid
[DATA3D.NODES_FACES,DATA3D.NODES_LINES] = NodesFacesLinesGID(SLICE.NAME) ;
ndim = size(DATA3D.COOR,2);
if ndim==2
    DATA3D.NODES_FACES = DATA3D.NODES_LINES ;
end
%%%%%%
[dummy1 DATA3D.NAME_SLICE    dummy2]= fileparts(SLICE.NAME) ;
% ------------------------------------------------------
%%% COMPUTING NORMALS ALL BOUNDARY ELEEMENTS
% ------------------------------------------
NORMALSv = NormalsBoundary(DATA3D.COOR,DATA3D.CNb)  ;
% ----------------------------------------------
% Identificataion nodes faces
[NODES_face1, NODES_face2, CENTRf1, CENTRf2, ROTATIONf1f2,CENTRf1_real,CENTRf2_real] =...
    SliceNodesFacesIdentification(DATA3D,NORMALSv,DATAIN) ;

DATA3D.NODES_FACES{1} = NODES_face1 ;
DATA3D.NODES_FACES{2} = NODES_face2 ;
DATAIN = DefaultField(DATAIN,'USE_CONSISTENT_CENTROIDS',1) ;   % Irrelevant
if DATAIN.USE_CONSISTENT_CENTROIDS == 0
    DATA3D.CENTRf1 = CENTRf1 ;
    DATA3D.CENTRf2 = CENTRf2 ;
else
    DATA3D.CENTRf1 = CENTRf1_real ;
    DATA3D.CENTRf2 = CENTRf2_real ;
end
DATA3D.CENTRf1_real = CENTRf1_real ;
DATA3D.CENTRf2_real = CENTRf2_real ;
DATA3D.ROTATIONf1f2 = ROTATIONf1f2 ;  % Rotation of face f2 with respect face f1
%DATA3D.CenterRotation = Crotation ;



% DATAIN = DefaultField(DATAIN,'GenerateLinesFaces',0) ; % = 1;

% if DATAIN.GenerateLinesFaces == 1
%     % Face f1
%     % -------
%     % Associated boundary elements
%    [CNbF1 setBelem] =  ElemBnd(CNb,DATA3D.NODES_face1) ;
%    % Writting GID's batch
%    BatchGID_FaceNodes(COOR,CNbF1,NameWS,'face1') ;
%        % Face f2
%     % -------
%     % Associated boundary elements
%    [CNbF2 setBelem] =  ElemBnd(CNb,DATA3D.NODES_face2) ;
%    % Writting GID's batch
%    BatchGID_FaceNodes(COOR,CNbF2,NameWS,'face2') ;
% end








%
%
% %% REMAINING FACES. Determining local axes for remaining   faces
% % -----------------------------------------------------------------------------------
%
% % First tangential vector is obtained as n times xaxis
% xaxis = repmat([1,0,0]',1,size(NORMALSv,2)) ;
% t1 = cross(NORMALSv,xaxis) ;
% t1norm = sqrt(sum(t1.^2,1)) ;
% for idim = 1:size(t1,1)
%     t1(idim,:)  = t1(idim,:)./t1norm ;
% end
% % Second tangential vector
% t2 = cross(NORMALSv,t1) ;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % CLASSIFICATION OF NORMALS.
% BoundaryElements = {} ;
% NormalElements = {} ;
% remELEM = 1:size(NORMALSv,2) ;
% %      remELEM(elements_f12) = [] ;
% EXISTelem = 1 ;
% ibound = 1;
% while length(remELEM)>1
%     refELEM = remELEM(1) ;
%     normalREF = NORMALSv(:,refELEM) ; % Normal reference element
%     if sum(abs(abs(normalREF)-[1,0,0]'))<TOL
%         nodDATA3D = CNb(refELEM,1) ;
%         if  COOR(nodDATA3D,1) == xmin
%             if1 = ibound ;
%         else
%             if2 = ibound ;
%         end
%     end
%     % Find boundary elements with similar normals
%     DIFFNORMAL = bsxfun(@minus,NORMALSv(:,remELEM),normalREF) ;
%     DIFFNORMAL = (sum(abs(DIFFNORMAL),1)) ;
%     elemsSIMnorm = find(DIFFNORMAL<TOL) ;
%     BoundaryElements{ibound} = remELEM(elemsSIMnorm) ;
%     NormalElements{ibound} = normalREF ;
%     remELEM(elemsSIMnorm)  = [] ;
%     ibound = ibound + 1;
% end
%
% BoundaryElementsNEW = cell(length(BoundaryElements),1) ;
% BoundaryElementsNEW{1} = unique(BoundaryElements{if1}) ;
% BoundaryElementsNEW{2} = unique(BoundaryElements{if2}) ;
% iacum = 3;
% for iface = 1:length(BoundaryElementsNEW)
%     if iface ~=if1 & iface ~=if2
%         BoundaryElementsNEW{iacum} =  unique(BoundaryElements{iface})  ;
%         iacum = iacum + 1;
%     end
%
% end
%
% posgp = [] ;
% DATA.NameWS = NameWS;
% GidPostProcess_NORMALS(COOR,CNb,TypeElementB,NORMALSv,t1,t2,posgp,BoundaryElementsNEW,DATA);
%
%
%
%
%
%
%
%
%
