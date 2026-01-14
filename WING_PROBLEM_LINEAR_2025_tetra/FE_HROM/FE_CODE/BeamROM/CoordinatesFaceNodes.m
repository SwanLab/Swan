function  [COORtransf CNbF1 COORface NODES_face ]= CoordinatesFaceNodes(MESH1D,node1D,elem1D,SLICES)
% Global coordinates of the face corresponding to 1D node1D "node1D". To
% determine it, we use the mesh information of the elem1D "elem1D" in contact
% with the domain 
% JAHO, 11-January-2017
% ----------------------

coorCENT = MESH1D.COOR(node1D,:);
% Type elem1D
typeSlice = MESH1D.MaterialType(elem1D) ;
% Initial and final node1Ds 1D
NODES = MESH1D.CN(elem1D,:) ;
xINI = MESH1D.COOR(NODES(1),:); xFIN = MESH1D.COOR(NODES(2),:) ;


% Rotation matrix
AngleRotationXaxis = SLICES(typeSlice).ROTATION_AROUND_LOCALaxis ;
R = RotationMatrix_gen(xINI,xFIN,AngleRotationXaxis)  ;
% Determining which face is node1D "node1D", as well as the associated 3D
% node1Ds, and the centroid in local coordinates
if node1D == NODES(1) ;        %
    NODES_face = SLICES(typeSlice).DATA3D.NODES_face1 ;
    CENTR = SLICES(typeSlice).DATA3D.CENTRf1 ;
    CENTR_opp = SLICES(typeSlice).DATA3D.CENTRf2 ;

else
    NODES_face = SLICES(typeSlice).DATA3D.NODES_face2 ;
    CENTR = SLICES(typeSlice).DATA3D.CENTRf2 ;
    CENTR_opp = SLICES(typeSlice).DATA3D.CENTRf1 ;

end
% Coordinates of 3D node1Ds
COOR3D = SLICES(typeSlice).DATA3D.COOR ;
CNb3D = SLICES(typeSlice).DATA3D.CNb ;
% Associated boundary elements
[CNbF1 setBelem] =  ElemBnd(CNb3D,NODES_face) ;
% Writting GID's batch
% Transformation of coordinates
%
TranslationVector = coorCENT-CENTR_opp' ;  % Translation vector (for GID's batch file)
COORtransf = CoordinatesChange(COOR3D',CENTR,coorCENT',R) ;

COORface = COORtransf(:,NODES_face); 