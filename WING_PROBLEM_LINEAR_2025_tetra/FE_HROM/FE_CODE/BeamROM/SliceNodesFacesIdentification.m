function [NODES_face1,NODESface2_ord,CENTRf1,CENTRf2,ROTATIONf1f2,CENTRf1_real,CENTRf2_real] = ...
    SliceNodesFacesIdentification(DATA3D,NORMALSv,DATAIN)
% IDENTIFICATION OF NODES FACE F1 AND F2
% ---------------------------------------
% % Face "f1"  is assumed to be normal to the x-axis . The geometry of the
% SLICE should be drawn in GID so that face f1 is parallel to the plane
% x=0.  Furthermore, the centroid of face f2 should have a positive x
% coordinate.
% JAHO 22TH jANUARY 2018
% -------------------------
if nargin == 0
    load('tmp0.mat')
end
% First we check whether the set of faces have been previously defined from
% GID's pre-processesor
% ----------------------------------------------------------------------------
% Local tolerance
% --------------------
TOL = 1e-1 ; % Local tolerance
DATAIN = DefaultField(DATAIN,'TOL_face_identification_NODES',TOL) ; 
TOL = DATAIN.TOL_face_identification_NODES ; 
CNb1 = DATA3D.CNb(1,1:2) ;
COOR1  = DATA3D.COOR(CNb1,:) ;
dCOOR = norm(COOR1(1,:)-COOR1(2,:)) ;
TOL = TOL*dCOOR ; % Abs. tolerance
% ------------------------------------------------------------------------------
% Identifying nodes face f1
% -----------------------------------------------------------------------------
% -----------------------------------------------------------------------------
if isempty(DATA3D.NODES_FACES) || isempty(DATA3D.NODES_FACES{1})
    MANUAL = 1;
    
    % Elements faces f1
    elements_f12 = find(abs(abs(NORMALSv(1,:))-1)<TOL)  ;
    NODES_faces12 = unique(DATA3D.CNb(elements_f12,:)) ;
    x  = DATA3D.COOR(NODES_faces12,1) ;
    xmin = min(x) ;
    INDMIN =     find(abs(DATA3D.COOR(NODES_faces12,1)-xmin) <TOL  ) ;  % Face xMIN, f1
    NODES_face1 = NODES_faces12(INDMIN) ; % Nodes face xMIN, f1
else
    MANUAL = 0 ;
    NODES_face1 = DATA3D.NODES_FACES{1} ;
end

% Check the proviso that face F1 should be parallel to the YZ plane
% Pick an  element of face f1
[CNloc1 ELEMloc] = ElemBnd(DATA3D.CNb,NODES_face1) ;



NORMAL_f1 = NORMALSv(:,ELEMloc(1)) ;
if abs(abs(NORMAL_f1(1))-1)>TOL
    error(['Face F1 should be parallel to the YZ plane. Rotate GID s geometry'])
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Identifying nodes fACE f2
% --------


if isempty(DATA3D.NODES_FACES) || isempty(DATA3D.NODES_FACES{2})
    if MANUAL == 1
        xmax = max(x) ;
        if xmin==xmax
            error('Face f2 is inclined with respect to face f1. Define FACE f2 using GID facilities ')
        else
            INDMAX =     find(abs(COOR(NODES_faces12,1)-xmax) <TOL  ) ;
            NODES_face2 = NODES_faces12(INDMAX) ;
        end
    else
        error(' Define FACE f2 using GID facilities ')
    end
else
    NODES_face2 = DATA3D.NODES_FACES{2} ;
end

%%%%%%%%%%%%%%%%%%%%%%5


DATAIN = DefaultField(DATAIN,'InterfacesDoNotMatch',0) ;

if (length(NODES_face1) ~= length(NODES_face2)) &  DATAIN.InterfacesDoNotMatch==0
    disp('Non-conforming meshes. Face f1 and f2 has different number of nodes')
    error('Check GID geometry. Perhaps faces have not been correctly assigned')
end

% Computing the centroid

NODESfaces{1} = NODES_face1;
NODESfaces{2} = NODES_face2;
CENTROID = ComputeCentroid(DATA3D,NODESfaces) ;
CENTRf1_real = CENTROID{1}(:) ;
CENTRf2_real = CENTROID{2}(:) ;

CENTRf1 =CENTRf1_real ; % sum(DATA3D.COOR(NODES_face1,:)',2)/length(NODES_face1) ; % (psesudo-)Centroid face1
CENTRf2 = CENTRf2_real ; %sum(DATA3D.COOR(NODES_face2,:)',2)/length(NODES_face2) ;  % Centroid face2

% PSEUDO_CENTROID =1;
%
% if PSEUDO_CENTROID == 1

%% THIS WORK IS DONE TWICE !!! RE-PROGRAM THIS PART ---2-July-2018


%end



if CENTRf2(1)-CENTRf1(1) <0
    error('Face f2 should have a   x coordinate greater than f1')
end

% Order NODES_face1 and NODES_face2 so that the mapping is
% COOR(NODES_face1(1)) -->  COOR(NODES_face2(1))
TRANSLATION = CENTRf2-CENTRf1 ;
COORf2 = DATA3D.COOR(NODES_face2,:) ;
COORf1 = DATA3D.COOR(NODES_face1,:) ;

% DETERMINATION OF ROTATION MATRIX  BETWEEN FACES(CURVED ELEMENTS)
% Pick an element of face 2
[CNloc, ELEMloc] = ElemBnd(DATA3D.CNb,NODES_face2) ;

normalF2 = NORMALSv(:,ELEMloc(1)) ;
if normalF2(1)<0
    normalF2 = -normalF2 ;
end
% Rotation matrix
% First vector  --> normalF2
r1 = normalF2 ;

%r1 = (CENTRf2-CENTRf1)/norm(CENTRf2-CENTRf1) ;

if size(DATA3D.COOR,2) == 3
    r3 = [0 0 1]' ;
    % Second vector  is perpendicular to r1 and nr
    r2 = cross(r3,r1) ;
    % Third vector
    r3 = cross(r1,r2) ;
    ROTATIONf1f2 = [r1,r2,r3] ; % Rotation matrix
else
    r2 = r1 ;
    r2(1) = -r1(2) ;
    r2(2) = r1(1) ;
    ROTATIONf1f2 = [r1,r2] ; % Rotation matrix
    
    
end

%%%% Rotation center face2 with respect to face 1 (empty for straight slices)
% if abs(abs(r1(1))-1)>0
% n1 = [1 0]' ;
% n2 = r1(1:2) ;
% c1 = CENTRf1(1:2) ;
% c2 = CENTRf2(1:2) ;
% A = [n1';n2'] ;
% b = [n1'*c1;n2'*c2] ;
% Crotation = A\b ;
% else
%     Crotation = [] ;
% end
%Crotation = [] ;

% -----------------------------------

%%%%%
%%% Relative coordinates
COORf1_rel = zeros(size(COORf1)) ;
for idim = 1:length(TRANSLATION)
    COORf1_rel(:,idim) = COORf1(:,idim) - CENTRf1(idim) ;
end
% Rotation
COORf1_rel = (ROTATIONf1f2*COORf1_rel')' ;
% TRANSLATION
COORf1_trans = zeros(size(COORf1)) ;
for idim = 1:length(TRANSLATION)
    COORf1_trans(:,idim) = COORf1_rel(:,idim) + CENTRf1(idim) + TRANSLATION(idim) ;
end

if DATAIN.InterfacesDoNotMatch==0
    [IDX DISTANCES]= knnsearch(COORf2,COORf1_trans) ;
    if any(find(DISTANCES > TOL))
        error('Non-conforming meshes. Check tolerance TOL')
    end
    NODESface2_ord = NODES_face2(IDX) ;
else
    NODESface2_ord = NODES_face2 ;
end


% %%%% DETECTING BOUNDARY NODES OF either NODES_face1
% % Local connectivity matrix --> CNloc1
% CN_BoundaryFace1 =  IdentifyBoundaryMesh(CNloc1) ;
% % FAce 2  CNloc
% CN_BoundaryFace2 =  IdentifyBoundaryMesh(CNloc) ;


