function [DATA3D] =    SliceNodesFacesIdentificationRVE2D(DATA3D,NORMALSv,DATAIN)
% IDENTIFICATION OF NODES FACE F1, F2, F3 and F4
% ---------------------------------------
% % Face "f1"  is assumed to be normal to the x-axis . The geometry of the
% SLICE should be drawn in GID so that face f1 is parallel to the plane
% x=0.  Furthermore, the centroid of face f2 should have a positive x
% coordinate.
% The same applies for faces 3 and f
% JAHO 22TH jANUARY 2018
% -------------------------
if nargin == 0
    load('tmp1.mat')
end
% First we check whether the set of faces have been previously defined from
% GID's pre-processesor
% ----------------------------------------------------------------------------
% Local tolerance
% --------------------

%DATAIN = DefaultField(DATAIN,'TOLERANCE_PERIODIC_NODES',1e-1) ; 
%TOL = DATAIN.TOLERANCE_PERIODIC_NODES  ; 
TOL = 1e-1 ; % Local tolerance
CNb1 = DATA3D.CNb(1,1:2) ;
COOR1  = DATA3D.COOR(CNb1,:) ;
dCOOR = norm(COOR1(1,:)-COOR1(2,:)) ;
TOL = TOL*dCOOR ; % Abs. tolerance

DATAIN = DefaultField(DATAIN,'TOLERANCE_PERIODIC_NODES_GIVEN',TOL) ; 
TOL = DATAIN.TOLERANCE_PERIODIC_NODES_GIVEN ; 

% Computing the exact  centroids
CENTRf = ComputeCentroid(DATA3D,DATA3D.NODES_FACES) ;
DATA3D.CENTRf = CENTRf ; 
% FACES 1 and 2 
% --------------
face1 = 1; face2  =3 ;  
NODESface2_ord = SortingNodesFaces(face1,face2,DATA3D,CENTRf,NORMALSv,TOL,DATAIN) ; 
DATA3D.NODES_FACES{face2} = NODESface2_ord ; 


% FACES 3 and 4
% ----------------
face1 = 2; face2  =4 ;  
NODESface2_ord = SortingNodesFaces(face1,face2,DATA3D,CENTRf,NORMALSv,TOL,DATAIN) ; 
DATA3D.NODES_FACES{face2} = NODESface2_ord ; 




 
 



 

end


function NODESface2_ord = SortingNodesFaces(face1,face2,DATA3D,CENTRf,NORMALSv,TOL,DATAIN)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the proviso that face F1 should be parallel to the YZ plane
if face1 == 1 
    LEG = 'YZ' ;
    LEG2  = 'x' ;
    DIRNORMAL = 1; 
else
    LEG = 'XZ' ; 
    LEG2 = 'y' ; 
    DIRNORMAL = 2; 
end
[CNloc1 ELEMloc] = ElemBnd(DATA3D.CNb, DATA3D.NODES_FACES{face1}) ;
NORMAL_f1 = NORMALSv(:,ELEMloc(1)) ;

if abs(abs(NORMAL_f1(DIRNORMAL))-1)>TOL
    error(['Face ',num2str(face1),'should be parallel to the ',LEG,' plane. Rotate GID s geometry'])
end
% Face 1 and 2 
if CENTRf{2}(1)-CENTRf{1}(1) <0
    error(['Face f',num2str(face2),' should have a   ',LEG2,' coordinate greater than f',num2str(face1)])
end
DATAIN = DefaultField(DATAIN,'InterfacesDoNotMatch',0) ;
if (length(DATA3D.NODES_FACES{face1}) ~= length(DATA3D.NODES_FACES{face2})) &  DATAIN.InterfacesDoNotMatch==0
    error('Non-conforming meshes. Faces ',  num2str(face1),  ' and ',  num2str(face2),' has different number of nodes')
end

%%% SORTING NODES 
% -------------------------- 

% Order NODES_face1 and NODES_face2 so that the mapping is
% COOR(NODES_face1(1)) -->  COOR(NODES_face2(1))
TRANSLATION = CENTRf{face2}-CENTRf{face1} ;
COORf2 = DATA3D.COOR(DATA3D.NODES_FACES{face2},:) ;
COORf1 = DATA3D.COOR(DATA3D.NODES_FACES{face1},:) ;

 %%%%%
%%% Relative coordinates
COORf1_rel = zeros(size(COORf1)) ;
for idim = 1:length(TRANSLATION)
    COORf1_rel(:,idim) = COORf1(:,idim) - CENTRf{face1}(idim) ;
end
% TRANSLATION
COORf1_trans = zeros(size(COORf1)) ;
for idim = 1:length(TRANSLATION)
    COORf1_trans(:,idim) = COORf1_rel(:,idim) + CENTRf{face1}(idim) + TRANSLATION(idim) ;
end

if DATAIN.InterfacesDoNotMatch==0
    [IDX DISTANCES]= knnsearch(COORf2,COORf1_trans) ;
    if any(find(DISTANCES > TOL))
        disp('Non-conforming meshes. Check tolerance TOL')
        
        error(['Try with a DATAIN.TOLERANCE_PERIODIC_NODES_GIVEN  ',num2str(max(DISTANCES))])
    end
    NODESface2_ord = DATA3D.NODES_FACES{face2}(IDX) ;
else
    NODESface2_ord = DATA3D.NODES_FACES{face2} ;
end

end
