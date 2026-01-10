function   DATA3D  = CurvedDomains_beams(DATA3D,angDOM)

if nargin == 0
    load('tmp2.mat')
end

%% Coordinates centroids f1 and f3
Cf1 = DATA3D.CENTRf1 ;
Cf2= DATA3D.CENTRf2;
% Coordinates with respect Cf1
% Rotation angle face 1 with respect to face 3

angDOM  = angDOM/2 ;

% LENGTH
L = Cf2(1)-Cf1(1) ;

% Curved coordinates  (all nodes)
% ----------------------------------------------------
ISBEAM = 1 ; 

DATA3D.COOR  = CurvedCoordinates(DATA3D.COOR',Cf1(:),L,angDOM,ISBEAM)' ;

% Curved coordinates  (centroids)
% ----------------------------------------------------
DATA3D.CENTRf1  = CurvedCoordinates( DATA3D.CENTRf1(:),Cf1(:),L,angDOM,ISBEAM)' ;
DATA3D.CENTRf2  = CurvedCoordinates( DATA3D.CENTRf2(:),Cf1(:),L,angDOM,ISBEAM)' ;
DATA3D.CENTRf1_real  = CurvedCoordinates( DATA3D.CENTRf1_real(:),Cf1(:),L,angDOM,ISBEAM)' ;
DATA3D.CENTRf2_real  = CurvedCoordinates( DATA3D.CENTRf2_real(:),Cf1(:),L,angDOM,ISBEAM)' ;

DATA3D.rotDOMfacesLOC = cell(1,2) ;

ifaceROT = 2;

if length(Cf1) == 3
R = [cos(2*angDOM)   sin(2*angDOM)  0
    -sin(2*angDOM) cos(2*angDOM)   0
    0            0        1 ] ;
else
    R = [cos(2*angDOM) sin(2*angDOM)
     -sin(2*angDOM)       cos(2*angDOM)] ;
end

DATA3D.rotDOMfacesLOC{ifaceROT} = R ;

DATA3D.ROTATIONf1f2= R ;



