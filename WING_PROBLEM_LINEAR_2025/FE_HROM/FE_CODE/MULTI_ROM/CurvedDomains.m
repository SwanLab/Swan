function   DATA3D  = CurvedDomains(DATA3D,angDOM)

if nargin == 0
    load('tmp.mat')
end

%% Coordinates centroids f1 and f3
Cf1 = DATA3D.COOR_MIDSIDE_FACE{1}(:) ;
Cf3 = DATA3D.COOR_MIDSIDE_FACE{3}(:) ;
% Coordinates with respect Cf1
% Rotation angle face 1 with respect to face 3

angDOM  = angDOM/2 ; 

% LENGTH
L = Cf3(1)-Cf1(1) ;

% Curved coordinates  (all nodes)
% ----------------------------------------------------
DATA3D.COOR  = CurvedCoordinates(DATA3D.COOR',Cf1,L,angDOM)' ;

% Curved coordinates  (centroids)
% ----------------------------------------------------
for iface  =1:length(DATA3D.CENTRf)
    DATA3D.CENTRf{iface}  = CurvedCoordinates( DATA3D.CENTRf{iface}',Cf1,L,angDOM)' ;
end
for iface  =1:length(DATA3D.COOR_CORNER)
    DATA3D.COOR_CORNER{iface}  = CurvedCoordinates( DATA3D.COOR_CORNER{iface},Cf1,L,angDOM) ;
end
for iface  =1:length(DATA3D.COOR_MIDSIDE_FACE)
    DATA3D.COOR_MIDSIDE_FACE{iface}  = CurvedCoordinates( DATA3D.COOR_MIDSIDE_FACE{iface},Cf1,L,angDOM) ;
end

DATA3D.CORNER_1_CENTROID   = CurvedCoordinates( DATA3D.CORNER_1_CENTROID,Cf1,L,angDOM) ;

DATA3D.rotDOMfacesLOC = cell(1,length(DATA3D.COOR_MIDSIDE_FACE)) ;

ifaceROT = 3;
R = [cos(2*angDOM), 0 , sin(2*angDOM)
    0            1        0
    -sin(2*angDOM)  0     cos(2*angDOM)] ;

DATA3D.rotDOMfacesLOC{ifaceROT} = R ;

 
 


