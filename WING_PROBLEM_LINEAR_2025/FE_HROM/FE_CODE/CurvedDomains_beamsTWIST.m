function   DATA3D  = CurvedDomains_beamsTWIST(DATA3D,angDOM)

if nargin == 0
    load('tmp1.mat')
end

%% Coordinates centroids f1 and f3
Cf1 = DATA3D.CENTRf1 ;
Cf2= DATA3D.CENTRf2;
 
% LENGTH
L = Cf2(1)-Cf1(1) ;

% Curved coordinates  (all nodes)
% ----------------------------------------------------
ISBEAM = 1 ; 

DATA3D.COOR  = CurvedCoordinatesTWIST(DATA3D.COOR',Cf1,L,angDOM)' ;

% Curved coordinates  (centroids)
% ----------------------------------------------------
DATA3D.CENTRf1  = CurvedCoordinatesTWIST( DATA3D.CENTRf1,Cf1,L,angDOM)' ;
DATA3D.CENTRf2  = CurvedCoordinatesTWIST( DATA3D.CENTRf2,Cf1,L,angDOM)' ;
DATA3D.CENTRf1_real  = CurvedCoordinatesTWIST( DATA3D.CENTRf1_real,Cf1,L,angDOM)' ;
DATA3D.CENTRf2_real  = CurvedCoordinatesTWIST( DATA3D.CENTRf2_real,Cf1,L,angDOM)' ;

DATA3D.rotDOMfacesLOC = cell(1,2) ;

ifaceROT = 2;

R = [1    0   0 
     0   cos(angDOM)  -sin(angDOM)  
     0  sin(angDOM) cos(angDOM)    ] ;
 

DATA3D.rotDOMfacesLOC{ifaceROT} = R ;

DATA3D.ROTATIONf1f2= R ;



