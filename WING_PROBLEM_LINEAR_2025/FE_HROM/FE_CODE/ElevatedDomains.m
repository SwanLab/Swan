function   DATA3D  = ElevatedDomains(DATA3D,DATA)

if nargin == 0
    load('tmp.mat')
end

%% Coordinates centroids f1 and f3
Cf1 = DATA3D.CENTRf1(:) ;
Cf2= DATA3D.CENTRf2(:);
 
% LENGTH
L = Cf2(1)-Cf1(1) ;

% Curved coordinates  (all nodes)
% ----------------------------------------------------
H = DATA.ELEVATION_Z ; 
DATA3D.COOR  = ElevatedZlocal(DATA3D.COOR',Cf1(:),L,H)' ;

% Curved coordinates  (centroids)
% ----------------------------------------------------
DATA3D.CENTRf1  = ElevatedZlocal( DATA3D.CENTRf1(:),Cf1(:),L,H)' ;
DATA3D.CENTRf2  = ElevatedZlocal( DATA3D.CENTRf2(:),Cf1(:),L,H)' ;
DATA3D.CENTRf1_real  = ElevatedZlocal( DATA3D.CENTRf1_real(:),Cf1,L,H)' ;
DATA3D.CENTRf2_real  = ElevatedZlocal( DATA3D.CENTRf2_real(:),Cf1,L,H)' ;
 

