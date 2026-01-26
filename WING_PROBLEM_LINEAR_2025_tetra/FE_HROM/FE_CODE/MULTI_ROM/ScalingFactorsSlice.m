function   DATA3D  = ScalingFactorsSlice(DATA3D,SCALE_FACTOR)

if nargin == 0
    load('tmp2.mat')
end

%% Coordinates centroids f1 
Cf1 = DATA3D.CENTRf1 ; % ---------

DATA3D.COOR= ScalingSliceLoc(DATA3D.COOR,Cf1,SCALE_FACTOR) ; 
DATA3D.CENTRf2= ScalingSliceLoc(DATA3D.CENTRf2(:)',Cf1,SCALE_FACTOR) ; 
DATA3D.CENTRf1_real= ScalingSliceLoc(DATA3D.CENTRf1_real(:)',Cf1,SCALE_FACTOR) ; 
DATA3D.CENTRf2_real= ScalingSliceLoc(DATA3D.CENTRf2_real(:)',Cf1,SCALE_FACTOR) ; 
 

%  
% % Curved coordinates  (centroids)
% % ----------------------------------------------------
% DATA3D.CENTRf1  = CurvedCoordinates( DATA3D.CENTRf1,Cf1,L,angDOM,ISBEAM)' ;
% DATA3D.CENTRf2  = CurvedCoordinates( DATA3D.CENTRf2,Cf1,L,angDOM,ISBEAM)' ;
% DATA3D.CENTRf1_real  = CurvedCoordinates( DATA3D.CENTRf1_real,Cf1,L,angDOM,ISBEAM)' ;
% DATA3D.CENTRf2_real  = CurvedCoordinates( DATA3D.CENTRf2_real,Cf1,L,angDOM,ISBEAM)' ;
% 
% DATA3D.rotDOMfacesLOC = cell(1,2) ;
% 
% ifaceROT = 2;
% 
% if length(Cf1) == 3
% R = [cos(2*angDOM)   sin(2*angDOM)  0
%     -sin(2*angDOM) cos(2*angDOM)   0
%     0            0        1 ] ;
% else
%     R = [cos(2*angDOM) sin(2*angDOM)
%      -sin(2*angDOM)       cos(2*angDOM)] ;
% end
% 
% DATA3D.rotDOMfacesLOC{ifaceROT} = R ;
% 
% DATA3D.ROTATIONf1f2= R ;
% 
% 
% 
