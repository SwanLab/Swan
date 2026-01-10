function  [COORtransf,CENTRf1,CENTRf2 ]= TransformToCubicRefDomain3D(COORref,dP_1,R_G1,ay,az,P1_G,ax,DATAcubic) 
% See NEW_IDEAS.pdf, 
% Cubic mapping of a rectangular domain to a cubically curved domain 
% 3D version. 2D version --> 
% ------------------------------------------------------------------
%  JAHO, 14-OCT-2020/25-Oct-2020 (KAFETERIA, BELGRADE)
if nargin == 0
    load('tmp2.mat')
end

% Coordinate reference domain (relative centroid face 1)
X0 = COORref(:,1); % linspace(0,dX,nstep) ;
Y0 = COORref(:,2);
Z0 = COORref(:,3);
% --------------------------------------
% Shrink/elongates so that max(X0) = dX
xmin  = min(X0) ; 
xmax  = max(X0)  ; 
L = xmax-xmin ;  % Length of the reference domain ---along the x-coordindate. 
dX = dP_1(1) ; % Projection of the midline of the deformed slice along the x axis
X0 = dX/L*X0 ;  % % Stretch/elongate  the undeformed slice, so that the x- coordinates of the centroid of face

%%% CHANGE IN CROSS-SECTIONAL AREA
% -----------------------------------
% Transform Y0 > 0
eLOC_1 = DATAcubic.efact(1).yp ; 
eLOC_2 = DATAcubic.efact(2).yp ;
IND = find(Y0>0) ; % Indexes y >0
xLOC = X0(IND) ; 
eLOC = eLOC_1 + (eLOC_2-eLOC_1)/dX*xLOC ; 
Y0(IND) = Y0(IND).*eLOC ; 
% Transform Y0 < 0
eLOC_1 = DATAcubic.efact(1).yn ; 
eLOC_2 = DATAcubic.efact(2).yn ;
IND = find(Y0<0) ; % Indexes y >0
xLOC = X0(IND) ; 
eLOC = eLOC_1 + (eLOC_2-eLOC_1)/dX*xLOC ; 
Y0(IND) = Y0(IND).*eLOC ; 
% Transform Z0 > 0
eLOC_1 = DATAcubic.efact(1).zp ; 
eLOC_2 = DATAcubic.efact(2).zp ;
IND = find(Z0>0) ; % Indexes y >0
xLOC = X0(IND) ; 
eLOC = eLOC_1 + (eLOC_2-eLOC_1)/dX*xLOC ; 
Z0(IND) = Z0(IND).*eLOC ; 
% Transform Y0 < 0
eLOC_1 = DATAcubic.efact(1).zn ; 
eLOC_2 = DATAcubic.efact(2).zn ;
IND = find(Z0<0) ; % Indexes y >0
xLOC = X0(IND) ; 
eLOC = eLOC_1 + (eLOC_2-eLOC_1)/dX*xLOC ; 
Z0(IND) = Z0(IND).*eLOC ;




 % 2 is   dP_1(1)
% -----------------------------------------------
% Points of the deformed midline 
yX0 = polyval(ay,X0)   ;  % Position of the points after applying the displacement dictated by the midline 
zX0 = polyval(az,X0)   ;
% ----------------------------------
CENTRf1 = [0,0,0] ; 
CENTRf2 = [dX,polyval(ay,dX),polyval(az,dX)] ;   % Position of the centroids 
% Now we have to calculate the rotations of each point. Since the
% cross-sections are assumed to be remain normal to the midline, the
% rotation should be given by the normal to the curve 
%
%  x = t 
%  y = T*ay
%  z = T*az
% where  T = [t^3 t^2 t 1] ; 
% Therefore 
% x' =  1 
% y' = T'*ay  = polyval(ay,X0) 
% z' = T'*az  = polyval(az,X0)
%
% We know that %-----------------------------------------------------------
ayD = polyder(ay) ; 
azD = polyder(az) ; 
yD = polyval(ayD,X0) ; 
zD = polyval(azD,X0) ;  
% Vector normal to the curve = [1 yD zD]'
% Rotation matrix determined from this vector 
% We know that (see SymbolicRotations.m) that 
%R =
 
% R =[cos(ay)*cos(az), -cos(ay)*sin(az), -sin(ay)]
%[        sin(az),          cos(az),        0]
%[cos(az)*sin(ay), -sin(ay)*sin(az),  cos(ay)]
% and that 
% R(:,1) = [1 ayD azD]'/norm([1 ayD azD]')
% Thus 
MODnVECT = sqrt(1 + yD.^2 + zD.^2) ;  % Norm of the vector tangent to the curve 
nVECT{1} = 1./MODnVECT ;   % Components of the unit vector tangent to the curve 
nVECT{2} = yD./MODnVECT ; 
nVECT{3} = zD./MODnVECT ; 
% Therefore 
sinAZ = nVECT{2}  ; 
cosAZ = sqrt(nVECT{1}.^2 + nVECT{3}.^2)  ; %s cos(asin(sinAZ)) ; 
cosAY = nVECT{1}./cosAZ  ; 
sinAY =  nVECT{3}./cosAZ  ;
%%% The 6 entries of the rotation matrix (at each point of the curve) are therefore given by 
if ~isempty(ax)
    axLOC = X0/max(X0)*ax; 
    cosAX = cos(axLOC); 
    sinAX = sin(axLOC) ; 
else
     cosAX = ones(size(cosAZ))  ; 
 sinAX = zeros(size(cosAZ))  ; 
end

%   sin_az = nt2_1(2)  ;
%             cos_az = sqrt(nt2_1(1)^2 + nt2_1(3)^2)  ;
%             sin_ay =  nt2_1(3)/cos_az ;
%             cos_ay =  nt2_1(1)/cos_az ;
%             az = real(acos(cos_az)) ;
%             ay = real(acos(cos_ay)) ;
%             ax = 0 ;

% [cos(ay)*cos(az), - sin(ax)*sin(ay) - cos(ax)*cos(ay)*sin(az), cos(ay)*sin(ax)*sin(az) - cos(ax)*sin(ay)]
% [        sin(az),                             cos(ax)*cos(az),                          -cos(az)*sin(ax)]
% [cos(az)*sin(ay),   cos(ay)*sin(ax) - cos(ax)*sin(ay)*sin(az), cos(ax)*cos(ay) + sin(ax)*sin(ay)*sin(az)]
%  
R = cell(3,3) ; 
R{1,1} = cosAY.*cosAZ ; 
R{2,1} = sinAZ ;
R{3,1} = cosAZ.*sinAY ; 
%-----------------------
R{1,2} = -cosAY.*sinAZ.*cosAX - sinAX.*sinAY ; 
R{2,2} = cosAZ.*cosAX ;
R{3,2} = -sinAZ.*sinAY.*cosAX + sinAX.*cosAY ; 
%-----------------------
R{1,3} = -sinAY.*cosAX + cosAY.*sinAX.*sinAZ ; 
R{2,3} = -cosAZ.*sinAX ;
R{3,3} = cosAY.*cosAX + sinAX.*sinAY.*sinAZ ; 

%%%
% Coordinates domain (rotation around midline) -- R*[0,Y0,Z0]'; 
X_1 = X0  + R{1,2}.*Y0 + R{1,3}.*Z0  ; %   % EXPRESSED IN THE REFERENCE SYSTEM AT CENTROID 1 
Y_1 = yX0 + R{2,2}.*Y0 + R{2,3}.*Z0 ;  %
Z_1 = zX0   + R{3,2}.*Y0 + R{3,3}.*Z0 ;  % %
% --------------------------------------------
%X = X0 -sin(theta).*Y0_X0 ; 
%Y = yX0 +cos(theta).*Y0_X0 ; 


% Rotation using rotation matrix of the domain 
%R = [cos(ANGini) -sin(ANGini); sin(ANGini), cos(ANGini)] ;
COORtransf = R_G1*[X_1';Y_1';Z_1'] ;  % global coordinates 
% Translation 
 COORtransf = COORtransf' ; 

%  translation 
dim = 3 ; 
for idim = 1:dim
    COORtransf(:,idim) = COORtransf(:,idim)' + P1_G(idim) ; 
end
 