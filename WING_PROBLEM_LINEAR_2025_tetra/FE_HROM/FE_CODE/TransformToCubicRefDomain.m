function  COORtransf = TransformToCubicRefDomain(COORref,dX,ANGini,xini,yini,a) 
% See NEW_IDEAS.pdf, 
% Cubic mapping of a rectangular domain to a cubically curved domain 
% ------------------------------------------------------------------
%  JAHO, 14-OCT-2020

% Coordinate reference domain (relative centroid face 1)
X0 = COORref(:,1); % linspace(0,dX,nstep) ;
Y0_X0 = COORref(:,2);
% --------------------------------------
% Shrink/elongates so that max(X0) = dX
xmin  = min(X0) ; 
xmax = max(X0)  ; 
L = xmax-xmin ; 
X0 = dX/L*X0 ;  % 
% -----------------------------------------------
% Equation of the deformed midline 
yX0 = polyval(a,X0) ;
% Slope of the midline 
aDER = polyder(a) ;
theta = atan(polyval(aDER,X0)) ;
% --------------------------------------------
% Coordinates domain (rotation around midlines)
X = X0 -sin(theta).*Y0_X0 ; 
Y = yX0 +cos(theta).*Y0_X0 ; 
% Rotation using rotation matrix of the domain 
R = [cos(ANGini) -sin(ANGini); sin(ANGini), cos(ANGini)] ;
XYrot = R*[X';Y'] ;
 

%  translation 
X = XYrot(1,:) + xini ;
Y = XYrot(2,:) + yini ; %f(xini) ;

COORtransf=  [X',Y'] ; 
