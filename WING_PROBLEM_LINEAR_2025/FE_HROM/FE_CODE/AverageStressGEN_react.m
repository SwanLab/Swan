function  DATAOUT = AverageStressGEN_react(DATA,COOR,React,DATAOUT) ; 

dbstop('4')
if nargin == 0
    load('tmp.mat')
end

% Area 
xmin = min(COOR(:,1)) ; 
xmax = max(COOR(:,1)) ; 
dx = xmax(1)-xmin(1) ; 
ymin = min(COOR(:,2)) ; 
ymax = max(COOR(:,2)) ; 
zmin = min(COOR(:,3)) ; 
zmax = max(COOR(:,3)) ; 

% %Z%
% X = COOR(:,1) ; Y =COOR(:,2) ; 
% Z = COOR(:,3) ; 
% xREACT = React(1:3:end) ;
% yREACT = React(2:3:end) ;
% zREACT = React(3:3:end) ;
% % Mx 
% sum(Y.*zREACT-Z.*yREACT)
 


%% Resultant forces on plane x=xmax 
nodes = DATA.NODESpl{1,1} ; 
AREA = (ymax(1)-ymin(1))*(zmax(1)-zmin(1)) ;
Z = COOR(:,3) ; Z = Z(nodes);
%%%%
  
% Resultant force along x 
DOF = nodes*3-2 ; 
Nx = sum(React(DOF))/AREA; 
% Mx moment
Mx = sum(React(DOF).*Z)/AREA ; 
% Resultant force along y
DOF = nodes*3-1 ; 
Nxy = sum(React(DOF))/AREA; 
% Mx moment
Mxy = sum(React(DOF).*Z)/AREA ; 

% Resultant force along z
DOF = nodes*3 ; 
Qx = sum(React(DOF))/AREA;

%% Resultant forces on plane x=xmin 
nodes = DATA.NODESpl{1,2} ; 
AREA = (ymax(1)-ymin(1))*(zmax(1)-zmin(1)) ;
Z = COOR(:,3) ; Z = Z(nodes);
% Resultant force along x 
DOF = nodes*3-2 ; 
Nxmin = sum(React(DOF))/AREA; 
% Mx moment
Mxmin = sum(React(DOF).*Z)/AREA ; 
% Resultant force along y
DOF = nodes*3-1 ; 
Nxymin = sum(React(DOF))/AREA; 
% Mx moment
Mxymin = sum(React(DOF).*Z)/AREA ; 

% Resultant force along z
DOF = nodes*3 ; 
Qx = sum(React(DOF))/AREA;
 

%% Resultant forces on plane y=xmax 
nodes = DATA.NODESpl{2,1} ; 
Z = COOR(:,3) ; Z = Z(nodes);
AREA = (xmax(1)-xmin(1))*(zmax(1)-zmin(1)) ;
% % Resultant force along x 
% DOF = nodes*3-2 ; 
% Nxy_2 = sum(React(DOF))/AREA; 
% Resultant force along y
DOF = nodes*3-1 ; 
Ny = sum(React(DOF))/AREA; 
% My moment
My = sum(React(DOF).*Z)/AREA ; 
% Resultant force along z
DOF = nodes*3 ; 
Qy = sum(React(DOF))/AREA;


%% Resultant forces on plane z=zmax 
nodes = DATA.NODESpl{3,2} ; 
Z = COOR(:,3) ; Z = Z(nodes);
AREA = (xmax(1)-xmin(1))*(ymax(1)-ymin(1)) ;
% % Resultant force along x 
 DOF = nodes*3-2 ; 
 Qx_z = sum(React(DOF))/AREA; 
 %Resultant force along y
DOF = nodes*3-1 ; 
Qy_z = sum(React(DOF))/AREA; 
% My moment
Mz = sum(React(DOF).*Z)/AREA ; 
% Resultant force along z
DOF = nodes*3 ; 
Nz = sum(React(DOF))/AREA;


%% Resultant forces on plane z=zmin 
nodes = DATA.NODESpl{3,1} ; 
Z = COOR(:,3) ; Z = Z(nodes);
AREA = (xmax(1)-xmin(1))*(ymax(1)-ymin(1)) ;
% % Resultant force along x 
 DOF = nodes*3-2 ; 
 Qx_zmin = sum(React(DOF))/AREA; 
 %Resultant force along y
DOF = nodes*3-1 ; 
Qy_zmin = sum(React(DOF))/AREA; 
% My moment
Mz = sum(React(DOF).*Z)/AREA ; 
% Resultant force along z
DOF = nodes*3 ; 
Nz = sum(React(DOF))/AREA;

%dbstop('56')
DATAOUT.stressAVG = [Nx Ny Nxy Mx My Mxy Qy Qx]' ; 

