%==================================================================
%                        General Data File
% Title: Default_title
% Units: SI
% Dimensions: 2D
% Type of problem: Plane_Stress
% Type of Phisics: Stokes
% Micro/Macro: MACRO
%
%==================================================================

%% Data
order={'QUADRATIC','LINEAR'};
Data_prb = {
'TRIANGLE';
'SI';
'2D';
'Plane_Stress';
'Stokes';
'MACRO'
};
state = 'Transient';
dtime = 0.01;
finalTime = 1;

%% Coordinates
% Node                X                Y                Z

coord = [
1            0            0            0
2            0         0.25            0
3         0.25            0            0
4         0.25         0.25            0
5            0          0.5            0
6          0.5            0            0
7         0.25          0.5            0
8          0.5         0.25            0
9          0.5          0.5            0
10         0.75            0            0
11            0         0.75            0
12         0.75         0.25            0
13         0.25         0.75            0
14         0.75          0.5            0
15          0.5         0.75            0
16            1            0            0
17            0            1            0
18            1         0.25            0
19         0.25            1            0
20         0.75         0.75            0
21            1          0.5            0
22          0.5            1            0
23         0.75            1            0
24            1         0.75            0
25            1            1            0
];

%% Conectivities
% Element        Node(1)                Node(2)                Node(3)                Material

connec = [
1 6 8 3 0
2 8 9 4 0
3 3 4 1 0
4 8 4 3 0
5 9 7 4 0
6 7 5 2 0
7 4 2 1 0
8 7 2 4 0
9 16 18 10 0
10 18 21 12 0
11 10 12 6 0
12 18 12 10 0
13 21 14 12 0
14 14 9 8 0
15 12 8 6 0
16 14 8 12 0
17 9 15 7 0
18 15 22 13 0
19 7 13 5 0
20 15 13 7 0
21 22 19 13 0
22 19 17 11 0
23 13 11 5 0
24 19 11 13 0
25 21 24 14 0
26 24 25 20 0
27 14 20 9 0
28 24 20 14 0
29 25 23 20 0
30 23 22 15 0
31 20 15 9 0
32 23 15 20 0
];





%% Variable Prescribed
% Node            Dimension                Value

isLeft   = @(coor) (abs(coor(:,1) - min(coor(:,1)))   < 1e-12);
isRight  = @(coor) (abs(coor(:,1) - max(coor(:,1)))   < 1e-12);
isBottom = @(coor) (abs(coor(:,2) - min(coor(:,2)))   < 1e-12);
isTop    = @(coor) (abs(coor(:,2) - max(coor(:,2)))   < 1e-12);

% Dirichlet

% Velocity...
sDir{1}.domain    = @(coor) isLeft(coor) | isRight(coor) | isTop(coor) | isBottom(coor);
sDir{1}.direction = [1,2];
sDir{1}.value     = 0;
% Pressure...
sDir{2}.domain    = @(coor) isRight(coor) & isTop(coor);
sDir{2}.direction = 1;
sDir{2}.value     = 0;

velocityBC.domain = 'Border';
velocityBC.value  = 0;

velocity = [
1 1 0 
1 2 0 
5 1 0 
5 2 0 
6 1 0 
6 2 0 
16 1 0 
16 2 0 
17 1 0 
17 2 0 
21 1 0 
21 2 0 
22 1 0 
22 2 0 
25 1 0 
25 2 0 
];

pressure = [
25 1 0 
];

nu=1;

sMesh.coord = coord(:,2:3);
sMesh.connec = connec(:,2:4);
mesh = Mesh.create(sMesh);

sAF.fHandle = @(coor) [nu*4*(coor(1,:,:).^3.*((6 - 12*coor(2,:,:))) + coor(1,:,:).^4.*((-3 + 6*coor(2,:,:))) + coor(2,:,:).*(1 - 3*coor(2,:,:) + ...
    2*coor(2,:,:).^2)-6*coor(1,:,:).*coor(2,:,:).*(1 - 3*coor(2,:,:) + 2*coor(2,:,:).^2)...
+ 3*coor(1,:,:).^2.*((-1 + 4*coor(2,:,:) - 6*coor(2,:,:).^2 + 4*coor(2,:,:).^3)))+(2*coor(1,:,:)-1).*(coor(2,:,:)-1);
-4*nu*(-3*((-1 + coor(2,:,:)).^2).*coor(2,:,:).^2 - 3*coor(1,:,:).^2.*((1 - 6*coor(2,:,:) + 6*coor(2,:,:).^2))+2*coor(1,:,:).^3.*((1 - 6*coor(2,:,:) + 6*coor(2,:,:).^2)) +...
coor(1,:,:).*(1 - 6*coor(2,:,:) + 12*coor(2,:,:).^2 - 12*coor(2,:,:).^3 + 6*coor(2,:,:).^4))+coor(1,:,:).*(coor(1,:,:)-1)];
sAF.ndimf   = 2;
sAF.mesh    = mesh;
Vol_force = AnalyticalFunction(sAF);

% nu=1;
% Vol_force = @(x,y){nu*4*(x^3*((6 - 12*y)) + x^4*((-3 + 6*y)) + y*(1 - 3*y + 2*y^2)-6*x*y*(1 - 3*y + 2*y^2)...
% + 3*x^2*((-1 + 4*y - 6*y^2 + 4*y^3)))+nu*y*(1 - y)*(1 - 2*x);
% -4*nu*(-3*((-1 + y)^2)*y^2 - 3*x^2*((1 - 6*y + 6*y^2))+2*x^3*((1 - 6*y + 6*y^2)) +...
% x*(1 - 6*y + 12*y^2 - 12*y^3 + 6*y^4)+x*(1 - x)*(1 - 2*y))};





%% Group Elements
% Element        Group_num

Group = [
];

%% Initial Holes
% Elements that are considered holes initially
% Element

Initial_holes = [
];

%% Boundary Elements
% Elements that can not be removed
% Element

Boundary_elements = [
];

%% Micro gauss post
%
% Element

Micro_gauss_post = [
];


%% Micro Slave-Master
% Nodes that are Slaves
% Nodes             Value (1-Slave,0-Master)

Micro_slave = [
];

%% Nodes solid
% Nodes that must remain 
% Nodes

% nodesolid = unique(pointload_complete(:,1));

%% External border Elements
% Detect the elements that define the edge of the domain
% Element               Node(1)           Node(2)

External_border_elements = [
1 3 6
3 1 3
6 5 2
7 2 1
9 16 18
9 10 16
10 18 21
11 6 10
21 22 19
22 17 11
22 19 17
23 11 5
25 21 24
26 24 25
29 25 23
30 23 22
];

%% External border Nodes
% Detect the nodes that define the edge of the domain
% Node

External_border_nodes = [
];

%% Materials
% Materials that have been used
% Material_Num              Mat_density        Young_Modulus        Poisson

Materials = [
];
