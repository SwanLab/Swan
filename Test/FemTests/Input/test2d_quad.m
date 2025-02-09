%==================================================================
%                        General Data File
% Title: QUAD
% Units: SI
% Dimensions: 2D
% Type of problem: Plane_Stress
% Type of Phisics: ELASTIC
% Micro/Macro: MACRO
%
%==================================================================

%% Data

Data_prb = {
'QUAD';
'SI';
'2D';
'Plane_Stress';
'ELASTIC';
'MACRO';
};

%% Coordinates
% Node                X                Y                Z

coord = [
1            2            0            0
2          1.5            0            0
3            2          0.5            0
4          1.5          0.5            0
5            1            0            0
6            2            1            0
7          1.5            1            0
8            1          0.5            0
9            1            1            0
10          0.5            0            0
11          0.5          0.5            0
12          0.5            1            0
13            0            0            0
14            0          0.5            0
15            0            1            0
];

%% Conectivities
% Element        Node(1)                Node(2)                Node(3)                Node(4)                Material

connec = [
1 10 11 14 13 0
2 5 8 11 10 0
3 2 4 8 5 0
4 1 3 4 2 0
5 11 12 15 14 0
6 8 9 12 11 0
7 4 7 9 8 0
8 3 6 7 4 0
];

%% Variable Prescribed
% Node            Dimension                Value

dirichlet_data = [
13 1 0 
13 2 0 
14 1 0 
14 2 0 
15 1 0 
15 2 0 
];

%% Force Prescribed
% Node                Dimension                Value

pointload_complete = [
3 1 0
3 2 -1 
];


isLeft   = @(coor) (abs(coor(:,1) - min(coor(:,1)))   < 1e-12);
isRight  = @(coor) (abs(coor(:,1) - max(coor(:,1)))   < 1e-12);
isMiddle = @(coor) (abs(coor(:,2) - max(coor(:,2)/2)) == 0);

% Dirichlet
sDir{1}.domain    = @(coor) isLeft(coor);
sDir{1}.direction = [1,2];
sDir{1}.value     = 0;

% Point load
sPL{1}.domain    = @(coor) isMiddle(coor) & isRight(coor);
sPL{1}.direction = 2;
sPL{1}.value     = -1;

%% Volumetric Force
% Element        Dim                Force_Dim

Vol_force = [
];

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

nodesolid = unique(pointload_complete(:,1));

%% External border Elements
% Detect the elements that define the edge of the domain
% Element               Node(1)           Node(2)

External_border_elements = [
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
