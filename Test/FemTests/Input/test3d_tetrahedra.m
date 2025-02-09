%==================================================================
%                        General Data File
% Title: TETREAHEDRA
% Units: SI
% Dimensions: 3D
% Type of problem: Plane_Stress
% Type of Phisics: ELASTIC
% Micro/Macro: MACRO
%
%==================================================================

%% Data

Data_prb = {
'TETRAHEDRA';
'SI';
'3D';
'Plane_Stress';
'ELASTIC';
'MACRO';
};

%% Coordinates
% Node                X                Y                Z

coord = [
1            2            0            1
2            2            0          0.5
3            2          0.5            1
4            2          0.5          0.5
5            2            0            0
6            2            1            1
7            1            0            1
8            1            0          0.5
9            1          0.5            1
10            2            1          0.5
11            2          0.5            0
12            1          0.5          0.5
13            2            1            0
14            1            0            0
15            1            1            1
16            1          0.5            0
17            1            1          0.5
18            1            1            0
19            0            0            1
20            0            0          0.5
21            0          0.5            1
22            0          0.5          0.5
23            0            1            1
24            0            0            0
25            0          0.5            0
26            0            1          0.5
27            0            1            0
];

%% Conectivities
% Element        Node(1)                Node(2)                Node(3)                Node(4)                Material

connec = [
1 10 17 12 6 0
2 6 10 4 12 0
3 15 9 12 6 0
4 6 3 9 12 0
5 6 15 17 12 0
6 3 4 12 6 0
7 11 16 12 13 0
8 13 11 4 12 0
9 18 17 12 13 0
10 13 10 17 12 0
11 13 18 16 12 0
12 10 4 12 13 0
13 9 21 23 12 0
14 12 9 15 23 0
15 22 26 23 12 0
16 12 17 26 23 0
17 12 22 21 23 0
18 17 15 23 12 0
19 17 26 27 12 0
20 12 17 18 27 0
21 22 25 27 12 0
22 12 16 25 27 0
23 12 22 26 27 0
24 16 18 27 12 0
25 3 9 12 1 0
26 1 3 4 12 0
27 7 8 12 1 0
28 1 2 8 12 0
29 1 7 9 12 0
30 2 4 12 1 0
31 2 8 12 5 0
32 5 2 4 12 0
33 14 16 12 5 0
34 5 11 16 12 0
35 5 14 8 12 0
36 11 4 12 5 0
37 8 20 19 12 0
38 12 8 7 19 0
39 22 21 19 12 0
40 12 9 21 19 0
41 12 22 20 19 0
42 9 7 19 12 0
43 16 25 24 12 0
44 12 16 14 24 0
45 22 20 24 12 0
46 12 8 20 24 0
47 12 22 25 24 0
48 8 14 24 12 0
];

%% Variable Prescribed
% Node            Dimension                Value

dirichlet_data = [
19 1 0 
19 2 0 
19 3 0 
20 1 0 
20 2 0 
20 3 0 
21 1 0 
21 2 0 
21 3 0 
22 1 0 
22 2 0
22 3 0 
23 1 0 
23 2 0 
23 3 0 
24 1 0 
24 2 0 
24 3 0 
25 1 0 
25 2 0 
25 3 0 
26 1 0 
26 2 0 
26 3 0 
27 1 0 
27 2 0 
27 3 0 
];

%% Force Prescribed
% Node                Dimension                Value

pointload_complete = [
4 1 0 
4 2 -1 
];

isLeft   = @(coor) (abs(coor(:,1)) < 1e-12);
isRight  = @(coor) (abs(coor(:,1) - max(coor(:,1))) < 1e-12);
isMiddle = @(coor) (abs(coor(:,2) - max(coor(:,2)/2)) == 0);
isCenter = @(coor) (abs(coor(:,3) - max(coor(:,3)/2)) == 0);

% Dirichlet
sDir{1}.domain    = @(coor) isLeft(coor);
sDir{1}.direction = [1,2,3];
sDir{1}.value     = 0;

% Point load
sPL{1}.domain    = @(coor) isMiddle(coor) & isRight(coor) & isCenter(coor);
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
