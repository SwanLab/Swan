%==================================================================
%                        General Data File
% Title: HEXAHEDRA
% Units: SI
% Dimensions: 3D
% Type of problem: Plane_Stress
% Type of Phisics: ELASTIC
% Micro/Macro: MACRO
%
%==================================================================

%% Data

Data_prb = {
'HEXAHEDRA';
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
3          1.5            0            1
4            2          0.5            1
5            2          0.5          0.5
6          1.5          0.5            1
7          1.5            0          0.5
8          1.5          0.5          0.5
9            1            0            1
10            2            0            0
11            2            1            1
12            1            0          0.5
13            2            1          0.5
14          1.5            1            1
15            1          0.5            1
16            2          0.5            0
17          1.5            0            0
18          1.5            1          0.5
19            1          0.5          0.5
20          1.5          0.5            0
21            1            1            1
22            2            1            0
23            1            0            0
24            1            1          0.5
25          1.5            1            0
26            1          0.5            0
27          0.5            0            1
28          0.5            0          0.5
29          0.5          0.5            1
30          0.5          0.5          0.5
31            1            1            0
32          0.5            0            0
33          0.5            1            1
34          0.5            1          0.5
35          0.5          0.5            0
36            0            0            1
37            0          0.5            1
38            0            0          0.5
39          0.5            1            0
40            0          0.5          0.5
41            0            1            1
42            0            0            0
43            0            1          0.5
44            0          0.5            0
45            0            1            0
];

%% Conectivities
% Element        Node(1)                Node(2)                Node(3)                Node(4)                Node(5)                Node(6)                Node(7)                Node(8)                Material

connec = [
1 30 40 44 35 28 38 42 32 0
2 19 30 35 26 12 28 32 23 0
3 8 19 26 20 7 12 23 17 0
4 5 8 20 16 2 7 17 10 0
5 29 37 40 30 27 36 38 28 0
6 15 29 30 19 9 27 28 12 0
7 6 15 19 8 3 9 12 7 0
8 4 6 8 5 1 3 7 2 0
9 34 43 45 39 30 40 44 35 0
10 24 34 39 31 19 30 35 26 0
11 18 24 31 25 8 19 26 20 0
12 13 18 25 22 5 8 20 16 0
13 33 41 43 34 29 37 40 30 0
14 21 33 34 24 15 29 30 19 0
15 14 21 24 18 6 15 19 8 0
16 11 14 18 13 4 6 8 5 0
];

%% Variable Prescribed
% Node            Dimension                Value

dirichlet_data = [
36 1 0 
36 2 0
36 3 0
37 1 0 
37 2 0
37 3 0
38 1 0 
38 2 0
38 3 0
40 1 0 
40 2 0
40 3 0
41 1 0 
41 2 0
41 3 0
42 1 0 
42 2 0
42 3 0
43 1 0 
43 2 0
43 3 0
44 1 0 
44 2 0
44 3 0
45 1 0 
45 2 0
45 3 0
];

%% Force Prescribed
% Node                Dimension                Value

pointload_complete = [
5 1 0 
5 2 -1
5 3 0
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
