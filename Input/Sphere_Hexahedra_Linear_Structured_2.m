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
1            0            2            2
2            0            1            2
3            1            2            2
4            0            2            1
5            0            1            1
6            1            1            2
7            1            2            1
8            1            1            1
9            0            0            2
10            0            2            0
11            2            2            2
12            0            0            1
13            2            2            1
14            2            1            2
15            1            0            2
16            0            1            0
17            1            2            0
18            2            1            1
19            1            0            1
20            1            1            0
21            2            0            2
22            2            2            0
23            0            0            0
24            2            0            1
25            1            0            0
26            2            1            0
27            2            0            0
];

%% Conectivities
% Element        Node(1)                Node(2)                Node(3)                Node(4)                Node(5)                Node(6)                Node(7)                Node(8)                Material

connec = [
1 8 18 14 6 19 24 21 15 0
2 5 8 6 2 12 19 15 9 0
3 20 26 18 8 25 27 24 19 0
4 16 20 8 5 23 25 19 12 0
5 7 13 11 3 8 18 14 6 0
6 4 7 3 1 5 8 6 2 0
7 17 22 13 7 20 26 18 8 0
8 10 17 7 4 16 20 8 5 0
];

%% Variable Prescribed
% Node            Dimension                Value

dirichlet_data = [
1 1 0 
1 2 0
1 3 0
9 1 0 
9 2 0
9 3 0
11 1 0 
11 2 0
11 3 0
21 1 0 
21 2 0
21 3 0
];

%% Force Prescribed
% Node                Dimension                Value

pointload_complete = [
10 1 0 
10 2 0
10 3 -1
22 1 0 
22 2 0
22 3 -1
23 1 0 
23 2 0
23 3 -1
27 1 0 
27 2 0
27 3 -1
];

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
