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

gidcoord = [
1            0            2            0
2            0          1.5            0
3          0.5            2            0
4          0.5          1.5            0
5            0            1            0
6            1            2            0
7            1          1.5            0
8          0.5            1            0
9            1            1            0
10          1.5            2            0
11            0          0.5            0
12          1.5          1.5            0
13          0.5          0.5            0
14          1.5            1            0
15            1          0.5            0
16            2            2            0
17            0            0            0
18          0.5            0            0
19            2          1.5            0
20          1.5          0.5            0
21            2            1            0
22            1            0            0
23            2          0.5            0
24          1.5            0            0
25            2            0            0
];

%% Conectivities
% Element        Node(1)                Node(2)                Node(3)                Node(4)                Material

gidlnods = [
1 18 13 11 17 0
2 22 15 13 18 0
3 24 20 15 22 0
4 25 23 20 24 0
5 13 8 5 11 0
6 15 9 8 13 0
7 20 14 9 15 0
8 23 21 14 20 0
9 8 4 2 5 0
10 9 7 4 8 0
11 14 12 7 9 0
12 21 19 12 14 0
13 4 3 1 2 0
14 7 6 3 4 0
15 12 10 6 7 0
16 19 16 10 12 0
];

%% Variable Prescribed
% Node            Dimension                Value

lnodes = [
17 1 0 
17 2 0 
25 1 0 
25 2 0 
];

%% Force Prescribed
% Node                Dimension                Value

pointload_complete = [
1 1 0 
1 2 1 
16 1 0 
16 2 1 
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
