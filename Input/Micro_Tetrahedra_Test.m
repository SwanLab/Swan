%==================================================================
%                        General Data File
% Title: TETRAHEDRA
% Units: SI
% Dimensions: 3D
% Type of problem: Plane_Stress
% Type of Phisics: ELASTIC
% Micro/Macro: MICRO
%
%==================================================================

%% Data

Data_prb = {
'TETRAHEDRA';
'SI';
'3D';
'Plane_Stress';
'ELASTIC';
'MICRO';
};

%% Coordinates
% Node                X                Y                Z

gidcoord = [
1            0            1            1
2            0            1          0.5
3          0.5            1            1
4            0          0.5            1
5            0          0.5          0.5
6          0.5            1          0.5
7          0.5          0.5            1
8          0.5          0.5          0.5
9            0            0            1
10            1            1            1
11            0            1            0
12          0.5            1            0
13            1            1          0.5
14            1          0.5            1
15          0.5            0            1
16            0            0          0.5
17            0          0.5            0
18          0.5          0.5            0
19            1          0.5          0.5
20          0.5            0          0.5
21            1            0            1
22            1            1            0
23            0            0            0
24            1            0          0.5
25            1          0.5            0
26          0.5            0            0
27            1            0            0
];

%% Conectivities
% Element        Node(1)                Node(2)                Node(3)                Node(4)                Material

gidlnods = [
1 16 5 8 23 0
2 23 16 20 8 0
3 17 18 8 23 0
4 23 26 18 8 0
5 23 17 5 8 0
6 26 20 8 23 0
7 15 7 8 9 0
8 9 15 20 8 0
9 4 5 8 9 0
10 9 16 5 8 0
11 9 4 7 8 0
12 16 20 8 9 0
13 18 12 11 8 0
14 8 18 17 11 0
15 6 2 11 8 0
16 8 5 2 11 0
17 8 6 12 11 0
18 5 17 11 8 0
19 5 2 1 8 0
20 8 5 4 1 0
21 6 3 1 8 0
22 8 7 3 1 0
23 8 6 2 1 0
24 7 4 1 8 0
25 26 18 8 27 0
26 27 26 20 8 0
27 25 19 8 27 0
28 27 24 19 8 0
29 27 25 18 8 0
30 24 20 8 27 0
31 24 19 8 21 0
32 21 24 20 8 0
33 14 7 8 21 0
34 21 15 7 8 0
35 21 14 19 8 0
36 15 20 8 21 0
37 19 13 22 8 0
38 8 19 25 22 0
39 6 12 22 8 0
40 8 18 12 22 0
41 8 6 13 22 0
42 18 25 22 8 0
43 7 3 10 8 0
44 8 7 14 10 0
45 6 13 10 8 0
46 8 19 13 10 0
47 8 6 3 10 0
48 19 14 10 8 0
];

%% Variable Prescribed
% Node            Dimension                Value

lnodes = [
1 1 0 
1 2 0
1 3 0
9 1 0 
9 2 0
9 3 0
10 1 0 
10 2 0
10 3 0
11 1 0 
11 2 0
11 3 0
21 1 0 
21 2 0
21 3 0
22 1 0 
22 2 0
22 3 0
23 1 0 
23 2 0
23 3 0
27 1 0 
27 2 0
27 3 0
];

%% Force Prescribed
% Node                Dimension                Value

pointload_complete = [
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
Master_slave=[
   5    19
    20     6
    18     7
    26    15
    26     3
    26    12
    17     4
    17    14
    17    25
    16     2
    16    13
    16    24
];


%% Nodes solid
% Nodes that must remain 

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
