%==================================================================
%                        General Data File
% Title: TETRAHEDRA
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

gidcoord = [
1           20            5            5
2           20            5          2.5
3           20          2.5            5
4           20          2.5          2.5
5           20            0            5
6           20            5            0
7           20            0          2.5
8           20          2.5            0
9       13.333            5            5
10           20            0            0
11       13.333          2.5            5
12       13.333            5          2.5
13       13.333          2.5          2.5
14       13.333            5            0
15       13.333            0            5
16       13.333            0          2.5
17       13.333          2.5            0
18       13.333            0            0
19       6.6667            5            5
20       6.6667          2.5            5
21       6.6667            5          2.5
22       6.6667          2.5          2.5
23       6.6667            5            0
24       6.6667            0            5
25       6.6667          2.5            0
26       6.6667            0          2.5
27       6.6667            0            0
28            0            5            5
29            0            5          2.5
30            0          2.5            5
31            0          2.5          2.5
32            0            0            5
33            0            5            0
34            0          2.5            0
35            0            0          2.5
36            0            0            0
];

%% Conectivities
% Element        Node(1)                Node(2)                Node(3)                Node(4)                Material

gidlnods = [
1 23 25 22 33 0
2 33 23 21 22 0
3 34 31 22 33 0
4 33 29 31 22 0
5 33 34 25 22 0
6 29 21 22 33 0
7 14 17 13 25 0
8 25 22 13 14 0
9 12 13 22 14 0
10 14 12 21 22 0
11 14 23 25 22 0
12 23 21 22 14 0
13 2 4 13 6 0
14 6 2 12 13 0
15 8 17 13 6 0
16 6 14 17 13 0
17 6 8 4 13 0
18 14 12 13 6 0
19 31 35 36 22 0
20 22 31 34 36 0
21 26 27 36 22 0
22 22 25 27 36 0
23 22 26 35 36 0
24 25 34 36 22 0
25 22 26 27 13 0
26 13 22 25 27 0
27 16 18 27 13 0
28 13 17 18 27 0
29 13 16 26 27 0
30 17 25 27 13 0
31 17 18 10 13 0
32 13 17 8 10 0
33 16 7 10 13 0
34 13 4 7 10 0
35 13 16 18 10 0
36 4 8 10 13 0
37 29 31 22 28 0
38 28 29 21 22 0
39 30 20 22 28 0
40 28 19 20 22 0
41 28 30 31 22 0
42 19 21 22 28 0
43 19 21 12 22 0
44 22 13 12 19 0
45 9 12 13 19 0
46 19 9 11 13 0
47 19 20 22 13 0
48 20 11 13 19 0
49 9 11 13 1 0
50 1 9 12 13 0
51 3 4 13 1 0
52 1 2 4 13 0
53 1 3 11 13 0
54 2 12 13 1 0
55 20 24 32 22 0
56 22 20 30 32 0
57 26 35 32 22 0
58 22 31 35 32 0
59 22 26 24 32 0
60 31 30 32 22 0
61 11 15 24 13 0
62 13 11 20 24 0
63 16 26 24 13 0
64 13 22 26 24 0
65 13 16 15 24 0
66 22 20 24 13 0
67 4 7 5 13 0
68 13 4 3 5 0
69 16 15 5 13 0
70 13 11 15 5 0
71 13 16 7 5 0
72 11 3 5 13 0
];

%% Variable Prescribed
% Node            Dimension                Value

lnodes = [
28 1 0 
28 2 0
28 3 0
32 1 0 
32 2 0
32 3 0
33 1 0 
33 2 0
33 3 0
36 1 0 
36 2 0
36 3 0
];

%% Force Prescribed
% Node                Dimension                Value

pointload_complete = [
1 1 2 
1 2 0
1 3 0
5 1 2 
5 2 0
5 3 0
6 1 2 
6 2 0
6 3 0
10 1 2 
10 2 0
10 3 0
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
