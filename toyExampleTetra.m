%==================================================================
%                        General Data File
% Title: Default_title
% Units: SI
% Dimensions: 3D
% Type of problem: Plane_Stress
% Type of Phisics: ELASTIC
% Micro/Macro: MACRO
%
%==================================================================

%% Data

Data_prb = {
'Default_title';
'SI';
'3D';
'Plane_Stress';
'ELASTIC';
'MACRO';
};

%% Coordinates
% Node                X                Y                Z

gidcoord = [
1            0            0            0
2          0.5            0            0
3            0          0.5            0
4            0            0          0.5
5      0.35569      0.33625      0.34788
6          0.5            0          0.5
7          0.5          0.5            0
8            0          0.5          0.5
9      0.71695      0.28611      0.29082
10      0.48605      0.64783      0.34894
11      0.65624      0.34802      0.64402
12      0.30434      0.63189      0.68914
13            1            0            0
14            0            1            0
15            0            0            1
16            1          0.5            0
17          0.5            1            0
18          0.5            0            1
19            1            0          0.5
20            0            1          0.5
21            0          0.5            1
22            1          0.5          0.5
23          0.5            1          0.5
24          0.5          0.5            1
25        0.688      0.76418      0.68791
26            1            1            0
27            1            0            1
28            0            1            1
29            1          0.5            1
30            1            1          0.5
31          0.5            1            1
32            1            1            1
];

%% Conectivities
% Element        Node(1)                Node(2)                Node(3)                Node(4)                Material

gidlnods = [
1 1 2 3 4 0
2 13 19 16 2 0
3 14 20 3 17 0
4 26 30 17 16 0
5 15 21 18 4 0
6 3 2 7 5 0
7 2 7 5 9 0
8 7 5 9 10 0
9 7 5 10 3 0
10 9 7 10 16 0
11 9 7 16 2 0
12 5 9 10 11 0
13 10 5 11 12 0
14 10 5 12 8 0
15 5 12 8 4 0
16 11 10 12 25 0
17 10 12 25 23 0
18 25 10 23 30 0
19 10 12 23 20 0
20 10 23 30 17 0
21 12 25 23 31 0
22 25 23 31 30 0
23 12 25 31 24 0
24 25 31 24 29 0
25 12 25 24 11 0
26 24 12 11 18 0
27 25 24 11 29 0
28 24 11 29 18 0
29 10 23 17 20 0
30 11 10 25 22 0
31 11 10 22 9 0
32 22 11 9 19 0
33 10 22 9 16 0
34 22 9 16 19 0
35 5 11 12 18 0
36 5 9 11 6 0
37 11 5 6 18 0
38 5 9 6 2 0
39 9 11 6 19 0
40 6 9 19 2 0
41 8 5 4 3 0
42 11 25 29 22 0
43 31 12 24 21 0
44 12 24 21 18 0
45 29 11 22 19 0
46 6 5 2 4 0
47 6 5 4 18 0
48 8 5 3 10 0
49 6 11 18 19 0
50 10 25 22 30 0
51 22 10 30 16 0
52 12 23 20 31 0
53 25 29 22 30 0
54 12 10 8 20 0
55 8 12 20 21 0
56 8 12 21 4 0
57 25 31 29 32 0
58 25 31 32 30 0
59 29 25 32 30 0
60 7 10 16 17 0
61 7 10 17 3 0
62 9 16 19 2 0
63 12 20 21 28 0
64 21 12 28 31 0
65 12 20 28 31 0
66 11 29 18 27 0
67 11 29 27 19 0
68 18 11 27 19 0
69 5 12 4 18 0
70 3 20 8 10 0
71 5 4 3 2 0
72 30 10 17 16 0
73 12 21 4 18 0
74 17 10 20 3 0
];

%% Variable Prescribed
% Node            Dimension                Value

lnodes = [
];

%% Force Prescribed
% Node                Dimension                Value

pointload_complete = [
];

pointload_adjoint = [
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

%nodesolid = unique(pointload_complete(:,1));

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
