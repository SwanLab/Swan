%==================================================================
%                        General Data File
% Title: tetra_y
% Units: SI
% Dimensions: 3D
% Type of problem: Plane_Stress
% Type of Phisics: HYPERELASTIC
% Micro/Macro: MACRO
%
%==================================================================

%% Data

Data_prb = {
'TETRAHEDRA';
'SI';
'3D';
'Plane_Stress';
'HYPERELASTIC';
'MACRO';
};

%% Coordinates
% Node                X                Y                Z

coord = [
1            0            0            0
2          0.5            0            0
3            0            0          0.5
4            0          0.5            0
5      0.30843       0.3144      0.30672
6      0.49495            0      0.50505
7      0.49495      0.50505            0
8            0      0.50505      0.50505
9      0.63092      0.32647      0.34224
10       0.3501      0.63394      0.32301
11      0.38579      0.41013      0.65156
12            0            0            1
13            1            0            0
14            0            1            0
15      0.72289      0.76028      0.25085
16      0.72171       0.2879      0.75617
17      0.27032      0.72723      0.75817
18            1            0          0.5
19            0            1          0.5
20          0.5            1            0
21            0          0.5            1
22            1          0.5            0
23          0.5            0            1
24      0.66835      0.66835       0.6484
25      0.49495            1      0.50505
26            1      0.50505      0.49495
27      0.50505      0.50505            1
28            0            1            1
29            1            1            0
30            1            0            1
31          0.5            1            1
32            1          0.5            1
33            1            1          0.5
34            1            1            1
];

%% Conectivities
% Element        Node(1)                Node(2)                Node(3)                Node(4)                Material

connec = [
1 32 27 31 24 0
2 32 27 24 16 0
3 27 24 16 11 0
4 16 27 11 23 0
5 24 16 11 9 0
6 16 11 9 6 0
7 11 24 9 10 0
8 9 11 10 5 0
9 10 9 5 7 0
10 9 11 5 6 0
11 11 10 5 8 0
12 24 9 10 15 0
13 9 10 15 7 0
14 10 24 15 25 0
15 24 9 15 26 0
16 9 15 26 22 0
17 11 24 10 17 0
18 10 11 17 8 0
19 24 10 17 25 0
20 10 17 25 19 0
21 11 24 17 27 0
22 17 24 25 31 0
23 15 24 26 33 0
24 24 16 9 26 0
25 17 11 27 21 0
26 15 10 25 20 0
27 16 9 26 18 0
28 5 10 7 4 0
29 5 11 8 3 0
30 5 9 6 2 0
31 24 15 25 33 0
32 24 17 27 31 0
33 9 5 7 2 0
34 11 5 6 3 0
35 10 15 7 20 0
36 10 5 8 4 0
37 24 16 26 32 0
38 9 16 6 18 0
39 11 17 8 21 0
40 25 17 31 19 0
41 16 11 6 23 0
42 15 9 7 22 0
43 17 10 8 19 0
44 26 15 33 22 0
45 32 27 16 23 0
46 17 27 31 21 0
47 15 25 33 20 0
48 16 26 32 18 0
49 7 5 4 2 0
50 6 16 23 18 0
51 24 26 33 32 0
52 5 8 4 3 0
53 5 6 3 2 0
54 7 15 22 20 0
55 25 24 33 31 0
56 10 25 20 19 0
57 8 17 19 21 0
58 9 26 18 22 0
59 11 27 21 23 0
60 7 10 20 4 0
61 8 11 21 3 0
62 6 9 18 2 0
63 11 6 23 3 0
64 9 7 22 2 0
65 10 8 19 4 0
66 15 33 22 29 0
67 15 33 29 20 0
68 22 15 29 20 0
69 17 31 19 28 0
70 17 31 28 21 0
71 19 17 28 21 0
72 16 32 23 30 0
73 16 32 30 18 0
74 23 16 30 18 0
75 5 4 2 1 0
76 5 4 1 3 0
77 2 5 1 3 0
78 34 32 31 24 0
79 34 32 24 33 0
80 24 34 33 31 0
81 13 22 2 9 0
82 13 22 9 18 0
83 9 13 18 2 0
84 12 21 23 3 0
85 14 4 20 10 0
86 20 14 10 19 0
87 14 10 19 4 0
88 21 11 23 3 0
];

%% Variable Prescribed
% Node            Dimension                Value

dirichlet_data = [
1 1 0 
1 2 0
1 3 0
3 1 0 
3 2 0
3 3 0
4 1 0 
4 2 0
4 3 0
8 1 0 
8 2 0
8 3 0
12 1 0 
12 2 0
12 3 0
14 1 0 
14 2 0
14 3 0
19 1 0 
19 2 0
19 3 0
21 1 0 
21 2 0
21 3 0
28 1 0 
28 2 0
28 3 0
];

%% Force Prescribed
% Node                Dimension                Value

pointload_complete = [
13 2 -1 
18 2 -1 
22 2 -1 
26 2 -1 
29 2 -1 
30 2 -1 
32 2 -1 
33 2 -1 
34 2 -1 
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
