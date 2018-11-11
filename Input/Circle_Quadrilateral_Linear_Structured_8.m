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
2            0         1.75            0
3         0.25            2            0
4         0.25         1.75            0
5          0.5            2            0
6            0          1.5            0
7          0.5         1.75            0
8         0.25          1.5            0
9          0.5          1.5            0
10         0.75            2            0
11            0         1.25            0
12         0.75         1.75            0
13         0.25         1.25            0
14          0.5         1.25            0
15         0.75          1.5            0
16            1            2            0
17            0            1            0
18            1         1.75            0
19         0.25            1            0
20         0.75         1.25            0
21          0.5            1            0
22            1          1.5            0
23            1         1.25            0
24         0.75            1            0
25            0         0.75            0
26         1.25            2            0
27         1.25         1.75            0
28         0.25         0.75            0
29         1.25          1.5            0
30          0.5         0.75            0
31            1            1            0
32         1.25         1.25            0
33         0.75         0.75            0
34          1.5            2            0
35            0          0.5            0
36         0.25          0.5            0
37          1.5         1.75            0
38          1.5          1.5            0
39          0.5          0.5            0
40            1         0.75            0
41         1.25            1            0
42          1.5         1.25            0
43         0.75          0.5            0
44            0         0.25            0
45         1.75            2            0
46         0.25         0.25            0
47         1.75         1.75            0
48         1.25         0.75            0
49            1          0.5            0
50          1.5            1            0
51          0.5         0.25            0
52         1.75          1.5            0
53         1.75         1.25            0
54         0.75         0.25            0
55          1.5         0.75            0
56         1.25          0.5            0
57            2            2            0
58            0            0            0
59            2         1.75            0
60            1         0.25            0
61         0.25            0            0
62         1.75            1            0
63            2          1.5            0
64          0.5            0            0
65          1.5          0.5            0
66         0.75            0            0
67            2         1.25            0
68         1.75         0.75            0
69         1.25         0.25            0
70            2            1            0
71            1            0            0
72          1.5         0.25            0
73         1.75          0.5            0
74            2         0.75            0
75         1.25            0            0
76         1.75         0.25            0
77            2          0.5            0
78          1.5            0            0
79            2         0.25            0
80         1.75            0            0
81            2            0            0
];

%% Conectivities
% Element        Node(1)                Node(2)                Node(3)                Node(4)                Material

gidlnods = [
1 61 46 44 58 0
2 64 51 46 61 0
3 66 54 51 64 0
4 71 60 54 66 0
5 75 69 60 71 0
6 78 72 69 75 0
7 80 76 72 78 0
8 81 79 76 80 0
9 46 36 35 44 0
10 51 39 36 46 0
11 54 43 39 51 0
12 60 49 43 54 0
13 69 56 49 60 0
14 72 65 56 69 0
15 76 73 65 72 0
16 79 77 73 76 0
17 36 28 25 35 0
18 39 30 28 36 0
19 43 33 30 39 0
20 49 40 33 43 0
21 56 48 40 49 0
22 65 55 48 56 0
23 73 68 55 65 0
24 77 74 68 73 0
25 28 19 17 25 0
26 30 21 19 28 0
27 33 24 21 30 0
28 40 31 24 33 0
29 48 41 31 40 0
30 55 50 41 48 0
31 68 62 50 55 0
32 74 70 62 68 0
33 19 13 11 17 0
34 21 14 13 19 0
35 24 20 14 21 0
36 31 23 20 24 0
37 41 32 23 31 0
38 50 42 32 41 0
39 62 53 42 50 0
40 70 67 53 62 0
41 13 8 6 11 0
42 14 9 8 13 0
43 20 15 9 14 0
44 23 22 15 20 0
45 32 29 22 23 0
46 42 38 29 32 0
47 53 52 38 42 0
48 67 63 52 53 0
49 8 4 2 6 0
50 9 7 4 8 0
51 15 12 7 9 0
52 22 18 12 15 0
53 29 27 18 22 0
54 38 37 27 29 0
55 52 47 37 38 0
56 63 59 47 52 0
57 4 3 1 2 0
58 7 5 3 4 0
59 12 10 5 7 0
60 18 16 10 12 0
61 27 26 16 18 0
62 37 34 26 27 0
63 47 45 34 37 0
64 59 57 45 47 0
];

%% Variable Prescribed
% Node            Dimension                Value

lnodes = [
58 1 0 
58 2 0 
81 1 0 
81 2 0 
];

%% Force Prescribed
% Node                Dimension                Value

pointload_complete = [
1 1 0 
1 2 1 
57 1 0 
57 2 1 
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
