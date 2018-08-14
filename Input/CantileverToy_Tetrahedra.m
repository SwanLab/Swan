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

coord = [
1           12            0            2
2           12            1            1
3           12           -1            1
4           12            0            0
5          9.6            0            2
6          9.6           -1            1
7          9.6            1            1
8           12            2            0
9           12           -2            0
10          9.6            0            0
11           12           -1           -1
12           12            1           -1
13          9.6            2            0
14          9.6           -2            0
15          9.6           -1           -1
16          9.6            1           -1
17           12            0           -2
18          9.6            0           -2
19          7.2            0            2
20          7.2           -1            1
21          7.2            1            1
22          7.2            0            0
23          7.2           -2            0
24          7.2            2            0
25          7.2           -1           -1
26          7.2            1           -1
27          7.2            0           -2
28          4.8            0            2
29          4.8            1            1
30          4.8           -1            1
31          4.8            0            0
32          4.8            2            0
33          4.8           -2            0
34          4.8           -1           -1
35          4.8            1           -1
36          4.8            0           -2
37          2.4            0            2
38          2.4           -1            1
39          2.4            1            1
40          2.4            0            0
41          2.4           -2            0
42          2.4            2            0
43          2.4           -1           -1
44          2.4            1           -1
45          2.4            0           -2
46            0            0            2
47            0           -1            1
48            0            1            1
49            0            0            0
50            0           -2            0
51            0            2            0
52            0            1           -1
53            0           -1           -1
54            0            0           -2
];

%% Conectivities
% Element        Node(1)                Node(2)                Node(3)                Node(4)                Material

connec = [
1 53 43 40 50 0
2 50 53 49 40 0
3 41 38 40 50 0
4 50 47 38 40 0
5 50 41 43 40 0
6 47 49 40 50 0
7 52 44 40 54 0
8 54 52 49 40 0
9 45 43 40 54 0
10 54 53 43 40 0
11 54 45 44 40 0
12 53 49 40 54 0
13 43 34 31 41 0
14 41 43 40 31 0
15 33 30 31 41 0
16 41 38 30 31 0
17 41 33 34 31 0
18 38 40 31 41 0
19 44 35 31 45 0
20 45 44 40 31 0
21 36 34 31 45 0
22 45 43 34 31 0
23 45 36 35 31 0
24 43 40 31 45 0
25 34 25 22 33 0
26 33 34 31 22 0
27 23 20 22 33 0
28 33 30 20 22 0
29 33 23 25 22 0
30 30 31 22 33 0
31 35 26 22 36 0
32 36 35 31 22 0
33 27 25 22 36 0
34 36 34 25 22 0
35 36 27 26 22 0
36 34 31 22 36 0
37 25 15 10 23 0
38 23 25 22 10 0
39 14 6 10 23 0
40 23 20 6 10 0
41 23 14 15 10 0
42 20 22 10 23 0
43 26 16 10 27 0
44 27 26 22 10 0
45 18 15 10 27 0
46 27 25 15 10 0
47 27 18 16 10 0
48 25 22 10 27 0
49 6 3 9 10 0
50 10 6 14 9 0
51 4 11 9 10 0
52 10 15 11 9 0
53 10 4 3 9 0
54 15 14 9 10 0
55 15 11 17 10 0
56 10 15 18 17 0
57 4 12 17 10 0
58 10 16 12 17 0
59 10 4 11 17 0
60 16 18 17 10 0
61 47 38 40 46 0
62 46 47 49 40 0
63 37 39 40 46 0
64 46 48 39 40 0
65 46 37 38 40 0
66 48 49 40 46 0
67 48 39 40 51 0
68 51 48 49 40 0
69 42 44 40 51 0
70 51 52 44 40 0
71 51 42 39 40 0
72 52 49 40 51 0
73 38 30 31 37 0
74 37 38 40 31 0
75 28 29 31 37 0
76 37 39 29 31 0
77 37 28 30 31 0
78 39 40 31 37 0
79 39 29 31 42 0
80 42 39 40 31 0
81 32 35 31 42 0
82 42 44 35 31 0
83 42 32 29 31 0
84 44 40 31 42 0
85 30 20 22 28 0
86 28 30 31 22 0
87 19 21 22 28 0
88 28 29 21 22 0
89 28 19 20 22 0
90 29 31 22 28 0
91 29 21 22 32 0
92 32 29 31 22 0
93 24 26 22 32 0
94 32 35 26 22 0
95 32 24 21 22 0
96 35 31 22 32 0
97 20 6 10 19 0
98 19 20 22 10 0
99 5 7 10 19 0
100 19 21 7 10 0
101 19 5 6 10 0
102 21 22 10 19 0
103 21 7 10 24 0
104 24 21 22 10 0
105 13 16 10 24 0
106 24 26 16 10 0
107 24 13 7 10 0
108 26 22 10 24 0
109 7 2 1 10 0
110 10 7 5 1 0
111 4 3 1 10 0
112 10 6 3 1 0
113 10 4 2 1 0
114 6 5 1 10 0
115 16 12 8 10 0
116 10 16 13 8 0
117 4 2 8 10 0
118 10 7 2 8 0
119 10 4 12 8 0
120 7 13 8 10 0
];

%% Variable Prescribed
% Node            Dimension                Value

dirichlet_data = [
46 1 0 
46 2 0
46 3 0
47 1 0 
47 2 0
47 3 0
48 1 0 
48 2 0
48 3 0
49 1 0 
49 2 0
49 3 0
50 1 0 
50 2 0
50 3 0
51 1 0 
51 2 0
51 3 0
52 1 0 
52 2 0
52 3 0
53 1 0 
53 2 0
53 3 0
54 1 0 
54 2 0
54 3 0
];

%% Force Prescribed
% Node                Dimension                Value

pointload_complete = [
4 1 0 
4 2 -1
4 3 0
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
