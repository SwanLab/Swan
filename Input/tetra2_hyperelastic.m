%==================================================================
%                        General Data File
% Title: quad_3D
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
1            1            1            1
2          0.5            1            1
3            1          0.5            1
4            1            1          0.5
5          0.5          0.5            1
6          0.5            1          0.5
7            1          0.5          0.5
8          0.5          0.5          0.5
9            1            1            0
10            1            0            1
11            0            1            1
12            1          0.5            0
13            0            1          0.5
14          0.5            1            0
15            1            0          0.5
16          0.5            0            1
17            0          0.5            1
18            0          0.5          0.5
19          0.5          0.5            0
20          0.5            0          0.5
21            0            1            0
22            1            0            0
23            0            0            1
24          0.5            0            0
25            0            0          0.5
26            0          0.5            0
27            0            0            0
];

%% Conectivities
% Element        Node(1)                Node(2)                Node(3)                Node(4)                Material

connec = [
1 2 6 8 11 0
2 11 2 5 8 0
3 13 18 8 11 0
4 11 17 18 8 0
5 11 13 6 8 0
6 17 5 8 11 0
7 3 7 8 1 0
8 1 3 5 8 0
9 4 6 8 1 0
10 1 2 6 8 0
11 1 4 7 8 0
12 2 5 8 1 0
13 18 26 21 8 0
14 8 18 13 21 0
15 19 14 21 8 0
16 8 6 14 21 0
17 8 19 26 21 0
18 6 13 21 8 0
19 6 14 9 8 0
20 8 6 4 9 0
21 19 12 9 8 0
22 8 7 12 9 0
23 8 19 14 9 0
24 7 4 9 8 0
25 17 18 8 23 0
26 23 17 5 8 0
27 25 20 8 23 0
28 23 16 20 8 0
29 23 25 18 8 0
30 16 5 8 23 0
31 16 20 8 10 0
32 10 16 5 8 0
33 15 7 8 10 0
34 10 3 7 8 0
35 10 15 20 8 0
36 3 5 8 10 0
37 20 24 27 8 0
38 8 20 25 27 0
39 19 26 27 8 0
40 8 18 26 27 0
41 8 19 24 27 0
42 18 25 27 8 0
43 7 12 22 8 0
44 8 7 15 22 0
45 19 24 22 8 0
46 8 20 24 22 0
47 8 19 12 22 0
48 20 15 22 8 0
];

%% Variable Prescribed
% Node            Dimension                Value

dirichlet_data = [
11 1 0 
11 2 0
11 3 0
13 1 0 
13 2 0
13 3 0
17 1 0 
17 2 0
17 3 0
18 1 0 
18 2 0
18 3 0
21 1 0 
21 2 0
21 3 0
23 1 0 
23 2 0
23 3 0
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
1 1 1 
3 1 1 
4 1 1 
7 1 1 
9 1 1 
10 1 1 
12 1 1 
15 1 1 
22 1 1 
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
