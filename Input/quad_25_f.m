%==================================================================
%                        General Data File
% Title: Default_title
% Units: SI
% Dimensions: 2D
% Type of problem: Plane_Stress
% Type of Phisics: HYPERELASTIC
% Micro/Macro: MACRO
%
%==================================================================

%% Data

Data_prb = {
'QUAD';
'SI';
'2D';
'Plane_Stress';
'HYPERELASTIC';
'MACRO';
};

%% Coordinates
% Node                X                Y                Z

coord = [
1            1            1            0
2          0.8            1            0
3            1          0.8            0
4      0.79908      0.79932            0
5          0.6            1            0
6            1          0.6            0
7      0.59864      0.79845            0
8      0.79766      0.59882            0
9      0.59886      0.59853            0
10          0.4            1            0
11            1          0.4            0
12      0.79734      0.39957            0
13      0.39909      0.79766            0
14      0.40063      0.59593            0
15      0.59951      0.39726            0
16            1          0.2            0
17          0.2            1            0
18      0.79817       0.2005            0
19      0.20051       0.7971            0
20      0.40163      0.39679            0
21      0.59885       0.2003            0
22      0.20081      0.59657            0
23            1            0            0
24            0            1            0
25      0.20081      0.39764            0
26      0.40032      0.19858            0
27          0.8            0            0
28            0          0.8            0
29            0          0.6            0
30          0.6            0            0
31       0.2008      0.19722            0
32            0          0.4            0
33          0.4            0            0
34          0.2            0            0
35            0          0.2            0
36            0            0            0
];

%% Conectivities
% Element        Node(1)                Node(2)                Node(3)                Node(4)                Material

connec = [
1 35 36 34 31 0
2 27 23 16 18 0
3 3 1 2 4 0
4 17 24 28 19 0
5 32 35 31 25 0
6 30 27 18 21 0
7 6 3 4 8 0
8 10 17 19 13 0
9 29 32 25 22 0
10 28 29 22 19 0
11 34 33 26 31 0
12 33 30 21 26 0
13 16 11 12 18 0
14 11 6 8 12 0
15 2 5 7 4 0
16 5 10 13 7 0
17 25 31 26 20 0
18 4 7 9 8 0
19 18 12 15 21 0
20 12 8 9 15 0
21 26 21 15 20 0
22 13 19 22 14 0
23 7 13 14 9 0
24 15 9 14 20 0
25 14 22 25 20 0
];

%% Variable Prescribed
% Node            Dimension                Value

dirichlet_data = [
24 1 0 
24 2 0 
28 1 0 
28 2 0 
29 1 0 
29 2 0 
32 1 0 
32 2 0 
35 1 0 
35 2 0 
36 1 0 
36 2 0 
];

%% Force Prescribed
% Node                Dimension                Value

pointload_complete = [
1 2 1 
3 2 1 
6 2 1 
11 2 1 
16 2 1 
23 2 1 
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
