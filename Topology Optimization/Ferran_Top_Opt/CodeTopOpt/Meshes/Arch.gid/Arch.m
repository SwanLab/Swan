<<<<<<< HEAD
%==================================================================
%                        General Data File
% Title: Default_title
% Units: SI
% Dimensions: 2D
% Type of problem: Plane_Stress
% Type of Phisics: ELASTIC
% Micro/Macro: MACRO
%
%==================================================================

%% Data

Data_prb = {
'Default_title';
'SI';
'2D';
'Plane_Stress';
'ELASTIC';
'MACRO';
};

%% Coordinates
% Node                X                Y                Z

gidcoord = [
1            0            1            0
2         0.25         0.75            0
3            0          0.5            0
4          0.5            1            0
5          0.5          0.5            0
6         0.25         0.25            0
7         0.75         0.75            0
8            0            0            0
9            1            1            0
10         0.75         0.25            0
11            1          0.5            0
12          0.5            0            0
13         1.25         0.75            0
14            1            0            0
15         1.25         0.25            0
16          1.5            1            0
17          1.5          0.5            0
18         1.75         0.75            0
19          1.5            0            0
20         1.75         0.25            0
21            2            1            0
22            2          0.5            0
23            2            0            0
];

%% Conectivities
% Element        Node(1)                Node(2)                Node(3)                Material

gidlnods = [
1 3 2 1 0
2 5 2 3 0
3 4 2 5 0
4 1 2 4 0
5 8 6 3 0
6 12 6 8 0
7 5 6 12 0
8 3 6 5 0
9 5 7 4 0
10 11 7 5 0
11 9 7 11 0
12 4 7 9 0
13 12 10 5 0
14 14 10 12 0
15 11 10 14 0
16 5 10 11 0
17 22 20 23 0
18 17 20 22 0
19 19 20 17 0
20 23 20 19 0
21 21 18 22 0
22 16 18 21 0
23 17 18 16 0
24 22 18 17 0
25 17 15 19 0
26 11 15 17 0
27 14 15 11 0
28 19 15 14 0
29 16 13 17 0
30 9 13 16 0
31 11 13 9 0
32 17 13 11 0
];

%% Variable Prescribed
% Node            Dimension                Value

lnodes = [
8 1 0 
8 2 0 
23 1 0 
23 2 0 
];

%% Force Prescribed
% Node                Dimension                Value

pointload_complete = [
14 1 0 
14 2 -1 
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
=======
%==================================================================
%                        General Data File
% Title: Default_title
% Units: SI
% Dimensions: 2D
% Type of problem: Plane_Stress
% Type of Phisics: ELASTIC
% Micro/Macro: MACRO
%
%==================================================================

%% Data

Data_prb = {
'Default_title';
'SI';
'2D';
'Plane_Stress';
'ELASTIC';
'MACRO';
};

%% Coordinates
% Node                X                Y                Z

coord = [
1            0            1            0
2         0.25         0.75            0
3            0          0.5            0
4          0.5            1            0
5          0.5          0.5            0
6         0.25         0.25            0
7         0.75         0.75            0
8            0            0            0
9            1            1            0
10         0.75         0.25            0
11            1          0.5            0
12          0.5            0            0
13         1.25         0.75            0
14            1            0            0
15         1.25         0.25            0
16          1.5            1            0
17          1.5          0.5            0
18         1.75         0.75            0
19          1.5            0            0
20         1.75         0.25            0
21            2            1            0
22            2          0.5            0
23            2            0            0
];

%% Conectivities
% Element        Node(1)                Node(2)                Node(3)                Material

connec = [
1 3 2 1 0
2 5 2 3 0
3 4 2 5 0
4 1 2 4 0
5 8 6 3 0
6 12 6 8 0
7 5 6 12 0
8 3 6 5 0
9 5 7 4 0
10 11 7 5 0
11 9 7 11 0
12 4 7 9 0
13 12 10 5 0
14 14 10 12 0
15 11 10 14 0
16 5 10 11 0
17 22 20 23 0
18 17 20 22 0
19 19 20 17 0
20 23 20 19 0
21 21 18 22 0
22 16 18 21 0
23 17 18 16 0
24 22 18 17 0
25 17 15 19 0
26 11 15 17 0
27 14 15 11 0
28 19 15 14 0
29 16 13 17 0
30 9 13 16 0
31 11 13 9 0
32 17 13 11 0
];

%% Variable Prescribed
% Node            Dimension                Value

dirichlet_data = [
8 1 0 
8 2 0 
23 1 0 
23 2 0 
];

%% Force Prescribed
% Node                Dimension                Value

pointload_complete = [
14 1 0 
14 2 -1 
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
>>>>>>> refs/remotes/origin/master
