%==================================================================
%                        General Data File
% Title: TRIANGLE
% Units: SI
% Dimensions: 2D
% Type of problem: Plane_Stress
% Type of Phisics: ELASTIC
% Micro/Macro: MACRO
%
%==================================================================

%% Data

Data_prb = {
'TRIANGLE';
'SI';
'2D';
'Plane_Stress';
'ELASTIC';
'MACRO';
};

%% Coordinates
% Node                X                Y                Z

coord = [
1            2            0            0
2            2         -0.2            0
3            2         -0.4            0
4          1.5            0            0
5          1.5         -0.2            0
6          1.5         -0.4            0
7            1            0            0
8            1         -0.2            0
9            1         -0.4            0
10          0.5            0            0
11          0.5         -0.2            0
12          0.5         -0.4            0
13            0            0            0
14            0         -0.2            0
15            0         -0.4            0
];

%% Conectivities
% Element        Node(1)                Node(2)                Node(3)                Material

connec = [
1 14 11 13 0
2 11 10 13 0
3 15 12 14 0
4 12 11 14 0
5 11 8 10 0
6 8 7 10 0
7 12 9 11 0
8 9 8 11 0
9 8 5 7 0
10 5 4 7 0
11 9 6 8 0
12 6 5 8 0
13 5 2 4 0
14 2 1 4 0
15 6 3 5 0
16 3 2 5 0
];

%% Variable Prescribed
% Node            Dimension                Value

dirichlet_data = [
   13 1 0
   13 2 0
   14 1 0
   14 2 0
   15 1 0
   15 2 0
];

%% Force Prescribed
% Node                Dimension                Value

pointload_complete = [
    2 2 -1
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
