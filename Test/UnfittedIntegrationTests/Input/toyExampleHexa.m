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
2      0.33333            0            0
3      0.66667            0            0
4            1            0            0
5            0            1            0
6            0            0            1
7      0.33333            1            0
8      0.33333            0            1
9      0.66667            0            1
10      0.66667            1            0
11            1            1            0
12            1            0            1
13            0            1            1
14      0.33333            1            1
15      0.66667            1            1
16            1            1            1
];

%% Conectivities
% Element        Node(1)                Node(2)                Node(3)                Node(4)                Node(5)                Node(6)                Node(7)                Node(8)                Material

gidlnods = [
1 14 13 5 7 8 6 1 2 0
2 15 14 7 10 9 8 2 3 0
3 16 15 10 11 12 9 3 4 0
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
