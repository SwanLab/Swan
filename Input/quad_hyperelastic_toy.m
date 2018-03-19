%==================================================================
%                        General Data File
% Title: QUADRILAT
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
2            1          0.5            0
3          0.5            1            0
4          0.5          0.5            0
5            0            1            0
6            1            0            0
7            0          0.5            0
8          0.5            0            0
9            0            0            0
];

%% Conectivities
% Element        Node(1)                Node(2)                Node(3)                Node(4)                Material

connec = [
1 7 9 8 4 0
2 5 7 4 3 0
3 1 3 4 2 0
4 2 4 8 6 0
];

%% Variable Prescribed
% Node            Dimension                Value

dirichlet_data = [
5 1 0 
5 2 0 
7 1 0 
7 2 0 
9 1 0 
9 2 0 
];

%% Force Prescribed
% Node                Dimension                Value

pointload_complete = [
1 1 1 
2 1 1 
6 1 1 
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
