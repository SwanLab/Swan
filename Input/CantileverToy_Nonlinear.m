%==================================================================
%                        General Data File
% Title: TRIANGLE
% Units: SI
% Dimensions: 2D
% Type of problem: Plane_Stress
% Type of Phisics: HYPERELASTIC
% Micro/Macro: MACRO
%
%==================================================================

%% Data

Data_prb = {
'TRIANGLE';
'SI';
'2D';
'Plane_Stress';
'HYPERELASTIC';
'MACRO';
};

%% Coordinates
% Node                X                Y                Z

gidcoord = [
1            0            0            0
2            1            0            0
3            0            1            0
4            1            1            0
];

%% Conectivities
% Element        Node(1)                Node(2)                Node(3)                Material

gidlnods = [
1 1 2 3 0
2 2 4 3 0
];

%% Variable Prescribed
% Node            Dimension                Value

lnodes = [
1 1 0 
1 2 0 
3 1 0 
3 2 0
];

%% Force Prescribed
% Node                Dimension                Value

pointload_complete = [
2 1 1
4 1 1
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
