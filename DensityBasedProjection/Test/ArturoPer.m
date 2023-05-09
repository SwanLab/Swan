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
'QUADS';
'SI';
'2D';
'Plane_Stress';
'ELASTIC';
'MACRO';
};

%% Coordinates
% Node                X                Y                Z
load('coordsMatrix.mat');
% gidcoord = [
% 
% 
% ];

%% Conectivities
% Element        Node(1)                Node(2)                Node(3)                Node(4)                Material
load('connectivityMatrix.mat');
% gidlnods = [
% 
% 
% ];

%% Variable Prescribed
% Node            Dimension                Value
load('prescribedNode.mat');
% lnodes = [
% 1 1 0 
% 1 2 0 
% 3 1 0 
% 3 2 0 
% 6 1 0 
% 6 2 0 
% 11 1 0 
% 11 2 0 
% 17 1 0 
% 17 2 0 
% 26 1 0 
% 26 2 0 
% 35 1 0 
% 35 2 0 
% 44 1 0 
% 44 2 0 
% 57 1 0 
% 57 2 0 
% 72 1 0 
% 72 2 0 
% 88 1 0 
% 88 2 0 
% 105 1 0 
% 105 2 0 
% 123 1 0 
% 123 2 0 
% 146 1 0 
% 146 2 0 
% 168 1 0 
% 168 2 0 
% 193 1 0 
% 193 2 0 
% 216 1 0 
% 216 2 0 
% 243 1 0 
% 243 2 0 
% 271 1 0 
% 271 2 0 
% 301 1 0 
% 301 2 0 
% 334 1 0 
% 334 2 0 
% 364 1 0 
% 364 2 0 
% 401 1 0 
% 401 2 0 
% 436 1 0 
% 436 2 0 
% 473 1 0 
% 473 2 0 
% 514 1 0 
% 514 2 0 
% 557 1 0 
% 557 2 0 
% 600 1 0 
% 600 2 0 
% 642 1 0 
% 642 2 0 
% 687 1 0 
% 687 2 0 
% 736 1 0 
% 736 2 0 
% 782 1 0 
% 782 2 0 
% 834 1 0 
% 834 2 0 
% 885 1 0 
% 885 2 0 
% 940 1 0 
% 940 2 0 
% 999 1 0 
% 999 2 0 
% 1050 1 0 
% 1050 2 0 
% 1111 1 0 
% 1111 2 0 
% 1167 1 0 
% 1167 2 0 
% 1233 1 0 
% 1233 2 0 
% 1297 1 0 
% 1297 2 0 
% 1357 1 0 
% 1357 2 0 
% 1424 1 0 
% 1424 2 0 
% 1491 1 0 
% 1491 2 0 
% 1563 1 0 
% 1563 2 0 
% 1635 1 0 
% 1635 2 0 
% 1702 1 0 
% 1702 2 0 
% 1778 1 0 
% 1778 2 0 
% 1852 1 0 
% 1852 2 0 
% 1931 1 0 
% 1931 2 0 
% 2011 1 0 
% 2011 2 0 
% 2095 1 0 
% 2095 2 0 
% 2176 1 0 
% 2176 2 0 
% 2256 1 0 
% 2256 2 0 
% 2341 1 0 
% 2341 2 0 
% 2425 1 0 
% 2425 2 0 
% 2518 1 0 
% 2518 2 0 
% 2605 1 0 
% 2605 2 0 
% 2698 1 0 
% 2698 2 0 
% 2788 1 0 
% 2788 2 0 
% 2883 1 0 
% 2883 2 0 
% 2981 1 0 
% 2981 2 0 
% 3078 1 0 
% 3078 2 0 
% 3177 1 0 
% 3177 2 0 
% 3278 1 0 
% 3278 2 0 
% 3384 1 0 
% 3384 2 0 
% 3484 1 0 
% 3484 2 0 
% 3585 1 0 
% 3585 2 0 
% 3693 1 0 
% 3693 2 0 
% 3807 1 0 
% 3807 2 0 
% 3914 1 0 
% 3914 2 0 
% 4025 1 0 
% 4025 2 0 
% 4133 1 0 
% 4133 2 0 
% 4255 1 0 
% 4255 2 0 
% 4373 1 0 
% 4373 2 0 
% 4491 1 0 
% 4491 2 0 
% 4608 1 0 
% 4608 2 0 
% 4729 1 0 
% 4729 2 0 
% 4856 1 0 
% 4856 2 0 
% 4974 1 0 
% 4974 2 0 
% 5101 1 0 
% 5101 2 0 
% 5230 1 0 
% 5230 2 0 
% 5357 1 0 
% 5357 2 0 
% 5491 1 0 
% 5491 2 0 
% 5618 1 0 
% 5618 2 0 
% 5758 1 0 
% 5758 2 0 
% 5891 1 0 
% 5891 2 0 
% 6030 1 0 
% 6030 2 0 
% 6167 1 0 
% 6167 2 0 
% 6301 1 0 
% 6301 2 0 
% 6452 1 0 
% 6452 2 0 
% 6591 1 0 
% 6591 2 0 
% 6734 1 0 
% 6734 2 0 
% 6879 1 0 
% 6879 2 0 
% 7027 1 0 
% 7027 2 0 
% 7181 1 0 
% 7181 2 0 
% 7326 1 0 
% 7326 2 0 
% 7479 1 0 
% 7479 2 0 
% 7636 1 0 
% 7636 2 0 
% 7789 1 0 
% 7789 2 0 
% 7953 1 0 
% 7953 2 0 
% ];

%% Force Prescribed
% Node                Dimension                Value
load('imposedForce.mat');
% pointload_complete = [
% 19634 2 -1 
% 19650 2 -1 
% 19661 2 -1 
% 19675 2 -1 
% 19688 2 -1 
% 19698 2 -1 
% 19714 2 -1 
% 19725 2 -1 
% 19741 2 -1 
% 19752 2 -1 
% 19769 2 -1 
% 19778 2 -1 
% 19797 2 -1 
% 19807 2 -1 
% 19824 2 -1 
% 19832 2 -1 
% 19850 2 -1 
% 19863 2 -1 
% 19878 2 -1 
% 19892 2 -1 
% 19901 2 -1 
% 19919 2 -1 
% 19929 2 -1 
% 19945 2 -1 
% 19959 2 -1 
% 19970 2 -1 
% 19986 2 -1 
% 19998 2 -1 
% 20009 2 -1 
% 20023 2 -1 
% 20035 2 -1 
% ];

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
