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
1          102            0            0
2          102         0.78            0
3          102         1.56            0
4          102         2.34            0
5          102         3.12            0
6         91.8            0            0
7         91.8         0.78            0
8         91.8         1.56            0
9         91.8         2.34            0
10         91.8         3.12            0
11         81.6            0            0
12         81.6         0.78            0
13         81.6         1.56            0
14         81.6         2.34            0
15         81.6         3.12            0
16         71.4            0            0
17         71.4         0.78            0
18         71.4         1.56            0
19         71.4         2.34            0
20         71.4         3.12            0
21         61.2            0            0
22         61.2         0.78            0
23         61.2         1.56            0
24         61.2         2.34            0
25         61.2         3.12            0
26           51            0            0
27           51         0.78            0
28           51         1.56            0
29           51         2.34            0
30           51         3.12            0
31         40.8            0            0
32         40.8         0.78            0
33         40.8         1.56            0
34         40.8         2.34            0
35         40.8         3.12            0
36         30.6            0            0
37         30.6         0.78            0
38         30.6         1.56            0
39         30.6         2.34            0
40         30.6         3.12            0
41         20.4            0            0
42         20.4         0.78            0
43         20.4         1.56            0
44         20.4         2.34            0
45         20.4         3.12            0
46         10.2            0            0
47         10.2         0.78            0
48         10.2         1.56            0
49         10.2         2.34            0
50         10.2         3.12            0
51            0            0            0
52            0         0.78            0
53            0         1.56            0
54            0         2.34            0
55            0         3.12            0
];

%% Conectivities
% Element        Node(1)                Node(2)                Node(3)                Node(4)                Material

gidlnods = [
1 46 47 52 51 0
2 41 42 47 46 0
3 36 37 42 41 0
4 31 32 37 36 0
5 26 27 32 31 0
6 21 22 27 26 0
7 16 17 22 21 0
8 11 12 17 16 0
9 6 7 12 11 0
10 1 2 7 6 0
11 47 48 53 52 0
12 42 43 48 47 0
13 37 38 43 42 0
14 32 33 38 37 0
15 27 28 33 32 0
16 22 23 28 27 0
17 17 18 23 22 0
18 12 13 18 17 0
19 7 8 13 12 0
20 2 3 8 7 0
21 50 55 54 49 0
22 45 50 49 44 0
23 40 45 44 39 0
24 35 40 39 34 0
25 30 35 34 29 0
26 25 30 29 24 0
27 20 25 24 19 0
28 15 20 19 14 0
29 10 15 14 9 0
30 5 10 9 4 0
31 49 54 53 48 0
32 44 49 48 43 0
33 39 44 43 38 0
34 34 39 38 33 0
35 29 34 33 28 0
36 24 29 28 23 0
37 19 24 23 18 0
38 14 19 18 13 0
39 9 14 13 8 0
40 4 9 8 3 0
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

if ~isempty(pointload_complete)
    nodesolid = unique(pointload_complete(:,1));
end

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
