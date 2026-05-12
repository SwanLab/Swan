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
1            0            0            0
2            0      0.00078            0
3            0      0.00156            0
4            0      0.00234            0
5            0      0.00312            0
6       0.0102            0            0
7       0.0102      0.00078            0
8       0.0102      0.00156            0
9       0.0102      0.00234            0
10       0.0102      0.00312            0
11       0.0204            0            0
12       0.0204      0.00078            0
13       0.0204      0.00156            0
14       0.0204      0.00234            0
15       0.0204      0.00312            0
16       0.0306            0            0
17       0.0306      0.00078            0
18       0.0306      0.00156            0
19       0.0306      0.00234            0
20       0.0306      0.00312            0
21       0.0408            0            0
22       0.0408      0.00078            0
23       0.0408      0.00156            0
24       0.0408      0.00234            0
25       0.0408      0.00312            0
26        0.051            0            0
27        0.051      0.00078            0
28        0.051      0.00156            0
29        0.051      0.00234            0
30        0.051      0.00312            0
31       0.0612            0            0
32       0.0612      0.00078            0
33       0.0612      0.00156            0
34       0.0612      0.00234            0
35       0.0612      0.00312            0
36       0.0714            0            0
37       0.0714      0.00078            0
38       0.0714      0.00156            0
39       0.0714      0.00234            0
40       0.0714      0.00312            0
41       0.0816            0            0
42       0.0816      0.00078            0
43       0.0816      0.00156            0
44       0.0816      0.00156            0
45       0.0816      0.00234            0
46       0.0816      0.00312            0
47       0.0918            0            0
48       0.0918      0.00078            0
49       0.0918      0.00156            0
50       0.0918      0.00156            0
51       0.0918      0.00234            0
52       0.0918      0.00312            0
53        0.102            0            0
54        0.102      0.00078            0
55        0.102      0.00156            0
56        0.102      0.00156            0
57        0.102      0.00234            0
58        0.102      0.00312            0
];

%% Conectivities
% Element        Node(1)                Node(2)                Node(3)                Node(4)                Material

gidlnods = [
1 6 7 2 1 0
2 11 12 7 6 0
3 16 17 12 11 0
4 21 22 17 16 0
5 26 27 22 21 0
6 31 32 27 26 0
7 36 37 32 31 0
8 41 42 37 36 0
9 47 48 42 41 0
10 53 54 48 47 0
11 7 8 3 2 0
12 12 13 8 7 0
13 17 18 13 12 0
14 22 23 18 17 0
15 27 28 23 22 0
16 32 33 28 27 0
17 37 38 33 32 0
18 42 43 38 37 0
19 48 49 43 42 0
20 54 55 49 48 0
21 8 9 4 3 0
22 13 14 9 8 0
23 18 19 14 13 0
24 23 24 19 18 0
25 28 29 24 23 0
26 33 34 29 28 0
27 38 39 34 33 0
28 44 45 39 38 0
29 50 51 45 44 0
30 56 57 51 50 0
31 9 10 5 4 0
32 14 15 10 9 0
33 19 20 15 14 0
34 24 25 20 19 0
35 29 30 25 24 0
36 34 35 30 29 0
37 39 40 35 34 0
38 45 46 40 39 0
39 51 52 46 45 0
40 57 58 52 51 0
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
