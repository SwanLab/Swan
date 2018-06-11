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
1            1            0            0
2            1      0.14286            0
3      0.85714            0            0
4      0.85625      0.14462            0
5            1      0.28571            0
6      0.71429            0            0
7      0.85604      0.28799            0
8      0.71244      0.14631            0
9       0.7123      0.29003            0
10            1      0.42857            0
11      0.57143            0            0
12      0.57021      0.14589            0
13      0.85563      0.43159            0
14      0.56999       0.2891            0
15      0.71238      0.43201            0
16            1      0.57143            0
17      0.42857            0            0
18      0.42794      0.14486            0
19      0.85527      0.57478            0
20      0.56953      0.43206            0
21      0.42792      0.28844            0
22      0.71205      0.57434            0
23            1      0.71429            0
24      0.28571            0            0
25      0.42686      0.42976            0
26      0.56941      0.57309            0
27      0.28507      0.14513            0
28      0.85557      0.71657            0
29      0.28394      0.28752            0
30      0.71234      0.71652            0
31      0.42654      0.57222            0
32      0.28289      0.42882            0
33      0.56898      0.71622            0
34            1      0.85714            0
35      0.14286            0            0
36      0.14196      0.14392            0
37      0.85575      0.85867            0
38       0.1417      0.28659            0
39      0.71276      0.85897            0
40      0.42648      0.71477            0
41       0.2834      0.57126            0
42       0.1409       0.4292            0
43      0.56973      0.85866            0
44            0            0            0
45            1            1            0
46      0.85714            1            0
47            0      0.14286            0
48      0.28279      0.71474            0
49      0.14083      0.57174            0
50      0.42707      0.85839            0
51      0.71429            1            0
52            0      0.28571            0
53            0      0.42857            0
54      0.57143            1            0
55      0.14186      0.71522            0
56      0.28389      0.85826            0
57      0.42857            1            0
58            0      0.57143            0
59      0.14215      0.85766            0
60      0.28571            1            0
61            0      0.71429            0
62      0.14286            1            0
63            0      0.85714            0
64            0            1            0
];

%% Conectivities
% Element        Node(1)                Node(2)                Node(3)                Node(4)                Material

connec = [
1 47 44 35 36 0
2 3 1 2 4 0
3 34 45 46 37 0
4 62 64 63 59 0
5 52 47 36 38 0
6 6 3 4 8 0
7 23 34 37 28 0
8 60 62 59 56 0
9 53 52 38 42 0
10 35 24 27 36 0
11 11 6 8 12 0
12 2 5 7 4 0
13 16 23 28 19 0
14 46 51 39 37 0
15 57 60 56 50 0
16 63 61 55 59 0
17 58 53 42 49 0
18 61 58 49 55 0
19 24 17 18 27 0
20 17 11 12 18 0
21 5 10 13 7 0
22 10 16 19 13 0
23 51 54 43 39 0
24 54 57 50 43 0
25 4 7 9 8 0
26 38 36 27 29 0
27 27 18 21 29 0
28 37 39 30 28 0
29 19 28 30 22 0
30 13 19 22 15 0
31 7 13 15 9 0
32 39 43 33 30 0
33 18 12 14 21 0
34 12 8 9 14 0
35 42 38 29 32 0
36 49 42 32 41 0
37 59 55 48 56 0
38 55 49 41 48 0
39 50 56 48 40 0
40 43 50 40 33 0
41 9 15 20 14 0
42 22 30 33 26 0
43 15 22 26 20 0
44 29 21 25 32 0
45 21 14 20 25 0
46 48 41 31 40 0
47 40 31 26 33 0
48 20 26 31 25 0
49 32 25 31 41 0
];

%% Variable Prescribed
% Node            Dimension                Value

dirichlet_data = [
44 1 0 
44 2 0 
47 1 0 
47 2 0 
52 1 0 
52 2 0 
53 1 0 
53 2 0 
58 1 0 
58 2 0 
61 1 0 
61 2 0 
63 1 0 
63 2 0 
64 1 0 
64 2 0 
];

%% Force Prescribed
% Node                Dimension                Value

pointload_complete = [
1 1 1 
2 1 1 
5 1 1 
10 1 1 
16 1 1 
23 1 1 
34 1 1 
45 1 1 
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
