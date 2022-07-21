%% Data file

Data_prb = {
'TRIANGLE' ;
'SI' ;
'2D' ;
'Plane_Stress' ;
'ELASTIC' ;
'MICRO' ;
};


coord = [
1 0 0 0
2 1 0 0
3 1.5 0.866025 0
4 1 1.73205 0
5 0 1.73205 0
6 -0.5 0.866025 0
7 0.5 0 0
8 1.25 0.433013 0
9 1.25 1.29904 0
10 0.5 1.73205 0
11 -0.25 1.29904 0
12 -0.25 0.433013 0
13 0.5 0.866025 0
14 0.25 0.433013 0
15 0.75 0.433013 0
16 1 0.866025 0
17 0.75 1.29904 0
18 0.25 1.29904 0
19 0 0.866025 0
20 0.875 1.08253 0
21 0.125 0.649519 0
];

pointload = [
];



connec = [
1 19 6 12
2 14 13 21
3 6 19 11
4 19 12 21
5 14 15 13
6 14 1 7
7 21 12 14
8 12 1 14
9 14 7 15
10 21 13 19
11 17 4 10
12 13 18 19
13 19 18 11
14 11 18 5
15 13 15 16
16 18 10 5
17 18 13 17
18 18 17 10
19 13 20 17
20 20 9 17
21 17 9 4
22 20 13 16
23 16 15 8
24 15 7 2
25 20 16 9
26 15 2 8
27 16 8 3
28 16 3 9
];

%% Variable Prescribed
% Node 	 Dimension 	 Value

dirichlet_data = [
1 1 0
1 2 0
2 1 0
2 2 0
3 1 0
3 2 0
4 1 0
4 2 0
5 1 0
5 2 0
6 1 0
6 2 0
];


%% Force Prescribed
% Node 	 Dimension 	 Value

pointload_complete = [
];


%% Volumetric Force
% Element 	 Dimension 	 Force_Dim

Vol_force = [
];


%% Group Elements
% Element 	 Group_num

Group = [
];


%% Group Elements
% Elements that are considered holes initially
% Element

Initial_holes = [
];


%% Boundary Elements
% Elements that can not be removed
% Element

Boundary_elements = [
];


%% Micro Gauss Post
% Element

Micro_gauss_post = [
];


%% Master-Slave

Master_slave = [
7 10
8 11
9 12
];

%% Nodes Solid
% Nodes that must remain
% Nodes

% nodesolid = 1;


%% External Border Elements
% Detect the elements that define the edge of the domain
% Element 	 Node(1) 	 Node(2)

External_border_elements = [
];


%% External Border Nodes
% Detect the nodes that define the edge of the domain
% Node 	 

External_border_nodes = [
];


%% Materials
% Materials that have been used
% Material_Num 	 Mat_density 	 Young_Modulus 	 Poisson

Materials = [
];