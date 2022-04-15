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
7 0.333333 0 0
8 0.666667 0 0
9 1.16667 0.288675 0
10 1.33333 0.57735 0
11 1.33333 1.1547 0
12 1.16667 1.44338 0
13 0.666667 1.73205 0
14 0.333333 1.73205 0
15 -0.166667 1.44338 0
16 -0.333333 1.1547 0
17 -0.333333 0.57735 0
18 -0.166667 0.288675 0
19 0.5 0.866025 0
20 0.166667 0.288675 0
21 0.833333 0.288675 0
22 1.16667 0.866025 0
23 0.833333 1.44338 0
24 0.166667 1.44338 0
25 -0.166667 0.866025 0
26 0.388889 0.288675 0
27 0.611111 0.288675 0
28 1 0.57735 0
29 1 1.1547 0
30 0.611111 1.44338 0
31 0.388889 1.44338 0
32 -2.77556e-17 1.1547 0
33 0 0.57735 0
34 0.333333 0.57735 0
35 0.666667 0.57735 0
36 0.833333 0.866025 0
37 0.666667 1.1547 0
38 0.333333 1.1547 0
39 0.166667 0.866025 0
40 0.5 0.57735 0
41 0.5 1.1547 0
];

pointload = [
];



connec = [
1 25 33 39
2 19 35 36
3 8 2 21
4 20 7 26
5 34 19 39
6 19 36 37
7 29 37 36
8 35 27 21
9 12 11 4
10 9 21 2
11 9 28 21
12 28 35 21
13 21 27 8
14 29 36 22
15 17 18 33
16 26 7 8
17 40 34 26
18 20 1 7
19 34 20 26
20 33 18 20
21 20 34 33
22 26 27 40
23 25 6 17
24 18 1 20
25 33 25 17
26 16 15 6
27 39 32 25
28 24 14 5
29 25 16 6
30 32 24 15
31 15 5 6
32 31 38 41
33 15 24 5
34 15 16 32
35 25 32 16
36 38 24 32
37 34 39 33
38 38 32 39
39 31 30 14
40 23 37 29
41 31 41 30
42 41 19 37
43 24 31 14
44 24 38 31
45 30 13 14
46 30 23 13
47 13 23 4
48 30 37 23
49 23 12 4
50 23 29 12
51 22 28 10
52 29 11 12
53 29 22 11
54 36 35 28
55 11 3 4
56 11 22 3
57 39 19 38
58 34 40 19
59 30 41 37
60 38 19 41
61 35 40 27
62 35 19 40
63 26 8 27
64 10 28 9
65 22 36 28
66 22 10 3
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
7 14
8 13
9 16
10 15
11 18
12 17
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