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
2 0 0.166667 0
3 0 0.333333 0
4 0 0.5 0
5 0 0.666667 0
6 0 0.833333 0
7 0 1 0
8 0.166667 0 0
9 0.166667 0.166667 0
10 0.166667 0.333333 0
11 0.166667 0.5 0
12 0.166667 0.666667 0
13 0.166667 0.833333 0
14 0.166667 1 0
15 0.333333 0 0
16 0.333333 0.166667 0
17 0.333333 0.333333 0
18 0.333333 0.5 0
19 0.333333 0.666667 0
20 0.333333 0.833333 0
21 0.333333 1 0
22 0.5 0 0
23 0.5 0.166667 0
24 0.5 0.333333 0
25 0.5 0.5 0
26 0.5 0.666667 0
27 0.5 0.833333 0
28 0.5 1 0
29 0.666667 0 0
30 0.666667 0.166667 0
31 0.666667 0.333333 0
32 0.666667 0.5 0
33 0.666667 0.666667 0
34 0.666667 0.833333 0
35 0.666667 1 0
36 0.833333 0 0
37 0.833333 0.166667 0
38 0.833333 0.333333 0
39 0.833333 0.5 0
40 0.833333 0.666667 0
41 0.833333 0.833333 0
42 0.833333 1 0
43 1 0 0
44 1 0.166667 0
45 1 0.333333 0
46 1 0.5 0
47 1 0.666667 0
48 1 0.833333 0
49 1 1 0
];

pointload = [
];



connec = [
1 10 4 3
2 3 2 9
3 11 4 10
4 13 14 7
5 26 27 20
6 20 14 13
7 20 21 14
8 28 21 27
9 16 15 22
10 34 27 33
11 18 17 24
12 2 1 8
13 39 40 33
14 47 48 41
15 40 41 34
16 38 45 39
17 41 35 34
18 9 2 8
19 30 31 24
20 9 8 15
21 39 33 32
22 16 9 15
23 24 17 23
24 17 11 10
25 29 30 23
26 25 26 19
27 10 3 9
28 16 10 9
29 17 18 11
30 17 16 23
31 17 10 16
32 12 5 11
33 26 20 19
34 5 4 11
35 7 6 13
36 5 12 6
37 25 18 24
38 12 19 13
39 27 21 20
40 18 12 11
41 13 6 12
42 13 19 20
43 12 18 19
44 18 25 19
45 33 27 26
46 40 34 33
47 35 28 34
48 32 31 38
49 37 36 43
50 34 28 27
51 41 42 35
52 39 46 40
53 49 42 48
54 47 41 40
55 48 42 41
56 45 46 39
57 46 47 40
58 38 37 44
59 30 24 23
60 44 45 38
61 31 25 24
62 32 26 25
63 38 31 37
64 38 39 32
65 26 32 33
66 25 31 32
67 23 16 22
68 37 30 36
69 23 22 29
70 36 30 29
71 37 31 30
72 43 44 37
];

%% Variable Prescribed
% Node 	 Dimension 	 Value

dirichlet_data = [
1 1 0
1 2 0
7 1 0
7 2 0
43 1 0
43 2 0
49 1 0
49 2 0
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
2 44
3 45
4 46
5 47
6 48
8 14
15 21
22 28
29 35
36 42
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