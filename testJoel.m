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
2 2 0 0
3 0.75253 5.86889 0
4 -1.24747 5.86889 0
5 1 0 0
6 1.79209 0.978148 0
7 1.58418 1.9563 0
8 1.37626 2.93444 0
9 1.16835 3.91259 0
10 0.960442 4.89074 0
11 -0.24747 5.86889 0
12 -1.03956 4.89074 0
13 -0.831647 3.91259 0
14 -0.623735 2.93444 0
15 -0.415823 1.9563 0
16 -0.207912 0.978148 0
17 0.792088 0.978148 0
18 0.584177 1.9563 0
19 0.376265 2.93444 0
20 0.168353 3.91259 0
21 -0.0395585 4.89074 0
];

pointload = [
];



connec = [
1 14 15 19
2 20 19 9
3 20 14 19
4 13 14 20
5 17 1 5
6 12 13 21
7 17 16 1
8 15 18 19
9 15 16 18
10 11 12 21
11 20 10 21
12 13 20 21
13 21 3 11
14 21 10 3
15 19 18 8
16 10 20 9
17 9 19 8
18 7 17 6
19 18 16 17
20 6 17 5
21 18 7 8
22 18 17 7
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
5 11
6 16
7 15
8 14
9 13
10 12
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