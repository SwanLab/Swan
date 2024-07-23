%==================================================================
%                        General Data File
% Title: TETRAHEDRA
% Units: SI
% Dimensions: 2D
% Type of problem: Plane_Stress
% Type of Phisics: ELASTIC
% Micro/Macro: MACRO
%
%==================================================================

%% Data

Data_prb = {
'TRIANGLE';
'SI';
'2D';
'Plane_Stress';
'ELASTIC';
'MACRO';
};

%% Coordinates
% Node                X                Y                Z

coord = [
   1        00000        00000 
     2  1.00000e-02  1.00000e-02 
     3        00000  2.00000e-02 
     4  2.00000e-02        00000 
     5  2.00000e-02  2.00000e-02 
     6  1.00000e-02  3.00000e-02 
     7  3.00000e-02  1.00000e-02 
     8        00000  4.00000e-02 
     9  4.00000e-02        00000 
    10  3.00000e-02  3.00000e-02 
    11  4.00000e-02  2.00000e-02 
    12  2.00000e-02  4.00000e-02 
    13  4.00000e-02  4.00000e-02 
];

%% Conectivities
% Element        Node(1)                Node(2)                Node(3)                Node(4)                Material

connec = [
     1      5      2      4   
     2      3      2      5   
     3      1      2      3   
     4      4      2      1   
     5     12      6      5   
     6      8      6     12   
     7      3      6      8  
     8      5      6      3   
     9     11      7      9   
    10      5      7     11   
    11      4      7      5   
    12      9      7      4   
    13     13     10     11   
    14     12     10     13   
    15      5     10     12   
    16     11     10      5  
];

%% Variable Prescribed
% Node            Dimension                Value

dirichlet_data = [
    1 1 1
    1 2 2
    3 1 1
    3 2 0
    8 1 1
    8 2 1
];

%% Force Prescribed
% Node                Dimension                Value

pointload_complete = [
11 2 -1
];

pointload_complete = [
11 2 -1
];

isLeft   = @(coor) (abs(coor(:,1) - min(coor(:,1)))   < 1e-12);
isRight  = @(coor) (abs(coor(:,1) - max(coor(:,1)))   < 1e-12);
isMiddle = @(coor) (abs(coor(:,2) - max(coor(:,2)/2)) == 0);
isTop    = @(coor) (abs(coor(:,2) - max(coor(:,2))) == 0);
isBottom = @(coor) (abs(coor(:,2) - min(coor(:,2))) == 0);

% Dirichlet
sDir{1}.domain    = @(coor) isBottom(coor) & isLeft(coor);
sDir{1}.direction = [1,2];
sDir{1}.value     = [1,2];

sDir{2}.domain    = @(coor) isMiddle(coor) & isLeft(coor);
sDir{2}.direction = [1,2];
sDir{2}.value     = [1,0];

sDir{3}.domain    = @(coor) isTop(coor) & isLeft(coor);
sDir{3}.direction = [1,2];
sDir{3}.value     = 1;

% Point load
sPL{1}.domain    = @(coor) isMiddle(coor) & isRight(coor);
sPL{1}.direction = 2;
sPL{1}.value     = -1;

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
