%% Proof of concept
% Using functions!
clc; clear; close all

% Create the data container for the FEM problem
% a.fileName = 'holeinclusion3d';
a.fileName = 'test3d_micro_cube';
m = FemDataContainer(a);

m.bc.dirichlet = [641 1 0;
    641 2 0;
    641 3 0;
    1250 1 0;
    1250 2 0;
    1250 3 0;
    1251 1 0;
    1251 2 0;
    1251 3 0;
    1331 1 0;
    1331 2 0;
    1331 3 0;
    ];

m.bc.pointload = [245 1 1000];

% Create the characteristic function (1 inside circle, 0 outside)
s.mesh    = m.mesh;
s.fxy     = @(x,y,z) (y-0.5).^2+(z-0.5).^2 -0.1.^2;
circleFun = CharacteristicFunction(s);

% Project the function to P0. Useful later on
x.mesh   = m.mesh;
x.connec = m.mesh.connec;
projP0 = Projector_toP0(x);
p0c = projP0.project(circleFun);

% Generate the hole in the material using the values we just found
fV = squeeze(p0c.fValues);
holeNodes = find(fV==1);
m.material.C(:,:,holeNodes) = m.material.C(:,:, holeNodes)*1e3;


% Solve the problem
m.scale = 'MACRO';
fem = FEM.create(m);
fem.solve();