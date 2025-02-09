%% Generating a 2D mesh with a hole inclusion
% Using functions!
clear; close all

% Create the data container for the FEM problem
a.fileName = 'test2d_micro';
m = FemDataContainer(a);

% Create the characteristic function (1 inside circle, 0 outside)
s.mesh    = m.mesh;
s.fxy     = @(x,y) (x-0.5).^2+(y-0.5).^2-0.3.^2;
circleFun = CharacteristicFunction(s);

% Project the function to P0. Useful later on, also can be plotted
x.mesh   = m.mesh;
x.connec = m.mesh.connec;
projP0 = Projector_toP0(x);
p0c = projP0.project(circleFun);
p0c.plot();

% Generate the hole in the material using the values we just found
fV = squeeze(p0c.fValues);
holeNodes = find(fV==1);
% m.material.C(:,:,holeNodes) = zeros(3,3, length(holeNodes));
m.material.C(:,:,holeNodes) = m.material.C(:,:, holeNodes)*1e-10;
% m.material.C(:,:,holeNodes) = repmat([10*10e-3, 10e-3, 0; 10e-3 10*10e-3 0; 0 0 10e-3], [1 1 length(holeNodes)]);


% sUM.backgroundMesh = m.mesh;
% sUM.boundaryMesh = m.mesh.createBoundaryMesh();
% uM = UnfittedMesh(sUM);

% Solve the problem
fem = ElasticProblemMicro(m);
fem.computeChomog();

fem.uFun{1}.plot
fem.uFun{2}.plot
fem.uFun{3}.plot

%% Generating a 3D mesh with a hole inclusion
% Using functions!
clc; clear; close all

% Create the data container for the FEM problem
% a.fileName = 'holeinclusion3d';
% a.fileName = 'test3d_micro_cube';
a.fileName = 'test3d_micro_cube_v2';
m = FemDataContainer(a);


% Create the characteristic function (1 inside circle, 0 outside)
s.mesh    = m.mesh;
s.fxy     = @(x,y,z) (x-0.5).^2+(y-0.5).^2+(z-0.5).^2 -0.3.^2;
circleFun = CharacteristicFunction(s);

% Project the function to P0. Useful later on
x.mesh   = m.mesh;
x.connec = m.mesh.connec;
projP0 = Projector_toP0(x);
p0c = projP0.project(circleFun);

% Generate the hole in the material using the values we just found
fV = squeeze(p0c.fValues);
holeNodes = find(fV==1);
m.material.C(:,:,holeNodes) = m.material.C(:,:, holeNodes)*1e-3;


% % Solve the problem
fem = ElasticProblemMicro_Fast(m);
fem.solve();
% sss.filename = 'fluctHoleGID';
sss.filename = 'fluctHoleMaterial';
fem.uFun{3}.print(sss);
