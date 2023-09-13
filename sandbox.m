%% This is a sandbox file!
% % Feel free to test anything here :)
% % clc; clear; close all;
% 
% % Load mesh
% file = 'test2d_micro';
% a.fileName = file;
% s = FemDataContainer(a);
% mesh = s.mesh;
% clear s;
% 
% % Create functions
% 
% sAF.fHandle = @(x) x(1,:,:);
% sAF.ndimf   = 1;
% sAF.mesh    = mesh;
% xFun = AnalyticalFunction(sAF);
% 
% p1trial = xFun.project('P1');
% % p1trial = P1Function.create(mesh, 1);
% p0test  = P0Function.create(mesh, 1);
% 
% % LHS integrator
% 
% s.type = 'MassTestTrial';
% s.mesh = mesh;
% s.test = p0test;
% s.trial = p1trial;
% lhs = LHSintegrator.create(s);
% LHS = lhs.compute();
% 
% % Mass P0
% 
% s.type = 'MassTestTrial';
% s.mesh = mesh;
% s.test = p0test;
% s.trial = p0test;
% mp0 = LHSintegrator.create(s);
% MP0 = mp0.compute();
% 
% % Testing
% 
% gj = MP0\(LHS*p1trial.fValues);
% 
% z.fValues = gj;
% z.mesh = mesh;
% p0_result = P0Function(z);

% %% Generating a 2D mesh with a hole inclusion
% % Using functions!
% clear; close all
% 
% % Create the data container for the FEM problem
% a.fileName = 'test2d_micro';
% m = FemDataContainer(a);
% 
% % Create the characteristic function (1 inside circle, 0 outside)
% s.mesh    = m.mesh;
% s.fxy     = @(x,y) (x-0.5).^2+(y-0.5).^2-0.3.^2;
% circleFun = CharacteristicFunction(s);

%% Generating a 3D mesh with a hole inclusion
% Using functions!
clc; clear; close all

% Create the data container for the FEM problem
% a.fileName = 'holeinclusion3d';
a.fileName = 'test3d_micro_cube';
m = FemDataContainer(a);


% Create the characteristic function (1 inside circle, 0 outside)
s.mesh    = m.mesh;
s.fxy     = @(x,y,z) (x-0.5).^2+(y-0.5).^2+(z-0.5).^2 -0.3.^2;
circleFun = CharacteristicFunction(s);
