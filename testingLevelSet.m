% %% Generating a 3D mesh with a hole inclusion
% % Using functions!
% clc; clear; close all
% 
% % Create the data container for the FEM problem
% % a.fileName = 'holeinclusion3d';
% a.fileName = 'test3d_micro_cube';
% m = FemDataContainer(a);
% 
% 
% % Create the characteristic function (1 inside circle, 0 outside)
% s.mesh    = m.mesh;
% s.fxy     = @(x,y,z) (x-0.5).^2+(y-0.5).^2+(z-0.5).^2 -0.3.^2;
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