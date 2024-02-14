%% Mesh refactoring
clc; clear; close all;

% 2D Mesh
file = 'test2d_triangle';
a.fileName = file;
s = FemDataContainer(a);
mesh2d = s.mesh;


% 3D Mesh
file = 'test3d_hexahedra';
a.fileName = file;
s = FemDataContainer(a);
mesh3d = s.mesh;
