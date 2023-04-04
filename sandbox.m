%% This is a sandbox file!
% Feel free to test anything here :)
clc; clear; close all;

% file = 'test2d_triangle';
% file = 'test2d_quad';
% file = 'test3d_hexahedra';
file = 'test2d_micro';
a.fileName = file;
s = FemDataContainer(a);
mesh = s.mesh;
fem = FEM.create(s);
fem.solve();
