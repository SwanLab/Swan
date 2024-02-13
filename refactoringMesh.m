%% Mesh refactoring
clc; clear; close all;

file = 'test3d_hexahedra';
a.fileName = file;
s = FemDataContainer(a);

b.coord = s.mesh.coord;
b.connec = s.mesh.connec;
mesh = Mesh.create(b);