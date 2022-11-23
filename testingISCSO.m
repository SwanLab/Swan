%% Testing ISCSO FEM
clc; clear;

run('ISCSOMesh.m')
s.coord  = coord;
s.connec = connec;
s.material = ISCSOMaterial(s);
% forces = [5 1 0;
%           5 2 10];
s.neumann = neumann;
s.dirichlet = dirichlet;

fem = StructuralTrussProblem(s);
fem.solve();