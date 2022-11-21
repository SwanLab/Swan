%% Testing ISCSO FEM
clc; clear;

run('ISCSOMesh.m')
s.coord  = coord;
s.connec = connec;
s.material = ISCSOMaterial(s);

fem = StructuralTrussProblem(s);
fem.solve();