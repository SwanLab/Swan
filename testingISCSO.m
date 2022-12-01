%% Testing ISCSO FEM
clc; clear;

run('ISCSOMesh3D.m')
s.coord  = coord;
s.connec = connec;
s.material = ISCSOMaterial(s);
s.neumann = neumann;
s.dirichlet = dirichlet;

fem = StructuralTruss3DProblem(s);
fem.solve();