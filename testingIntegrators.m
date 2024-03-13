clc; clear; close all;

% Test and trial
mesh = UnitTriangleMesh(10,10);
quad = Quadrature.set(mesh.type);
quad.computeQuadrature('LINEAR');

test  = TestFunction(mesh, 1, 'P1');
trial = TrialFunction(mesh, 1, 'P1');


u = LagrangianFunction.create(mesh, 1, 'P1');
grad = Grad(u).evaluate(quad.posgp);