%% Goal: delete RHS integrators
% Previous goal: delete BMatrixComputer
clc; clear; close all;

% Test and trial
mesh = UnitTriangleMesh(10,10);
quad = Quadrature.set(mesh.type);
quad.computeQuadrature('QUADRATIC');
xV = quad.posgp;

ndimf = 2;
trial = TrialFunction.create(mesh, ndimf, 'P1');
test  = TestFunction.create(mesh, ndimf, 'P1');

% Analytical function to integrate: y=x
sAF.fHandle = @(x) [x(1,:,:); 2*x(1,:,:)];
sAF.ndimf   = 2;
sAF.mesh    = mesh;
xFun = AnalyticalFunction(sAF); % 1 x 3 x 324
% 2 x 3 x 3 x 324
fun = xFun.project('P1');

dNdx = ShapeDer(test).evaluate(xV);
fG = xFun.evaluate(xV);
fP1g = fun.evaluate(xV);
dV(1,1,:,:) = mesh.computeDvolume(quad);



s.type = 'ShapeFunction';%'ShapeSymmetricDerivative';%'ShapeFunction';%'ShapeDerivative';
s.mesh = mesh;
s.quadratureOrder = 'QUADRATIC';
s.quadType = 'QUADRATIC';
rhs = RHSintegrator.create(s);
RHS = rhs.compute(xFun, test);