%% P0-P1 Mass matrix example
clc; clear; close all;

% Load mesh
file = 'test2d_micro';
a.fileName = file;
s = FemDataContainer(a);
mesh = s.mesh;
clear s;

% Create functions

sAF.fHandle = @(x) [x(1,:,:);x(2,:,:).*x(1,:,:)];
sAF.ndimf   = 2;
sAF.mesh    = mesh;
xFun = AnalyticalFunction(sAF);

p1trial = xFun.project('P1');
%p1trial = P1Function.create(mesh, 2);
p0test  = LagrangianFunction.create(mesh, 2, 'P0');

% LHS integrator

s.type = 'MassMatrix';
s.mesh = mesh;
s.test = p0test;
s.trial = p1trial;
lhs = LHSintegrator.create(s);
LHS = lhs.compute();

% Mass P0

s.type = 'MassMatrix';
s.mesh = mesh;
s.test = p0test;
s.trial = p0test;
s.quadratureOrder = 'QUADRATIC';
mp0 = LHSintegrator.create(s);
MP0 = mp0.compute();

% Result
p1V = p1trial.fValues';
p1V = p1V(:);
gj = MP0\(LHS*p1V);

gj = reshape(gj,2,[])';
% Plot
z.fValues = gj;
z.mesh = mesh;
z.order = 'P0';
p0_result = LagrangianFunction(z);
p0_result.plot();