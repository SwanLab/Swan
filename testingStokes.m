%% Testing other functions
% Create a Mesh FEM results
clear; close all;

file = 'test2d_triangle';
a.fileName = file;
s = FemDataContainer(a);
mesh = s.mesh;

%% Create functions
% AnalyticalFunction

sAF.fHandle = @(x) x(1,:,:);
sAF.ndimf   = 1;
sAF.mesh    = mesh;
xFun = AnalyticalFunction(sAF);

p1 = xFun.project('P1');
p1.plot

q = Quadrature.create(mesh,'LINEAR');
p1.evaluate(q.posgp)

% get dofs from logical analytical function?
leftCond = @(x) x(:,1) == 0;