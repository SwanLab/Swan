%% Testing characteristic function
% Get mesh
clear; close all;

file = 'RVE_Square_Triangle';
a.fileName = file;
m = FemDataContainer(a);

% Create characteristic function
s.mesh    = m.mesh;
s.fxy     = @(x,y) (x-0.5)^2+(y-0.5)^2-0.3^2;
circleFun = CharacteristicFunction(s);

% CharFun to P1 Function
x.mesh   = m.mesh;
x.connec = m.mesh.connec;
projP1 = Projector_toP1(x);
resCharFunInP1 = projP1.project(circleFun);
resCharFunInP1.plot(m.mesh);