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
% view(0,90)
title('p1')

resCharFunInP1part2 = projP1.project(resCharFunInP1);
resCharFunInP1part2.plot(m.mesh);
title('p1 part 2')

% CharFun to P0 Function
projP0 = Projector_toP0(x);
resCharFunInP0 = projP0.project(circleFun);
resCharFunInP0.plot(m.mesh);
title('p0')


% CharFun to P1 Discontinuous Function
projP1D = Projector_toP1Discontinuous(x);
resCharFunInP1D = projP1D.project(circleFun);
resCharFunInP1D.plot(m.mesh);
title('p1d')

resCharFunInP1dpart2 = projP1D.project(resCharFunInP1D);
resCharFunInP1dpart2.plot(m.mesh);
title('p1 part 2')