close all;clear all

x1 = linspace(0,1,5);
x2 = linspace(0,1,5);
[xv,yv] = meshgrid(x1,x2);
[F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');
sBg.coord  = V(:,1:2);
sBg.connec = F;
mesh = Mesh(sBg);
figure
mesh.plot()

sAF.fHandle = @(x) [x(1,:,:).^2;x(2,:,:).^2]; % f(x) = x
sAF.ndimf   = 2; % number of dimensions
sAF.mesh    = mesh;
xFun = AnalyticalFunction(sAF);
xV(1,:,:)=1/3;
xV(2,:,:)=1/3;
xFun.evaluate(xV);

refpoint = [0.5 0.5];

RBfun = xFun.project('RigidBody',refpoint);
P1=RBfun.project('P1');
P1.plot()
% P1.print('RB')

P1Fun=xFun.project('P1');
P1Fun.plot()