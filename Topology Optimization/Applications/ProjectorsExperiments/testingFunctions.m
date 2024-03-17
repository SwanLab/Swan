%% Testing other functions
% Create a Mesh FEM results
clear; close all;

x =linspace(-1,1,4);
y =linspace(-1,1,4);

[xv,yv] = meshgrid(x,y);
sM.coord(:,1) = xv(:);
sM.coord(:,2) = yv(:);
sM.connec = delaunay(sM.coord);

sM.coord = [0 0;1 0;0 1];
sM.connec = [1 2 3];

% sM.coord = [0 0;1 0;0 1];
% sM.connec = [2 3 1];

% sM.coord = [1 1;0 1;1 0];
% sM.connec = [1 2 3];

m = Mesh(sM);
mesh = m;

%% Create functions
% AnalyticalFunction


sAF.fHandle = @(x,y) x(1,:,:);
% sAF.fHandle = @(x,y) [-x(2,:,:)./10,x(1,:,:)./10];
sAF.ndimf   = 1;
sAF.mesh    = UnitQuadMesh(5,5);
xFun = AnalyticalFunction(sAF);

p0 = xFun.project('P0');
p1 = xFun.project('P1');
p2 = xFun.project('P2');

% %% Create projectors to P0, P1 and P1D
% % Projector to P0
% pp0.mesh   = mesh;
% pp0.connec = mesh.connec;
% projP0 = Projector_toP0(pp0);
% p0fun = projP0.project(xFun);
% p0fun.plot()
% title('P0')
% 
% % Projector to P1
% pp1.mesh   = mesh;
% pp1.connec = mesh.connec;
% projP1 = Projector_toP1(pp1);
% p1fun = projP1.project(xFun);
% p1fun.plot()
% title('P1 (quad linear)')
% 
% % Projector to P1 Discontinuous
% pp1d.mesh   = mesh;
% pp1d.connec = mesh.connec;
% projP1D = Projector_toP1Discontinuous(pp1d);
% p1dfun = projP1D.project(xFun);
% p1dfun.plot()
% title('P1 disc (quad quadratic)')

