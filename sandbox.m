%% This is a sandbox file!
% Feel free to test anything here :)
clc; clear; close all;

% file = 'test2d_triangle';
% file = 'test2d_quad';
file = 'test3d_hexahedra';
a.fileName = file;
s = FemDataContainer(a);
mesh = s.mesh;
fem = FEM.create(s);
fem.solve();

%% Create functions
% AnalyticalFunction

% sAF.fHandle = @(x) x(1,:,:);
% sAF.ndimf   = 1;
sAF.fHandle = @(x) [x(1,:,:).^2; x(2,:,:)];
sAF.ndimf   = 2;
sAF.mesh    = mesh;
xFun = AnalyticalFunction(sAF);

%% Create projectors to P0, P1 and P1D
% testingFunctions.m
% testingGradients.m

% Projector to P1
ppar.mesh   = mesh;
ppar.connec = mesh.connec;
projP1 = Projector_toP1(ppar);
p1fun = projP1.project(xFun);
% p1fun.plot(mesh)
% title('P1 (quad linear)')

% Projector to P0
projP0 = Projector_toP0(ppar);
p0fun = projP0.project(xFun);

% Projector to P1 Discontinuous
projP1D = Projector_toP1Discontinuous(ppar);
p1dfun = projP1D.project(xFun);

% FGaussDiscontinuousFunction
quad = Quadrature.set(mesh.type);
quad.computeQuadrature('LINEAR');
fgfun = p1fun.computeGradient(quad,mesh);
% fgp1 = projP1.project(fgfun);

%% Function printing
% aa.mesh = mesh;
% aa.filename = 'p1fun';
% p1fun.print(aa)
% 
% aa.filename = 'p0fun';
% p0fun.print(aa)
% 
% aa.filename = 'p1dfun';
% p1dfun.print(aa)
% 
% aa.filename = 'fgfun';
% fgfun.print(aa)
% 
% aa.filename = 'fgp1fun';
% fgp1.print(aa)

%% Multiple function printing

% bb.mesh     = mesh;
% bb.filename = 'funfunfun';
% bb.fun      = {p0fun, p1fun};
% bb.funNames = {'p0', 'p1'};
% fp = FunctionPrinter(bb);
% fp.print();

%% Paraview
zz.mesh     = mesh;
zz.filename = 'paraview2';
zz.fun      = {fgfun, p1fun, p0fun};
zz.funNames = {'fgfun', 'p1fun', 'p0fun'};
pvPst = ParaviewPostprocessor(zz);
% pvPst = ParaviewLegacyPostprocessor(zz);
