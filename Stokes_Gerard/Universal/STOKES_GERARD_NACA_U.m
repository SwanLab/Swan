clear
close all

%% Input data
% Backgorund mesh
m = QuadMesh(10,4,150,150*0.8);
s.type='Given';

% Ellipse
% dim_a = 0.5; % Semi-major axis
% dim_b = 0.08; % Semi-minor axis

% Airfoil NACA4
M=5/100;
p=5/10;
%t=12/100;

% Beam (inside the airfoil):
alt = 0.11;
ampl = 0.095;
x_pos = 0.3; %Centre de la biga respecte el LE

% General parameters
AOAd = 7; %deg
x_centr = 3.5;
y_centr = 2;

%% Airfoil creation

% Ellipse
% AOAr = deg2rad(AOAd);
% fH = @(x) -((((cos(AOAr)*(x(1,:,:)-x_centr)-(x(2,:,:)-y_centr)*sin(AOAr)).^2)/(dim_a^2)) ...
%     + ((((x(2,:,:)-y_centr)*cos(AOAr)+(x(1,:,:)-x_centr)*sin(AOAr)).^2)/(dim_b^2)) - 1);

% Airfoil NACA4
[t,Q] = findthickness(x_pos,ampl,alt,M,p);
fH = Find_fH_circles(M,p,t,x_centr,y_centr,AOAd);
% fH = Find_fH_points(M,p,t,x_centr,y_centr,AOAd);

%% Create mesh

s.fHandle = fH;
g = GeometricalFunction(s);
lsFun = g.computeLevelSetFunction(m); %D'aquí surt la malla de quadrats sense el forat
sUm.backgroundMesh = m;
sUm.boundaryMesh = m.createBoundaryMesh(); %sUm.boundaryMesh conté les mesh de les quatre fronteres del voltant. No té res del forat
uMesh = UnfittedMesh(sUm);
uMesh.compute(lsFun.fValues); % uMesh.boundaryCutMesh.mesh  és el forat
mesh = uMesh.createInnerMesh();
figure
plot(uMesh)
% hold on
% scatter(punts_rot(1,:,1),punts_rot(2,:,1))
%figure
%plot(lsFun)
e.type  = 'STOKES';
e.nelem = mesh.nelem;
material = Material.create(e);
dtime = Inf; %Estacionari

% VELOCITY AND PRESSURE FUNCTIONS
velocityFun = LagrangianFunction.create(mesh, 2, 'P2');
pressureFun = LagrangianFunction.create(mesh, 1, 'P1');
n_dofs = velocityFun.nDofs + pressureFun.nDofs;

%% Boudary conditions

[forcesFormula,dirichlet,dir_dofs,nodespresscyl] = boundary_conditions(mesh,uMesh,velocityFun,pressureFun);


%% Solver

[velocityFun,pressureFun] = solver_stokesG(forcesFormula,dirichlet,dir_dofs,velocityFun,pressureFun,dtime,mesh,material,n_dofs);


%% PLOT RESULTS
velocityFun.plot()
pressureFun.plot()
% caxis([-50 50]);


%% Lift and drag

[L,D] = aero_forces(nodespresscyl,pressureFun,mesh);

