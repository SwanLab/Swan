clear all
close all
clc

% % INPUT DATA
% file = 'test_gerard';
% a.fileName = file;
% f = StokesDataContainer(a);

dim_a = 0.2; % Semi-major axis 0.2
dim_b = 0.05; % Semi-minor axis 0.02
center_posx = 0.7; % x position of the ellipse center
center_posy = 0.5; % y position of the ellipse center
AOAd = -20; % Angle of attack of the semi-major axis (in degrees)


m = QuadMesh(2,1,100,100); % MESH
s.type='Given';
AOAr = -deg2rad(AOAd);

ellipse = calc_ellipse_classe(AOAr,center_posx,center_posy);
del_ab = ellipse.solvesys();

% del_ab=calc_ellipse(AOAr,center_posx,center_posy);

s.fHandle = @(x) -((((x(1,:,:)*cos(AOAr)+x(2,:,:)*sin(AOAr))-del_ab(1))/dim_a).^2+(((-x(1,:,:)*sin(AOAr)+x(2,:,:)*cos(AOAr))-del_ab(2))/dim_b).^2-1);
g = GeometricalFunction(s);
lsFun = g.computeLevelSetFunction(m);
sUm.backgroundMesh = m;
sUm.boundaryMesh = m.createBoundaryMesh();
uMesh = UnfittedMesh(sUm);
uMesh.compute(lsFun.fValues);
mesh = uMesh.createInnerMesh();

% mesh = TriangleMesh(1,1,40,40);

e.type  = 'STOKES';
e.nelem = mesh.nelem;
material = Material.create(e);
dtime = Inf;

% VELOCITY AND PRESSURE FUNCTIONS
velocityFun = LagrangianFunction.create(mesh, 2, 'P2');
pressureFun = LagrangianFunction.create(mesh, 1, 'P1');
n_dofs = velocityFun.nDofs + pressureFun.nDofs;

% DEFINE BOUNDARY CONDITIONS
isLeft   = @(coor) (abs(coor(:,1) - min(coor(:,1)))   < 1e-12);
isRight  = @(coor) (abs(coor(:,1) - max(coor(:,1)))   < 1e-12);
isBottom = @(coor) (abs(coor(:,2) - min(coor(:,2)))   < 1e-12);
isTop    = @(coor) (abs(coor(:,2) - max(coor(:,2)))   < 1e-12);
isCyl    = @(coor) (abs(abs(((coor(:,1)*cos(AOAr)+coor(:,2)*sin(AOAr))-del_ab(1))/dim_a).^2 + abs(((-coor(:,1)*sin(AOAr)+coor(:,2)*cos(AOAr))-del_ab(2))/dim_b).^2 - 1) < 3.5e-3); %3.5e-2 per l'el·lipse fina


%El problema és que li costa distingir els nodes a la frontera, perquè en
%alguns llocs si que els nodes cauen a prop de l'el·lipse, però en altres
%no.
% %Distribució
% isMiddle_a1 = @(coor) (abs(coor(:,2)) - 0.9   < 1e-12);
% isMiddle_a2 = @(coor) (0.1 - abs(coor(:,2))   < 1e-12);
% 
% isMiddle_b1 = @(coor) (abs(coor(:,2)) - 0.8   < 1e-12);
% isMiddle_b2 = @(coor) (0.2 - abs(coor(:,2))   < 1e-12);
% 
% isMiddle_c1 = @(coor) (abs(coor(:,2)) - 0.7   < 1e-12);
% isMiddle_c2 = @(coor) (0.3 - abs(coor(:,2))   < 1e-12);
% 
% isMiddle_d1 = @(coor) (abs(coor(:,2)) - 0.95   < 1e-12);
% isMiddle_d2 = @(coor) (0.05 - abs(coor(:,2))   < 1e-12);

%% Original (no-slip condition)
dir_vel{2}.domain    = @(coor) isTop(coor) | isBottom(coor) | isCyl(coor);
dir_vel{2}.direction = [1,2];
dir_vel{2}.value     = [0,0]; 

dir_vel{1}.domain    = @(coor) isLeft(coor) & not(isTop(coor) | isBottom(coor));
dir_vel{1}.direction = [1,2];
dir_vel{1}.value     = [1,0];

%% Prova distribució
% dir_vel{2}.domain    = @(coor) isTop(coor) | isBottom(coor) | isCyl(coor);
% dir_vel{2}.direction = [1,2];
% dir_vel{2}.value     = [0,0]; 
% 
% dir_vel{1}.domain    = @(coor) isLeft(coor) & not(isTop(coor) | isBottom(coor));
% dir_vel{1}.direction = [1,2];
% dir_vel{1}.value     = [0.2,0];
% 
% dir_vel{3}.domain    = @(coor) isMiddle_d1(coor) & isMiddle_d2(coor) & isLeft(coor);
% dir_vel{3}.direction = [1,2];
% dir_vel{3}.value     = [0.4,0];
% 
% dir_vel{4}.domain    = @(coor) isMiddle_a1(coor) & isMiddle_a2(coor) & isLeft(coor);
% dir_vel{4}.direction = [1,2];
% dir_vel{4}.value     = [0.6,0];
% 
% dir_vel{5}.domain    = @(coor) isMiddle_b1(coor) & isMiddle_b2(coor) & isLeft(coor);
% dir_vel{5}.direction = [1,2];
% dir_vel{5}.value     = [0.8,0];
% 
% dir_vel{6}.domain    = @(coor) isMiddle_c1(coor) & isMiddle_c2(coor) & isLeft(coor);
% dir_vel{6}.direction = [1,2];
% dir_vel{6}.value     = [1,0];

%% Modificat (free-slip condition)
% dir_vel{2}.domain    = @(coor) isTop(coor) | isBottom(coor);
% dir_vel{2}.direction = [1,2];
% dir_vel{2}.value     = [1,0]; 
% 
% dir_vel{1}.domain    = @(coor) isLeft(coor) & not(isTop(coor) | isBottom(coor));
% dir_vel{1}.direction = [1,2];
% dir_vel{1}.value     = [1,0];
% 
% dir_vel{3}.domain    = @(coor) isCyl(coor);
% dir_vel{3}.direction = [1,2];
% dir_vel{3}.value     = [0,0]; 

%% Altres
% dir_pre{1}.domain    = @(coor) isLeft(coor) & isTop(coor);
% dir_pre{1}.direction = 1;
% dir_pre{1}.value     = 0;

% dir_vel{2}.domain    = @(coor) isTop(coor) & not(isLeft(coor));
% dir_vel{2}.direction = [1,2];
% dir_vel{2}.value     = [0,0]; 
% 
% dir_vel{4}.domain    = @(coor) isBottom(coor) & not(isLeft(coor));
% dir_vel{4}.direction = [1,2];
% dir_vel{4}.value     = [0,0]; 
% 
% dir_vel{3}.domain    = @(coor) isCyl(coor);
% dir_vel{3}.direction = [1,2];
% dir_vel{3}.value     = [0,0]; 
% 
% dir_vel{1}.domain    = @(coor) isLeft(coor);
% dir_vel{1}.direction = [1,2];
% dir_vel{1}.value     = [1,0];
% 
dir_pre{1}.domain    = @(coor) isLeft(coor) & isTop(coor);
dir_pre{1}.direction = 1;
dir_pre{1}.value     = 0;

dirichlet = [];
dir_dofs = [];
for i = 1:length(dir_vel)
    dirDofs = velocityFun.getDofsFromCondition(dir_vel{i}.domain);
    nodes = 1 + (dirDofs(2:2:end)-2)/velocityFun.ndimf;
    nodes2 = repmat(nodes, [1 2]);
    iNod = sort(nodes2(:));
    mat12 = repmat([1;2], [length(iNod)/2 1]);
    valmat = repmat(dir_vel{i}.value', [length(iNod)/2 1]);
    dirichlet(size(dirichlet,1)+1:size(dirichlet,1)+length(iNod),:) = [iNod mat12 valmat];
    dir_dofs(size(dir_dofs,1)+1:size(dir_dofs,1)+length(iNod),1) = dirDofs;
end
% for i = 1:length(dir_pre)
%     dirDofs = pressureFun.getDofsFromCondition(dir_pre{i}.domain);
%     mat12 = ones(size(dirDofs));
%     valmat = ones(size(dirDofs)).*dir_pre{i}.value';
%     dirichlet(size(dirichlet,1)+1:size(dirichlet,1)+length(dirDofs),:) = [dirDofs+velocityFun.nDofs mat12 valmat];
%     dir_dofs(size(dir_dofs,1)+1:size(dir_dofs,1)+length(dirDofs),1) = dirDofs+velocityFun.nDofs;
% end

% DEFINE APPLIED FORCES
sAF.fHandle = @(coor) [0.*coor,0.*coor];
sAF.ndimf   = 2;
sAF.mesh    = mesh;
forcesFormula = AnalyticalFunction(sAF);

% CREATE SOLVER
b.type =  'DIRECT';
solver = Solver.create(b);

% COMPUTE LHS
c.type          = 'Stokes';
c.dt            = dtime;
c.mesh          = mesh;
c.material      = material;
c.velocityFun   = velocityFun;
c.pressureFun   = pressureFun;
LHSintegrator = LHSintegrator.create(c);
LHS = LHSintegrator.compute();
% LHS(end+1,:) = [0;ones(n_dofs-1,1)]';
% LHS(1:end,end+1) = [0;ones(n_dofs-1,1);0];

% COMPUTE RHS
d.type          = 'Stokes';
d.mesh          = mesh;
d.velocityFun   = velocityFun;
d.pressureFun   = pressureFun;
d.forcesFormula = forcesFormula;
RHSint = RHSintegrator.create(d);
F = RHSint.integrate();
% F(end+1) = 0;
uD = dirichlet(:,3);
R  = -LHS(:,dir_dofs)*uD;
RHS = F + R;

% SOLVE PROBLEM
% free_dofs_plus = setdiff(1:(n_dofs+1),dir_dofs);
free_dofs_plus = setdiff(1:n_dofs,dir_dofs);
LHSr = LHS(free_dofs_plus,free_dofs_plus); %Li treiem els nodes restringits per deixar la LHS només amb lliures i la RHS de la mateixa mida
RHSr = RHS(free_dofs_plus);
x = solver.solve(LHSr, RHSr);
% x(end)=[];

% ADD DDIRICHLET BOUNDARY CONDITIONS
uD  = dirichlet(:,3);
nsteps = length(x(1,:));
uD = repmat(uD,1,nsteps);
fullx = zeros(n_dofs,nsteps);
free_dofs = setdiff(1:(n_dofs),dir_dofs);
fullx(free_dofs,:) = x;
if ~isempty(dir_dofs)
    fullx(dir_dofs,:) = uD;
end

% SEPARATE VARIABLES FVALUES
ndofsV = velocityFun.nDofs;
vars.u = fullx(1:ndofsV,:);
vars.p = fullx(ndofsV+1:end,:);

% DEFINE VARIABLES
nu = velocityFun.ndimf;
nnode = round(length(vars.u)/nu);
nodes = 1:nnode;
velfval = zeros(nnode,nu);
for idim = 1:nu
    dofs = nu*(nodes-1)+idim;
    velfval(:,idim) = vars.u(dofs, end);
end
velocityFun.fValues = velfval;
pressureFun.fValues = vars.p(:,end);

% PLOT RESULTS
velocityFun.plot()
pressureFun.plot()

%% PRESURE AT SURFACE
boundary_mesh = uMesh.boundaryCutMesh.mesh;

pressure_boundary = uMesh.obtainFunctionAtUnfittedMesh(pressureFun);
%velocity_boundary = uMesh.obtainFunctionAtUnfittedMesh(velocityFun);
pressure_boundary.boundaryCutMeshFunction.plot() %Plot the ellipse

normal_vectors = zeros(boundary_mesh.nelem,boundary_mesh.ndim); %Assign the space for the vectors
length_element = zeros(boundary_mesh.nelem,1); %Assign the space for lenghts of the elements

centroid = mean(boundary_mesh.coord);
central_points = (boundary_mesh.coord(boundary_mesh.connec(:,1),:)+boundary_mesh.coord(boundary_mesh.connec(:,2),:))/2;
ref_vect = central_points - centroid;

for iE = 1:boundary_mesh.nelem
    node1 = boundary_mesh.coord(boundary_mesh.connec(iE,1),:);
    node2 = boundary_mesh.coord(boundary_mesh.connec(iE,2),:);
    nvect = (node2-node1)/(abs(norm(node2-node1)));
    nvect = nvect * [0 -1;1 0];
    if dot(ref_vect(iE,:),nvect)<0
        nvect = -nvect;
    end
    normal_vectors(iE,:) = nvect;
    length_element(iE) = abs(norm(node1-node2));
end

F_total = [0,0];
for iE = 1:boundary_mesh.nelem
    pressure_mean = (pressure_boundary.boundaryCutMeshFunction.fValues(boundary_mesh.connec(iE,1)) + pressure_boundary.boundaryCutMeshFunction.fValues(boundary_mesh.connec(iE,2)))/2;
    F_total = F_total + normal_vectors(iE,:)*length_element(iE)*(-pressure_mean);
end

quiver(central_points(:,1),central_points(:,2),normal_vectors(:,1),normal_vectors(:,2)) %Plot the vectors
hold on
quiver(center_posx,center_posy,F_total(1,1),F_total(1,2));
hold on
boundary_mesh.plot() %Plot mesh points
