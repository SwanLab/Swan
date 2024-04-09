clear
close all

% % INPUT DATA
% file = 'test_gerard';
% a.fileName = file;
% f = StokesDataContainer(a);



m = QuadMesh(3,1.5,100,100); %Per alguna raó, amb 150 la P no surt correcte.
s.type='CircleInclusion';%
s.radius = 0.25;
s.xCoorCenter = 0.5;
s.yCoorCenter = 0.5;
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
isCyl    = @(coor) (abs(coor(:,1) - 0.5).^2+abs(coor(:,2) - 0.5).^2-0.25^2 < 1e-5);

%% Original (no-slip condition)
dir_vel{2}.domain    = @(coor) isTop(coor) | isBottom(coor) | isCyl(coor);
dir_vel{2}.direction = [1,2];
dir_vel{2}.value     = [0,0]; 

dir_vel{1}.domain    = @(coor) isLeft(coor) & not(isTop(coor) | isBottom(coor));
dir_vel{1}.direction = [1,2];
dir_vel{1}.value     = [1,0];

%% Cavity:
% dir_vel{2}.domain    = @(coor) isTop(coor) | isBottom(coor) | isCyl(coor) | isRight;
% dir_vel{2}.direction = [1,2];
% dir_vel{2}.value     = [0,0]; 
% 
% dir_vel{1}.domain    = @(coor) isLeft(coor) & not(isTop(coor) | isBottom(coor));
% dir_vel{1}.direction = [1,2];
% dir_vel{1}.value     = [0,1];

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


%% 
% dir_vel{2}.domain    = @(coor) isLeft(coor) | isRight(coor) | isCyl(coor);
% dir_vel{2}.direction = [1,2];
% dir_vel{2}.value     = [0,0]; 

% dir_vel{2}.domain    = @(coor) isLeft(coor) | isRight(coor);
% dir_vel{2}.direction = [1,2];
% dir_vel{2}.value     = [0,0]; 
% 

% 
% dir_vel{1}.domain    = @(coor) isTop(coor) & not(isLeft(coor) | isRight(coor)); %Hem de posar el not per no incloure els nodes a les puntes
% dir_vel{1}.direction = [1,2];
% dir_vel{1}.value     = [0,1];
% 
% dir_pre{1}.domain    = @(coor) isLeft(coor) & isTop(coor);
% dir_pre{1}.direction = 1;
% dir_pre{1}.value     = 0;

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
