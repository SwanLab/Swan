clear
close all

%Input data
H=1;

% for S = 50:50:1000
S=54
mesh = QuadMesh(1,1,S,S); 

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


%% Original (no-slip condition)
dir_vel{2}.domain    = @(coor) isTop(coor);
dir_vel{2}.direction = [1,2];
dir_vel{2}.value     = [1,0]; 

dir_vel{1}.domain    = @(coor) (isLeft(coor) | isRight(coor) | isBottom(coor)) & not(isTop(coor));
dir_vel{1}.direction = [1,2];
dir_vel{1}.value     = [0,0];

dir_pre{1}.domain    = @(coor) isBottom(coor) & isLeft(coor);
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
for i = 1:length(dir_pre)
    dirDofs = pressureFun.getDofsFromCondition(dir_pre{i}.domain);
    mat12 = ones(size(dirDofs));
    valmat = ones(size(dirDofs)).*dir_pre{i}.value';
    dirichlet(size(dirichlet,1)+1:size(dirichlet,1)+length(dirDofs),:) = [dirDofs+velocityFun.nDofs mat12 valmat];
    dir_dofs(size(dir_dofs,1)+1:size(dir_dofs,1)+length(dirDofs),1) = dirDofs+velocityFun.nDofs;
end


% DEFINE APPLIED FORCES
sAF.fHandle = @(coor) [0.*coor,0.*coor];
sAF.ndimf   = 2;
sAF.mesh    = mesh;
forcesFormula = AnalyticalFunction(sAF);

%% Solver
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


% COMPUTE RHS
d.type          = 'Stokes';
d.mesh          = mesh;
d.velocityFun   = velocityFun;
d.pressureFun   = pressureFun;
d.forcesFormula = forcesFormula;
RHSint = RHSintegrator.create(d);
F = RHSint.integrate();
uD = dirichlet(:,3);
R  = -LHS(:,dir_dofs)*uD;
RHS = F + R;

% SOLVE PROBLEM
free_dofs_plus = setdiff(1:n_dofs,dir_dofs);
LHSr = LHS(free_dofs_plus,free_dofs_plus); 
RHSr = RHS(free_dofs_plus);
x = solver.solve(LHSr, RHSr);


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

%% PLOT RESULTS
velocityFun.plot()
pressureFun.plot()
% caxis([-3 3]);

x_target = 0.2;
y_target = 0.2;
x_line = 0.5;
tol = 1e-12;


isNode = @(coor) (sqrt((coor(:,1) - x_target).^2 + (coor(:,2) - y_target).^2) < tol);
isVerticalLine = @(coor) (abs(coor(:,1) - x_line) < tol);

dirDofs_vel = velocityFun.getDofsFromCondition(isVerticalLine);
dirDofs_pre = pressureFun.getDofsFromCondition(isVerticalLine);
node_vel = 1 + (dirDofs_vel(2:2:end)-2)/velocityFun.ndimf;

figure
xnod        = mesh.coord(dirDofs_pre,1);
ynod        = mesh.coord(dirDofs_pre,2);
scatter(xnod,ynod)
vels = velocityFun.fValues(dirDofs_pre,:); %Busquem la velocitat només en els nodes que comparteix amb la pressió
press = pressureFun.fValues(dirDofs_pre);

plot(vels(:,1),ynod,'LineWidth',1.1)
load('153x153_codiweb.mat')
hold on
plot(vel(:,1),xynodv(:,2),'r','LineWidth',1.1)
% legend('Swan code','LaCàN code','FontSize', 15)
% set(gca, 'FontSize', 15)
% ylabel('y position','FontSize', 15)
% xlabel('Velocity in x direction','FontSize', 15)
grid on

% figure
% plot(vels(:,2),ynod)
% hold on
% plot(vel(:,2),xynodv(:,2),'r')
% legend('Swan code','LaCàN code','FontSize', 16)
% set(gca, 'FontSize', 16)
% ylabel('y position','FontSize', 16)
% xlabel('Velocity in y direction','FontSize', 16)

figure
plot(press,ynod,'LineWidth',1.1)
hold on
plot(pre(:,1),xynodp(:,2),'r','LineWidth',1.1)
% legend('Swan code','LaCàN code','FontSize', 15)
% set(gca, 'FontSize', 15)
% ylabel('y position','FontSize', 15)
% xlabel('Pressure','FontSize', 15)
axis([0.3 0.4 0 1])
grid on

% H=H+1;
% clearvars('-except', 'vel', 'press','H');
% end
% % scatter3(velocityFun.coord(nodes(:),1),velocityFun.coord(nodes(:),2),1000,'r','x');



















