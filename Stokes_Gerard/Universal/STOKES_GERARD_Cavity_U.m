clear
close all

%% Input data and mesh

mesh = QuadMesh(1,1,54,54); 

e.type  = 'STOKES';
e.nelem = mesh.nelem;
material = Material.create(e);
dtime = Inf;

% VELOCITY AND PRESSURE FUNCTIONS
velocityFun = LagrangianFunction.create(mesh, 2, 'P2');
pressureFun = LagrangianFunction.create(mesh, 1, 'P1');
n_dofs = velocityFun.nDofs + pressureFun.nDofs;

%% Bonudary conditions
isLeft   = @(coor) (abs(coor(:,1) - min(coor(:,1)))   < 1e-12);
isRight  = @(coor) (abs(coor(:,1) - max(coor(:,1)))   < 1e-12);
isBottom = @(coor) (abs(coor(:,2) - min(coor(:,2)))   < 1e-12);
isTop    = @(coor) (abs(coor(:,2) - max(coor(:,2)))   < 1e-12);

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

[velocityFun,pressureFun] = solver_stokesG(forcesFormula,dirichlet,dir_dofs,velocityFun,pressureFun,dtime,mesh,material,n_dofs);


%% PLOT RESULTS
velocityFun.plot()
pressureFun.plot()
% caxis([-3 3]);

% Find pressure and velocity in a vertical line
x_line = 0.5;
tol = 1e-12;
isVerticalLine = @(coor) (abs(coor(:,1) - x_line) < tol);

dirDofs_vel = velocityFun.getDofsFromCondition(isVerticalLine);
dirDofs_pre = pressureFun.getDofsFromCondition(isVerticalLine);
node_vel = 1 + (dirDofs_vel(2:2:end)-2)/velocityFun.ndimf;

xnod  = mesh.coord(dirDofs_pre,1);
ynod  = mesh.coord(dirDofs_pre,2);
vels = velocityFun.fValues(dirDofs_pre,:);
press = pressureFun.fValues(dirDofs_pre);

figure
plot(vels(:,1),ynod,'LineWidth',1)
grid on
xlabel('Velocity in x direction');
ylabel('y position')

figure
plot(press,ynod,'LineWidth',1.1)
axis([0.3 0.4 0 1])
grid on
xlabel('Pressure');
ylabel('y position')



















