clear 
close all

% Prova per veure si es pot trobar els nodes de la frontera de manera diferent.
% % INPUT DATA

m = QuadMesh(2,1,300,300); % MESH
s.type='Given';

% sx = 0.2;
% sy = 0.2;
% x0 = 0.5;
% y0 = 0.5;
% w  = 30;
% R  = [cosd(w),sind(w);-sind(w),cosd(w)];
% xLoc = @(x) abs(pagemtimes(R,[x(1,:,:)-x0;x(2,:,:)-y0]));
% fH = @(x) max(pagemtimes([1,0],xLoc(x))/sx,pagemtimes([0,1],xLoc(x))/sy) - 0.5;

r = [0.4 0.4 0.2 0.2];
x_p = [1.5 0.4 0.9 0.6];
y_c = [0.5 0.5 0.5 0.6];

%fH = @(x) max(((x(1,:,:)-1).^2 + 0.2 - x(2,:,:)),(-(x(1,:,:)-1).^2 + 0.8 - x(2,:,:)));

%fH = @(x) -((x(1,:,:)-1).^2 + 0.2 - x(2,:,:)).*(-(x(1,:,:)-1).^2 + 0.8 - x(2,:,:)).*(x(1,:,:)-0.45228);

fH = @(x) min((min((((1-x(1,:,:)).^2 + 0.2 - x(2,:,:)).*((1-x(1,:,:)).^2 - 0.8 + x(2,:,:))),(x(1,:,:)-0.46))),(-x(1,:,:)+1.547722));

%fH = @(x) -min(((x(1,:,:)-1).^2 + 0.2 - x(2,:,:)),(-(x(1,:,:)-1).^2 + 0.8 - x(2,:,:))).*min(((x(1,:,:)-1).^2 + 0.2 - x(2,:,:)),(-(x(1,:,:)-1).^2 + 0.8 - x(2,:,:)));
%fH = @(x) -((x(1,:,:)-x_p(1)).^2+(x(2,:,:)-y_c(1)).^2-r(1)^2).*((x(1,:,:)-x_p(2)).^2+(x(2,:,:)-y_c(2)).^2-r(2)^2).*((x(1,:,:)-x_p(3)).^2+(x(2,:,:)-y_c(3)).^2-r(3)^2);

%fH = @(x) -min([(((x(1,:,:)-x_p(1)).^2+(x(2,:,:)-y_c(1)).^2)-(r(1)^2)), (((x(1,:,:)-x_p(2)).^2+(x(2,:,:)-y_c(2)).^2)-(r(2)^2))]);
%fH = @(x) -min(((x(1,:,:)-x_p(1)).^2+(x(2,:,:)-y_c(1)).^2-r(1)^2),((x(1,:,:)-x_p(2)).^2+(x(2,:,:)-y_c(2)).^2-r(2)^2)) && -min(((x(1,:,:)-x_p(3)).^2+(x(2,:,:)-y_c(3)).^2-r(3)^2),((x(1,:,:)-x_p(4)).^2+(x(2,:,:)-y_c(4)).^2-r(4)^2));

%fH = @(x) -min(((x(1,:,:)-1).^2 + 0.2 - x(2,:,:)),(-(x(1,:,:)-1).^2 + 0.8 - x(2,:,:)));

%% Create mesh and boundary conditions
s.fHandle = fH;
g = GeometricalFunction(s);
lsFun = g.computeLevelSetFunction(m); %D'aquí surt la malla de quadrats sense el forat
plot(lsFun)
sUm.backgroundMesh = m;
sUm.boundaryMesh = m.createBoundaryMesh(); %sUm.boundaryMesh conté les mesh de les quatre fronteres del voltant. No té res del forat
uMesh = UnfittedMesh(sUm);
uMesh.compute(lsFun.fValues); % uMesh.boundaryCutMesh.mesh  és el forat
mesh = uMesh.createInnerMesh();

figure
plot(uMesh)

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


% Original (no-slip condition)
dir_vel{2}.domain    = @(coor) isTop(coor) | isBottom(coor);
dir_vel{2}.direction = [1,2];
dir_vel{2}.value     = [0,0]; 

dir_vel{1}.domain    = @(coor) isLeft(coor) & not(isTop(coor) | isBottom(coor));
dir_vel{1}.direction = [1,2];
dir_vel{1}.value     = [1,0];

%Nodesnormals = uMesh.boundaryCutMesh.mesh

%Trobem els nodes de pressió al voltant de l'el·lipse
size_cutmesh = size(uMesh.boundaryCutMesh.mesh.coord,1);
dirDofspresscyl=zeros(2,size_cutmesh);
for i = 1:1:size_cutmesh
    isxcoord    = @(coor) coor(:,1) == uMesh.boundaryCutMesh.mesh.coord(i,1);
    isycoord    = @(coor) coor(:,2) == uMesh.boundaryCutMesh.mesh.coord(i,2);
    dircyl      = @(coor) isxcoord(coor) & isycoord(coor);

    dirDofspresscyl(:,i) = velocityFun.getDofsFromCondition(dircyl);

end

plot(uMesh);
dirDofspresscyl_bo = sort(reshape(dirDofspresscyl,size(dirDofspresscyl,2)*2,1));
nodespresscyl = 1 + (dirDofspresscyl_bo(2:2:end)-2)/velocityFun.ndimf;
scatter(velocityFun.coord(nodespresscyl(:),1),velocityFun.coord(nodespresscyl(:),2),'X','b');

%Trobem les coordenades dels nodes intermitjos
coor_occult=zeros(size_cutmesh,2);
for i = 1:1:size_cutmesh
    coor_occult(i,1)=(uMesh.boundaryCutMesh.mesh.coord(uMesh.boundaryCutMesh.mesh.connec(i,1),1)+uMesh.boundaryCutMesh.mesh.coord(uMesh.boundaryCutMesh.mesh.connec(i,2),1))/2;
    coor_occult(i,2)=(uMesh.boundaryCutMesh.mesh.coord(uMesh.boundaryCutMesh.mesh.connec(i,1),2)+uMesh.boundaryCutMesh.mesh.coord(uMesh.boundaryCutMesh.mesh.connec(i,2),2))/2;
end

dirDofsoccucyl=zeros(2,size_cutmesh);
for i = 1:1:size_cutmesh
    isxcoord    = @(coor) coor(:,1) == coor_occult(i,1);
    isycoord    = @(coor) coor(:,2) == coor_occult(i,2);
    dircyloccu  = @(coor) isxcoord(coor) & isycoord(coor);

    dirDofsoccucyl(:,i) = velocityFun.getDofsFromCondition(dircyloccu);

end

dirDofsoccucyl_bo = sort(reshape(dirDofsoccucyl,size(dirDofsoccucyl,2)*2,1));
nodesoccucyl = 1 + (dirDofsoccucyl_bo(2:2:end)-2)/velocityFun.ndimf;
scatter(velocityFun.coord(nodesoccucyl(:),1),velocityFun.coord(nodesoccucyl(:),2),'o','g');


dirichlet = [];
dir_dofs = [];
for i = 1:1:4
    if i == 1 || i == 2
        dirDofs = velocityFun.getDofsFromCondition(dir_vel{i}.domain);
    elseif i == 3
        dirDofs = dirDofspresscyl_bo;
        dir_vel{i}.value     = [0,0]; 
    elseif i == 4
        dirDofs = dirDofsoccucyl_bo;
        dir_vel{i}.value     = [0,0]; 
    end
    
    nodes = 1 + (dirDofs(2:2:end)-2)/velocityFun.ndimf;
    scatter(velocityFun.coord(nodes(:),1),velocityFun.coord(nodes(:),2),'y','g');
    nodes2 = repmat(nodes, [1 2]);
    iNod = sort(nodes2(:));
    mat12 = repmat([1;2], [length(iNod)/2 1]);
    valmat = repmat(dir_vel{i}.value', [length(iNod)/2 1]);
    dirichlet(size(dirichlet,1)+1:size(dirichlet,1)+length(iNod),:) = [iNod mat12 valmat];
    dir_dofs(size(dir_dofs,1)+1:size(dir_dofs,1)+length(iNod),1) = dirDofs;
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
% F(end+1) = 0;
uD = dirichlet(:,3);
R  = -LHS(:,dir_dofs)*uD;
RHS = F + R;

% SOLVE PROBLEM
free_dofs_plus = setdiff(1:n_dofs,dir_dofs);
LHSr = LHS(free_dofs_plus,free_dofs_plus); %Li treiem els nodes restringits per deixar la LHS només amb lliures i la RHS de la mateixa mida
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

% isEsquerra   = @(coor) (abs(coor(:,1) - min(coor(:,1)))   < 1e-2);
% 
% nodes_elim = pressureFun.getDofsFromCondition(isEsquerra);
% pressureFun_cut = pressureFun;



%% Lift and drag

nodesCyl    = nodespresscyl; 
presCylVals = pressureFun.fValues(nodesCyl,1);
xCyl        = mesh.coord(nodesCyl,1);
yCyl        = mesh.coord(nodesCyl,2);
mesh.computeEdges();
e  = mesh.edges.nodesInEdges;
bE = ismember(e,nodesCyl);
bE = find(prod(bE,2));
connec = e(bE,:);
ss.coord    = mesh.coord;
ss.connec   = connec;
ss.kFace    = -1;
bMesh       = Mesh.create(ss);
bMesh       = bMesh.computeCanonicalMesh();
presCyl     = LagrangianFunction.create(bMesh,1,pressureFun.order); 
presCyl.fValues = presCylVals;

presCyl.plot()

normal_vectors = zeros(bMesh.nelem,bMesh.ndim);
length_element = zeros(bMesh.nelem,1);

centroid = mean(bMesh.coord);
central_points = (bMesh.coord(bMesh.connec(:,1),:)+bMesh.coord(bMesh.connec(:,2),:))/2;
ref_vect = central_points - centroid;

for iE = 1:bMesh.nelem
    node1 = bMesh.coord(bMesh.connec(iE,1),:);
    node2 = bMesh.coord(bMesh.connec(iE,2),:);
    nvect = (node2-node1)/(abs(norm(node2-node1)));
    nvect = nvect * [0 -1;1 0];
    if dot(ref_vect(iE,:),nvect)<0
        nvect = -nvect;
    end
    normal_vectors(iE,:) = nvect;
    length_element(iE) = abs(norm(node1-node2));
end

nx = LagrangianFunction.create(bMesh,1,'P0');%Vectors normals
ny = LagrangianFunction.create(bMesh,1,'P0');
nx.fValues = normal_vectors(:,1);
ny.fValues = normal_vectors(:,2);
sss.operation = @(x) -presCyl.evaluate(x).*nx.evaluate(x);
pNx           = DomainFunction(sss);
D             = Integrator.compute(pNx,bMesh,'QUADRATIC');
sss.operation = @(x) -presCyl.evaluate(x).*ny.evaluate(x);
pNy           = DomainFunction(sss);
L             = Integrator.compute(pNy,bMesh,'QUADRATIC');

quiver(central_points(:,1),central_points(:,2),normal_vectors(:,1),normal_vectors(:,2)) %Plot the vectors
hold on
quiver(centroid(1,1),centroid(1,2),D,0);
hold on
quiver(centroid(1,1),centroid(1,2),0,L);
hold on
bMesh.plot() %Plot mesh points

disp(L/D);