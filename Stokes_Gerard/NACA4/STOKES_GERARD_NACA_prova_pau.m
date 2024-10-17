clear all
close all

% % INPUT DATA

m = QuadMesh(6,4,200,100); % MESH
s.type='Given';

% NACA 4
M=9/100;
p=4/10;
t=12/100;

AOAd = 30; %deg
x_centr = 2;
y_centr = 2;

%% Airfoil creation

pas=0.0001;

x_p=[0:pas:1];

yt = 5*t*(0.2969*sqrt(x_p)-0.1260*x_p-0.3516*x_p.^2+0.2843*x_p.^3-0.1036*x_p.^4); %-0.1015 el valor que hi ha ara és per tenir la punta tancada

for j=1:1:size(x_p,2)
    if x_p(j)<=p
        y_c(j)=(M/(p^2))*(2*p*x_p(j)-x_p(j)^2);
        theta(j)=atan((2*M/(p^2))*(p-x_p(j)));
    elseif x_p(j)>p
        y_c(j)=(M/(1-p)^2)*((1-2*p)+2*p*x_p(j)-x_p(j)^2);
        theta(j)=atan(((2*M/((1-p)^2))*(p-x_p(j))));
    end
end

x_u = x_p - yt.*sin(theta);
y_u = y_c + yt.*cos(theta);
x_l = x_p + yt.*sin(theta);
y_l = y_c - yt.*cos(theta);

x_coord = [x_u flip(x_l)];
y_coord = [y_u flip(y_l)];

plot(x_coord,y_coord);
hold on
plot(x_l,y_l);
axis equal


x_LE=x_centr-0.5;
AOA = -deg2rad(AOAd);
x_cnr = x_coord+x_LE;
y_cnr = y_coord+y_centr;

x_cn = (x_cnr-x_centr).*cos(AOA)-(y_cnr-y_centr).*sin(AOA)+x_centr;
y_cn = (x_cnr-x_centr).*sin(AOA)+(y_cnr-y_centr).*cos(AOA)+y_centr;

fH = @(x)  double(inpolygon(x(1,:,:),x(2,:,:),x_cn,y_cn))-0.5; 



%% Create mesh and boundary conditions
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
figure
plot(lsFun)
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
dir_vel{1}.value     = [1,0]; %Velocity on the inlet

%Nodesnormals = uMesh.boundaryCutMesh.mesh

%Trobem els nodes de "pressió" al voltant de l'el·lipse (nodes que no són
%intermitjos). Recordar que estem buscant els GDL de la velocitat per
%imposar les condicions de contorn, no els de la pressió (els noms de les
%variables poden confondre)
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

% Pressure bc
isHorizontal = @(coor) (abs(coor(:,2) - 2) < 1e-12); % La pressió es fixa a 

dir_pre{1}.domain    = @(coor) isRight(coor) & isHorizontal(coor);
dir_pre{1}.direction = 1;
dir_pre{1}.value     = 0;

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
caxis([-50 50]);
%caxis([-115 80]);

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

cont =1;

for iE = 1:bMesh.nelem
    node1 = bMesh.coord(bMesh.connec(iE,1),:);
    node2 = bMesh.coord(bMesh.connec(iE,2),:);

    if node1(1)<= 5
    nvect = (node2-node1)/(abs(norm(node2-node1)));
    nvect = -nvect * [0 -1;1 0];
%     if dot(ref_vect(iE,:),nvect)<0 %No cal
%         nvect = -nvect;
%     end
    normal_vectors(cont,:) = nvect;
    length_element(cont) = abs(norm(node1-node2));

    cont = cont +1;
    end

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

