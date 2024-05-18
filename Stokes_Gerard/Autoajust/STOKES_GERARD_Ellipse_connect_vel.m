clear all
close all

% % INPUT DATA

dim_a = 0.25; % Semi-major axis 0.2
dim_b = 0.24; % Semi-minor axis 0.02
center_posx = 0.5; % x position of the ellipse center
center_posy = 0.5; % y position of the ellipse center
AOAd = 0; % Angle of attack of the semi-major axis (in degrees)

metode = 1; %Mètode del marge


m = QuadMesh(1,1,100,100); % MESH
s.type='Given';
AOAr = -deg2rad(AOAd);

ellipse = calc_ellipse_classe(AOAr,center_posx,center_posy);
del_ab = ellipse.solvesys();

% del_ab=calc_ellipse(AOAr,center_posx,center_posy);

s.fHandle = @(x) -((((x(1,:,:)*cos(AOAr)+x(2,:,:)*sin(AOAr))-del_ab(1))/dim_a).^2+(((-x(1,:,:)*sin(AOAr)+x(2,:,:)*cos(AOAr))-del_ab(2))/dim_b).^2-1);
g = GeometricalFunction(s);
lsFun = g.computeLevelSetFunction(m); %D'aquí surt la malla de quadrats sense el forat
sUm.backgroundMesh = m;
sUm.boundaryMesh = m.createBoundaryMesh(); %sUm.boundaryMesh conté les mesh de les quatre fronteres del voltant. No té res del forat
uMesh = UnfittedMesh(sUm);
uMesh.compute(lsFun.fValues); % uMesh.boundaryCutMesh.mesh  és el forat
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

correct_margin = false;
connectat = true;

if metode == 1
%Franja alta:
k=1;
h=8;
margin = k*(10^(-h));

while correct_margin == false
    connectat = true;

    isCyl    = @(coor) (abs(abs(((coor(:,1)*cos(AOAr)+coor(:,2)*sin(AOAr))-del_ab(1))/dim_a).^2 + abs(((-coor(:,1)*sin(AOAr)+coor(:,2)*cos(AOAr))-del_ab(2))/dim_b).^2 - 1) < margin);

    nodesCyl    = pressureFun.getDofsFromCondition(isCyl);
%     xCyl        = mesh.coord(nodesCyl,1);
%     yCyl        = mesh.coord(nodesCyl,2);
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
    
    %Check that the connectivieties are correct (there is no more than 2 repetitions of each value):
    single_col = reshape(bMesh.connec, [], 1);
    [unique_el, ~, idx] = unique(single_col);
    count = accumarray(idx(:), 1);
    repeated_num = unique_el(count > 1);
    repetitions = zeros(size(unique_el));
    for i = 1:length(unique_el)
       idx = find(single_col == unique_el(i));
       if numel(idx)>2.5
          connectat = false;
       end
       repetitions(i) = numel(idx);
    end
    
    %Check the "interpolation" nodes of the surface of the object
    Nodoccult = @(coor) isCyl(coor);
    Dofscyl = velocityFun.getDofsFromCondition(Nodoccult);
    disp(size(Dofscyl,1) - bMesh.nnodes*4);

    if size(bMesh.coord,1) == size(bMesh.connec,1) && size(bMesh.coord,1)~=0 && connectat==true && size(Dofscyl,1) == bMesh.nnodes*4
        marginbo=margin;
        correct_margin = true;
        break
    else
        if k>9
            h=h-1;
            k=1;
        else
            k=k+0.1;
        end
        margin = k*(10^(-h));
    end

    disp(margin)


    if margin>1
        plot(bMesh);
        disp('MARGE NO TROBAT')
        return
    end


end

else

%Franja baixa:
k=9;
h=1;
margin = k*(10^-h);

while correct_margin == false
    connectat = true;
    isCyl    = @(coor) (abs(abs(((coor(:,1)*cos(AOAr)+coor(:,2)*sin(AOAr))-del_ab(1))/dim_a).^2 + abs(((-coor(:,1)*sin(AOAr)+coor(:,2)*cos(AOAr))-del_ab(2))/dim_b).^2 - 1) < margin);

    nodesCyl    = pressureFun.getDofsFromCondition(isCyl);
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
    
    %Check that the connectivieties are correct (there is no more than 2 repetitions of each value):
    single_col = reshape(bMesh.connec, [], 1);
    [unique_el, ~, idx] = unique(single_col);
    count = accumarray(idx(:), 1);
    repeated_num = unique_el(count > 1);
    repetitions = zeros(size(unique_el));
    for i = 1:length(unique_el)
       idx = find(single_col == unique_el(i));
       if numel(idx)>2.5
          connectat = false;
       end
       repetitions(i) = numel(idx);
    end


    if size(bMesh.coord,1)==0 && size(bMesh.connec,1)==0
        disp('MARGE NO TROBAT')
        return

    end

    
    %Check the "interpolation" nodes of the surface of the object
    Nodoccult = @(coor) isCyl(coor);
    Dofscyl = velocityFun.getDofsFromCondition(Nodoccult);
    
    if size(bMesh.coord,1) == size(bMesh.connec,1) && size(bMesh.coord,1)~=0 && connectat==true && size(Dofscyl,1) == bMesh.nnodes*4
        marginbo=margin;
        correct_margin = true; 
        break
    else
        if k<1
            h=h+1;
            k=9;
        else
            k=k-0.1;
        end
        margin = k*(10^(-h));
    end

    disp(margin)

end

end

disp(margin)
plot(bMesh);
disp(size(bMesh.coord,1) - size(bMesh.connec,1));

isCyl    = @(coor) (abs(abs(((coor(:,1)*cos(AOAr)+coor(:,2)*sin(AOAr))-del_ab(1))/dim_a).^2 + abs(((-coor(:,1)*sin(AOAr)+coor(:,2)*cos(AOAr))-del_ab(2))/dim_b).^2 - 1) < marginbo);

%% Original (no-slip condition)
dir_vel{2}.domain    = @(coor)isTop(coor) | isBottom(coor) |  isCyl(coor); 
dir_vel{2}.direction = [1,2];
dir_vel{2}.value     = [0,0]; 

dir_vel{1}.domain    = @(coor) isLeft(coor) & not(isTop(coor) | isBottom(coor));
dir_vel{1}.direction = [1,2];
dir_vel{1}.value     = [1,0];


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

% PLOT RESULTS
velocityFun.plot()
pressureFun.plot()

%% Lift and drag

nodesCyl    = pressureFun.getDofsFromCondition(isCyl); %Abans (per dirichlet) hem fet el mateix però pels nodes de la velocitat (velocityFun)
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

