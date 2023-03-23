% Testing domain decomposition for reduced order modeling

clear
clc


% STEP 1: Compute full domain

% Create geometry
fileName = 'CantileverAxialLoad';

s.testName    = [fileName,'.m'];
s.problemCase = 'cantilever';
s.x1          = 0.5;
s.y1          = 0.5;
s.N           = 20;
s.M           = 10;
s.P           = 1;
s.DoF         = 2;

CantileverAxial = FEMInputWriter(s);
CantileverAxial.createTest();

% Create Mesh
[XPL,YPL] = meshgrid(0:0.05:0.5, 0:0.05:0.5);
[XPR,YPR] = meshgrid(0.5:0.05:1, 0:0.05:0.5);
X = [XPL(:);XPR(:)];
Y = [YPL(:);YPR(:)];
ss.coord(:,1) = X;
ss.coord(:,2) = Y;
ss.connec = delaunay(ss.coord);
mTot = Mesh(ss);
clear ss
ss.coord(:,1) = XPL(:);
ss.coord(:,2) = YPL(:);
ss.connec = delaunay(ss.coord);
mLeft = Mesh(ss);

clear ss
ss.coord(:,1) = XPR(:);
ss.coord(:,2) = YPR(:);
ss.connec = delaunay(ss.coord);
mRight = Mesh(ss);

% Create elasticity problem
a.fileName = fileName;
s = FemDataContainer(a);
s.mesh = mTot;
s.material.C = s.material.C(:,:,1:mTot.nelem);
fem = FEM.create(s);
fem.solve();
% fem.print(fileName);

% Get main results
u = fem.uFun{1,1};
% u.plot();
e = fem.strainFun{1,1};
% e.plot();
sig = fem.stressFun{1,1};
% sig.plot();

% ...


% STEP 2: Mesh decomposition
% ...
mesh = s.mesh;
xV = mesh.computeBaricenter();

x_baricenter = xV(1,:);
y_baricenter = xV(2,:);

isLeft = x_baricenter<= 0.5;
leftConnec = mesh.connec(isLeft,:);

isRight = x_baricenter> 0.5;
rightConnec = mesh.connec(isRight,:);


nodesLeft = unique(mesh.connec(isLeft,:));
nodesRight = unique(mesh.connec(isRight,:));

coordsLeft = mesh.coord(nodesLeft,:);

ss.connec = leftConnec;
ss.coord = mesh.coord;
leftMesh = Mesh(ss);

tt.connec = rightConnec;
tt.coord = mesh.coord;
rightMesh = Mesh(tt);

figure
leftMesh.plot
figure
rightMesh.plot


fValuesLeft = u.fValues(nodesLeft,:);
fValuesRight = u.fValues(nodesRight,:);

zz.fValues = fValuesLeft;
zz.mesh    = mLeft;
p1left = P1Function(zz);

zz.fValues = fValuesRight;
zz.mesh    = mRight;
p1right = P1Function(zz);

% STEP 2: BoundaryMesh
bMeshes = mLeft.createBoundaryMesh();
bMeshRight = bMeshes{2};
figure
bMeshRight.mesh.plot

% STEP 3: Stiffness matrix
% ...

sLHS.type     = 'ElasticStiffnessMatrix';
sLHS.mesh     = leftMesh;
sLHS.fun      = p1left;
sLHS.material = s.material;
lhsL = LHSintegrator.create(sLHS);
Kleft = lhsL.compute();

sLHS.mesh     = rightMesh;
sLHS.fun      = p1right;
sLHS.material = s.material;
lhsR = LHSintegrator.create(sLHS);
Kright = lhsR.compute();


a = 1;