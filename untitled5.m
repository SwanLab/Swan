clc,clear,close all

fileName = 'CantileverAxialLoad';
s.testName    = [fileName,'.m'];
s.problemCase = 'cantilever';
s.x1          = 1;
s.y1          = 0.5;
s.N           = 799; %x dimension
s.M           = 399; %y dimension
s.P           = 1;
s.DoF         = 2;

CantileverAxial = FEMInputWriter(s);
CantileverAxial.createTest();

% Create elasticity problem
a.fileName = fileName;
data = FemDataContainer(a);
fem = FEM.create(data);
fem.solve();
fem.uFun.fValues = [fem.uFun.fValues zeros(size(fem.uFun.fValues,1),1)];
fem.uFun.ndimf = 3;
fem.print(fileName);

% Get main results
uTest = fem.uFun(1,1);
e_test = fem.strainFun(1,1);
sig_test = fem.stressFun(1,1);

%% Reactions plots fem

K = fem.LHS;
u = zeros(length(uTest.fValues)*2,1);
u(1:2:end-1) = uTest.fValues(:,1);
u(2:2:end) = uTest.fValues(:,2);
Reactions = K(fem.boundaryConditions.dirichlet,:)*u;
Reactions2plot = [Reactions(1:2:end-1) Reactions(2:2:end)];
idx = data.mesh.coord(:,1) == 0;

X = data.mesh.coord(idx,1);
Y = data.mesh.coord(idx,2);

i=5;
figure
plot(Reactions2plot(:,1),Y,'.b')
hold on
plot([0 0],[0 0.5],'-')
quiver(X(5:i:end-5),Y(5:i:end-5),Reactions2plot(5:i:end-5,1),zeros(size(Reactions2plot(5:i:end-5,1))),"off",'b')
title("Reaction X")
xlabel('Force magnitude')
ylabel('Section position')

figure
plot(Reactions2plot(:,2),Y,'.b')
hold on
plot([0 0],[0 0.5],'-')
quiver(X(5:i:end-5),Y(5:i:end-5),Reactions2plot(5:i:end-5,2),zeros(size(Reactions2plot(5:i:end-5,2))),"off",'b')
title("Reaction Y")
xlabel('Force magnitude')
ylabel('Section position')

RootMoment = sum(Reactions2plot(:,1)'*(Y-0.25));
TestTotalReact = [sum(Reactions2plot(:,1)) sum(Reactions2plot(:,2)) RootMoment]


idxTip = find(data.mesh.coord(:,1) == 1 & data.mesh.coord(:,2) == 0.25);
refU = uTest.fValues(idxTip,2);

Boundary = data.mesh.createBoundaryMesh();
BoundMesh = Boundary{1}.mesh;
% mass matrix
ndim  = data.mesh.ndim;
test  = P1Function.create(BoundMesh,ndim); 
trial = P1Function.create(BoundMesh,ndim); 
s.type    = 'MassMatrix';
s.mesh    = BoundMesh;
s.test    = test;
s.trial   = trial;
lhs       = LHSintegrator.create(s);
M         = lhs.compute();

Stress = M\Reactions;
Stress = [Stress(1:2:end-1) Stress(2:2:end)];
figure
plot(Stress,Y)
title('Stresses at root [x=0]')
xlabel('Stress magnitude')
ylabel('Section position')
legend('Normal stress $\sigma_x$','Shear stress $\tau_{xy}$','interpreter','latex')