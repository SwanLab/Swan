clc,clear,close all

matHomog = zeros(3,3,3,3,9);
density = zeros(1,9);
bulk = zeros(1,9);
shear = zeros(1,9);
err.e1 = zeros(1,9);
err.e2 = zeros(1,9);
err.e3 = zeros(1,9);
err.e4 = zeros(1,9);
err.e5 = zeros(1,9);

for i=1:9
%% Define mesh
file = ['C0',num2str(i)];
fprintf(['Start',file,'\n'])
TMC = TetradecahedronMeshComputer(file);
mesh = TMC.getMesh();
MS = TMC.getMasterSlave();
figure(i)
mesh.plot
drawnow

%% Create material
E  = 20.8;
nu = 0.3;
young   = ConstantFunction.create(E,mesh);
poisson = ConstantFunction.create(nu,mesh);
s.type    = 'ISOTROPIC';
s.ptype   = 'ELASTIC';
s.ndim    = mesh.ndim;
s.young   = young;
s.poisson = poisson;
tensor    = Material.create(s);
material  = tensor;

%% set boundary conditions
isV1X = @(coor) (abs(coor(:,1) - 0.15) < 1e-5);
isV1Y = @(coor) (abs(coor(:,2) - 0.3750)  < 1e-5);
isV1Z = @(coor) (abs(coor(:,3) - 0.9625)  < 1e-5);
isVertex1 = @(coor) isV1X(coor) & isV1Y(coor) & isV1Z(coor);

sDir{1}.domain    = @(coor) isVertex1(coor);
sDir{1}.direction = [1,2,3];
sDir{1}.value     = 0;

dirichletFun = [];
for j = 1:numel(sDir)
    dir = DirichletCondition(mesh, sDir{j});
    dirichletFun = [dirichletFun, dir];
end
s.dirichletFun = dirichletFun;
s.pointloadFun = [];
s.periodicFun  = 1; %Set to not be empty
s.mesh = mesh;
bc = BoundaryConditions(s);
bc.updatePeriodicConditions(MS);

%% Set micro problem
s.mesh     = mesh;
s.material = material;
s.scale    = 'MICRO';
s.dim      = '3D';
s.boundaryConditions = bc;

%% Continue problem definition
fprintf(['Solving step ',num2str(i),'\n'])
s.solverCase = DirectSolver();
s.solverType = 'REDUCED';
s.solverMode = 'FLUC';
fem = ElasticProblemMicro(s);
fem.solve();

%% Postprocess
% Save material properties
matHomog(:,:,:,:,i) = fem.Chomog;
density(i) = mesh.computeVolume(); 

C = convert2Voigt(matHomog(:,:,:,:,i),'Constitutive');
bulk(i) = (C(1,1)+C(2,2)+C(3,3)+2*(C(1,2)+C(1,3)+C(2,3)))/9;
shear(i) = (C(1,1)+C(2,2)+C(3,3)-C(1,2)-C(1,3)-C(2,3)+3*(C(4,4)+C(5,5)+C(6,6)))/15;

% Error computation
C1111 = matHomog(1,1,1,1,i);
C2222 = matHomog(2,2,2,2,i);
C2323 = matHomog(2,3,2,3,i);
C1313 = matHomog(1,3,1,3,i);
C1122 = matHomog(1,1,2,2,i);
C1123 = matHomog(1,1,2,3,i);
C1113 = matHomog(1,1,1,3,i);

ZenerRatio(i) = 2*C2323./(C1111-C1122);

% Measures of Agyekun
err.e1(i) = C1123/C1111;
err.e2(i) = C1113/C1111;
err.e3(i) = abs(C1111-C2222)/C1111;
err.e4(i) = abs(C2323-C1313)/C1111;
err.e5(i) = abs(2*C2323+C1122-C1111)/C1111;
end
save('TetradecahedronHomogenized','matHomog','density','err','bulk','shear')

%% Plots
figure(10)
hold on
plot(density,err.e1)
plot(density,err.e2)
lgd = legend('C1123/C1111','C1113/C1111')
xlabel('Density')
ylabel('Error')
lgd.FontSize = 30;
fontsize(gcf,40,'points')

figure(11)
hold on
plot(density,err.e3)
plot(density,err.e4)
plot(density,err.e5)
lgd = legend('abs(C1111-C2222)/C1111','abs(C2323-C1313)/C1111','abs(C2323+C1122-C1111)/C1111')
xlabel('Density')
ylabel('Error')
lgd.FontSize = 30;
fontsize(gcf,40,'points')


k3D  = @(nu) E/(3-6*nu);
mu3D = @(nu) E/(2*(1+nu));
l3D  = @(nu) E*nu/((1+nu)*(1-2*nu));
kUB3D   = @(phi,nu) k3D(nu)*(1-phi).*((2.*mu3D(nu)+l3D(nu)-k3D(nu))./(2.*mu3D(nu)+l3D(nu)-k3D(nu).*(1-phi)));
muUB3D  = @(phi,nu) mu3D(nu)*(1-phi).*((10.*(2.*mu3D(nu)+l3D(nu)) - 4.*(k3D(nu)+2.*mu3D(nu)))./(10.*(2.*mu3D(nu)+l3D(nu))-(1-phi).*4.*(k3D(nu)+2.*mu3D(nu))));

figure(12)
hold on
plot((1-density),bulk)
fplot(@(phi) kUB3D(phi,0.3),[0 1])
lgd = legend('Homog','H-S')
xlabel('Damage')
ylabel('Bulk')
lgd.FontSize = 30;
fontsize(gcf,40,'points')

figure(13)
hold on
plot((1-density),shear)
fplot(@(phi) muUB3D(phi,0.3),[0 1])
lgd = legend('Homog','H-S');
xlabel('Damage')
ylabel('Shear')
lgd.FontSize = 30;
fontsize(gcf,40,'points')

figure(14)
plot(density,ZenerRatio)
xlabel('Density')
ylabel('Zener Ratio')
fontsize(gcf,40,'points')