clc,clear,close all

%% Load base mesh
file = 'C04';

TMC = TetradecahedronMeshComputer(file);
mesh = TMC.getMesh();
MS = TMC.getMasterSlave();
mesh.plot

%% Create level set

%% set material
E  = 12;
nu = 0.35;
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

% isV2X = @(coor) (abs(coor(:,1) - 0.75) < 1e-5);
% isV2Y = @(coor) (abs(coor(:,2) - 0.5) < 1e-5);
% isV2Z = @(coor) (abs(coor(:,3) - 0) < 1e-5);
% isVertex2 = @(coor) isV2X(coor) & isV2Y(coor) & isV2Z(coor);
% 
% sDir{2}.domain    = @(coor) isVertex2(coor);
% sDir{2}.direction = [2];
% sDir{2}.value     = 0;

% isV3X = @(coor) (abs(coor(:,1) - 0.5) < 1e-5);
% isV3Y = @(coor) (abs(coor(:,2) - 0.25) < 1e-5);
% isV3Z = @(coor) (abs(coor(:,3) - 0) < 1e-5);
% isVertex3 = @(coor) isV3X(coor) & isV3Y(coor) & isV3Z(coor);
% 
% sDir{3}.domain    = @(coor) isVertex3(coor);
% sDir{3}.direction = [3];
% sDir{3}.value     = 0;

% sDir{1}.domain    = @(coor) abs(coor(:,1) - 100) < 1e-5;
% sDir{1}.direction = [1,2,3];
% sDir{1}.value     = 0;

dirichletFun = [];
for i = 1:numel(sDir)
    dir = DirichletCondition(mesh, sDir{i});
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

%% Set solver
% BCAp = BCApplier(s);
% 
% rigModes = RigidBodyFunction.create(mesh,[0.5,0.5,0.5]);
% RFun = rigModes.projectBasisFunctions('P1');
% for i = 1:length(RFun)
%     Rfull(:,i) = reshape(RFun{i}.fValues',[],1);
% end
% 
% for i = 1:size(Rfull,2)
%     R(:,i) = BCAp.fullToReducedVectorDirichlet(Rfull(:,i));
% end
% s.type = 'ELASTIC';
% s.nullSpace = R;
% s.nLevels = 5;
% s.tol = 1e-8;
% s.maxIter = 1;
% p     = pyAMG.create(s);
% 
% sS.preconditioner = p;
% sS.tol = 1e-5;
% solver = PCG(sS);

%% Continue problem definition

s.solverCase = DirectSolver();
s.solverType = 'REDUCED';
s.solverMode = 'FLUC';
fem = ElasticProblemMicro(s);
fem.solve();

%% Postprocess
%Cube
matHomog = fem.Chomog;
density = mesh.computeVolume(); 
% Tetradecaedron
% matHomog = fem.Chomog/0.5;
% density = mesh.computeVolume()/0.5; 

C11 = matHomog(1,1,1,1);
C12 = matHomog(1,1,2,2);
C44 = matHomog(2,3,2,3);
ZenerRatioTetra = 2*C44./(C11-C12)

meanC11 = (matHomog(1,1,1,1)+matHomog(2,2,2,2)+matHomog(3,3,3,3))/3;
meanC12 = (matHomog(1,1,2,2)+matHomog(1,1,3,3)+matHomog(2,2,1,1)+matHomog(2,2,3,3)+matHomog(3,3,1,1)+matHomog(3,3,2,2))/6;
meanC44 = (matHomog(1,2,1,2)+matHomog(1,2,2,1)+matHomog(1,3,1,3)+matHomog(1,3,3,1)+ ...
           matHomog(2,1,2,1)+matHomog(2,1,1,2)+matHomog(2,3,2,3)+matHomog(2,3,3,2)+ ...
           matHomog(3,1,3,1)+matHomog(3,1,1,3)+matHomog(3,2,3,2)+matHomog(3,2,2,3))/12;
meanZenerRatioTetra = 2*meanC44./(meanC11-meanC12);

% Measures of Agyekun
err1 = []; err2 = []; err3 = []; err4 = [];
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                err1 = [err1, abs(matHomog(i,j,k,l)-matHomog(k,l,i,j))];
                if j~=k
                    err2 = [err2, abs(matHomog(i,i,j,k))];
                end
                if i~=j && k~=l
                    err3 = [err3, abs(matHomog(i,i,j,j)-matHomog(k,k,l,l))];
                end
                if i~=j
                    err4 = [err4, abs(2*matHomog(i,j,i,j)+matHomog(i,i,j,j)-matHomog(i,i,i,i))];
                end
            end
        end
    end
end

% Measures of Perle (isotropy conditions)
% err1 = []; err2 = []; err3 = []; err4 = [];
% for i=1:3
%     for j=1:3
%         err1 = [err1, matHomog(i,i,i,i) - matHomog(j,j,j,j)];
%         for k=1:3
%             if i~=j && i~=k && j~=k
%                 err2 = [err2, matHomog(i,i,j,k)];
%             end
%             for l=1:3
%                 if i~=j && k~=l
%                     err3 = [err3, matHomog(i,i,j,j) - matHomog(k,k,l,l)];
%                 end
%                 for p = 1:3
%                     if j~=k && l~=p
%                         err4 = [err4, matHomog(i,i,i,i) - matHomog(j,k,j,k) - matHomog(l,l,p,p)];
%                     end
%                 end
%             end
%         end
%     end
% end

figure()
hold on
plot(err1/matHomog(1,1,1,1))
plot(err2/matHomog(1,1,1,1))
plot(err3/matHomog(1,1,1,1))
plot(err4/matHomog(1,1,1,1))
lgd = legend('abs(Cijkl - Cklij)/C1111','abs(Ciijk)/C1111','abs(Ciijj - Ckkll)/C1111','abs(Cijij+Ciijj-Ciiii)/C1111')
xlabel('Number of combinations of values')
ylabel('Error')
lgd.FontSize = 30;
fontsize(gcf,40,'points')

C = convert2Voigt(matHomog,'Constitutive')
bulk = (C(1,1) + 2*C(1,2))/3
shear = C(4,4)

bulkGeneral = (C(1,1)+C(2,2)+C(3,3)+2*(C(1,2)+C(1,3)+C(2,3)))/9
shearGeneral = (C(1,1)+C(2,2)+C(3,3)-C(1,2)-C(1,3)-C(2,3)+3*(C(4,4)+C(5,5)+C(6,6)))/15