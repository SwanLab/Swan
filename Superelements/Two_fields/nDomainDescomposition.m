clc,clear,close all 

Meshes = [7 15 31 61 121;
                  3 7 15 31  61 ];
TipError            = zeros(1,1,length(Meshes));
StressError         = zeros(2,2,length(Meshes));

%for k=1:length(Meshes)
%% CONVENTIONAL SOLVER
% Create geometry

fileName = 'CantileverAxialLoad';
s.testName    = [fileName,'.m'];
s.problemCase = 'cantilever';
s.x1          = 1;
s.y1          = 0.5;
s.N           = 121;%Meshes(1,k); %x dimension
s.M           = 51;%Meshes(2,k); %y dimension
s.P           = -1;
s.DoF         = 2;

CantileverAxial = FEMInputWriter(s);
CantileverAxial.createTest();

% Create elasticity problem
a.fileName = fileName;
data = FemDataContainer(a);
% fem = FEM.create(data);
% fem.solve();
% fem.print(fileName);

% Get main results
% uTest = fem.uFun(1,1);
% e_test = fem.strainFun(1,1);
% sig_test = fem.stressFun(1,1);

%% Reactions plots fem

% K = fem.LHS;
% u = zeros(length(uTest.fValues)*2,1);
% u(1:2:end-1) = uTest.fValues(:,1);
% u(2:2:end) = uTest.fValues(:,2);
% Reactions = K(fem.boundaryConditions.dirichlet,:)*u;
% Reactions = [Reactions(1:2:end-1) Reactions(2:2:end)];
% idx = data.mesh.coord(:,1) == 0;
% 
% X = data.mesh.coord(idx,1);
% Y = data.mesh.coord(idx,2);

refTestValues = load('TestReferenceData.mat');

%% DOMAIN DECOMPOSITION (N SUBDOMAINS)

in.subdomains = 3;
in.mesh       = data.mesh;
in.mat        = data.material;
in.type       = "Two";

%numberElements(end+1)=data.mesh.nelem;

% Mesh decomposition
[meshDecomposed] = MeshDecomposer(in);
plotMesh(meshDecomposed);

% Stiffness matrices computation
Kinputs.subMesh = meshDecomposed.subMesh;
Kinputs.material = data.material;
[subK] = SubdomainStiffnessMatrix(Kinputs);

% Lagrange matrices computation
Cinputs.subMesh = meshDecomposed.subMesh;
Cinputs.subBoundMesh = meshDecomposed.subBoundMesh;
Cinputs.lambdaType = "P1";
[subC] = SubdomainLagrangeMatrix(Cinputs);

subMatrices.K = subK;
subMatrices.C = subC;

% Global LHS matrix assembly
AKinputs.subK = subMatrices.K;
AKinputs.subdomains = in.subdomains;
[globalK] = AssemblyGlobalStiffness(AKinputs);

ACinputs.type = in.type;
ACinputs.subBoundMesh = meshDecomposed.subBoundMesh;
ACinputs.subC = subMatrices.C;
ACinputs.subK = subMatrices.K;
ACinputs.globalK = globalK;
[globalC] = AssemblyGlobalLagrange(ACinputs);

AGinputs.K = globalK;
AGinputs.C = globalC;
AGinputs.type = in.type;
[globalLHS] = AssemblyGlobalLHS(AGinputs);

% Global RHS matrix assembly
force.magnitude = s.P;
force.dim = 2; % 1 -> Axial force, 2-> Bending force

RHSinputs.subBoundMesh = meshDecomposed.subBoundMesh;
RHSinputs.force = force;
RHSinputs.subK = subMatrices.K;
RHSinputs.globalLHS = globalLHS;
[globalRHS] = AssemblyGlobalRHS(RHSinputs);

% Solver
globalMatrices.LHS = globalLHS;
globalMatrices.RHS = globalRHS;
[globalU] = SolverDomainDecomposition(globalMatrices);

% Decompose solution
decomposeInputs.subMatrices = subMatrices;
decomposeInputs.type = in.type;
decomposeInputs.super = false;
decomposeInputs.globalU = globalU;
[Fields] = DecomposeSolution(decomposeInputs);

% Compute stresses
stressInputs.Fields = Fields;
stressInputs.subMatrices = subMatrices;
stressInputs.meshDecomposed = meshDecomposed;
stressInputs.plot = true;
stressInputs.type = in.type;
[totalReact, stress] = ComputeReactionsAndStresses(stressInputs);

% Compute errors
ErrorConnInputs.globalU = globalU;
ErrorConnInputs.globalK = globalK;
ErrorConnInputs.globalLHS = globalLHS;
ErrorConnInputs.subC = subC;
[ErrorConnL2, ErrorConnLinf] = ErrorConnection(ErrorConnInputs);

Results.stress = stress;
Results.u = Fields.u;
super = false;
[DispError, StressError] = ErrorComputationL2(refTestValues,Results,meshDecomposed);
plotStresses(meshDecomposed.subBoundMesh,stress,in.type);

% Print results
PrintInputs.subMesh = meshDecomposed.subMesh ;
PrintInputs.U = Fields.u;
PrintInputs.type = in.type;
PrintResults(PrintInputs);








