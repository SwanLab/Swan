clc,clear,close all

%% CONVENTIONAL SOLVER
% Create geometry
fileName = 'CantileverAxialLoad';
s.testName    = [fileName,'.m'];
s.problemCase = 'cantilever';
s.x1          = 1;
s.y1          = 0.5;
s.N           = 121;
s.M           = 51;
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
% 
% % Get main results
% uTest = fem.uFun(1,1);
% e_test = fem.strainFun(1,1);
% sig_test = fem.stressFun(1,1);

refTestValues = load('TestReferenceData.mat');

%% DOMAIN DECOMPOSITION (N SUBDOMAINS)
in.subdomains = 3;
in.mesh       = data.mesh;
in.mat        = data.material;
in.type       = "Three";

% Mesh decomposition
[meshDecomposed] = MeshDecomposer(in);
plotMesh(meshDecomposed);

% Stiffness matrices computation
Kinputs.subMesh = meshDecomposed.subMesh;
Kinputs.material = data.material;
[subK] = SubdomainStiffnessMatrix(Kinputs);

% Mass matrices computation
Cinputs.subMesh = meshDecomposed.subMesh;
Cinputs.subBoundMesh = meshDecomposed.subBoundMesh;
Cinputs.lambdaType = "P1";
[subC] = SubdomainLagrangeMatrix(Cinputs);

% Interface matrices computation
Dinputs.subBoundMesh = meshDecomposed.subBoundMesh;
Dinputs.interMesh = meshDecomposed.interMesh;
Dinputs.lambdaType = Cinputs.lambdaType;
[subD] = SubdomainInterfaceMatrix(Dinputs);

subMatrices.K = subK;
subMatrices.C = subC;
subMatrices.D = subD;

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

ADinputs.subD = subMatrices.D;
[globalD] = AssemblyGlobalInterface(ADinputs);

LHS.K = globalK;
LHS.C = globalC;
LHS.D = globalD;
LHS.type = in.type;
[globalLHS] = AssemblyGlobalLHS(LHS);

% Global RHS matrix assembly
force.magnitude = s.P; 
force.dim = 2; % 1 -> Axial force, 2-> Bending force

RHS.subBoundMesh = meshDecomposed.subBoundMesh;
RHS.force = force;
RHS.subK = subMatrices.K;
RHS.globalLHS = globalLHS;
[globalRHS] = AssemblyGlobalRHS(RHS);

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
PrintInputs.interMesh = meshDecomposed.interMesh;
PrintInputs.U = Fields.u;
PrintInputs.interU = Fields.uI;
PrintInputs.type = in.type;
PrintResults(PrintInputs);

