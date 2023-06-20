clc,clear,close all 

%% CONVENTIONAL SOLVER
% Create geometry
fileName = 'CantileverBendingLoad';
s.testName    = [fileName,'.m'];
s.problemCase = 'cantilever';
s.x1          = 1;
s.y1          = 0.5;
s.N           = 121;
s.M           = 51;
s.P           = -1;
s.DoF         = 2;

CantileverBending = FEMInputWriter(s);
CantileverBending.createTest();

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
in.mesh = data.mesh;
in.type = 'Two';

% Mesh decomposition
[meshDecomposed] = MeshDecomposer(in);

% Stiffness submatrices computation
Kinputs.subMesh = meshDecomposed.subMesh;
Kinputs.material = data.material;
[subK] = SubdomainStiffnessMatrix(Kinputs);

% Static condensation
CondensInputs.subK = subK;
CondensInputs.subBoundMesh = meshDecomposed.subBoundMesh;
[SuperElemInfo] = SuperelementStiffnessMatrix(CondensInputs);

%Lagrange submatrices computation
Cinputs.subMesh = meshDecomposed.subMesh;
Cinputs.subBoundMesh = meshDecomposed.subBoundMesh;
Cinputs.lambdaType = "P1";
[subC] = SubdomainLagrangeMatrix(Cinputs);

subMatrices.K = SuperElemInfo.superK;
subMatrices.C = subC;

% Global LHS matrix assembly
AKinputs.subK = subMatrices.K;
AKinputs.subdomains = in.subdomains;
[globalK] = AssemblyGlobalStiffness(AKinputs);

ACinputs.subC = subMatrices.C;
[globalC] = AssemblyGlobalLagrangeCondensed(ACinputs);

AGinputs.K = globalK;
AGinputs.C = globalC;
AGinputs.type = in.type;
[globalLHS] = AssemblyGlobalLHS(AGinputs);

% Global RHS matrix assembly
force.magnitude = s.P;
force.dim = 2;% 1 -> Axial force, 2-> Bending force

RHSinputs.subK = SuperElemInfo.superK;
RHSinputs.force = force;
RHSinputs.dofsBoundary = SuperElemInfo.dofs.boundary;
RHSinputs.globalLHS = globalLHS;
[globalRHS] = AssemblyGlobalSuperelementsRHS(RHSinputs);

% Solver
GlobalMatrices.LHS = globalLHS;
GlobalMatrices.RHS = globalRHS;
[globalU] = SolverDomainDecomposition(GlobalMatrices);

% Decompose solution
decomposeInputs.dofsBoundary = SuperElemInfo.dofs.boundary;
decomposeInputs.subMatrices = subMatrices;
decomposeInputs.type = in.type;
decomposeInputs.super = true;
decomposeInputs.globalU = globalU;
[Fields] = DecomposeSolution(decomposeInputs);

% Reconstruction
Reconstruct.globalU = globalU;
Reconstruct.superElemInfo = SuperElemInfo;
[uFull] = FullMeshReconstruction(Reconstruct);


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
Results.u = uFull;
super = true;
[DispError, StressError] = Error. ComputationL2(refTestValues,Results,meshDecomposed);
plotStresses(meshDecomposed.subBoundMesh,stress,in.type);

% Print results
PrintInputs.subMesh = meshDecomposed.subMesh ;
PrintInputs.U = uFull;
PrintInputs.type = in.type;
PrintResults(PrintInputs);