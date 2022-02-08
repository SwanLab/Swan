%Mesh

%BC(dirichlet, newuman)

% Material

s.kappa = ...
s.mu = ...
s.mesh = 
%s.E   = 
%s.n_u = ...
Material(s)

% FEM
s.mesh = mesh;
s.boundaryConditions = BC;
s.material = material;
s.type = 'ELASTIC';
s.dim = '2D';
fem = FEM.create(s);
fem().solve();
fem.plot();

% Inputa data
% DOF + bondary
% Stifness with integrator
% Forces * BC inside + integrator

%Top uses NewFem instaead of FEM

% In NewElastic replace integrator by LHSintegrator
% Create LHSintergrator_StiffnessElasticOld child of LHSintergrator_Stiffness
 %with  its own compute wiht info from computeKgenerator
% Create  LHSintergrator_StiffnessElastic 
% LHS integrator is composed by integrator_simple for assemnbly