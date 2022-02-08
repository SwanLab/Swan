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

% 1. Vectorized test fix
% 2. Delete KgeneratorWithFullStoredB
% 3. Clean LHSintegratorStiffnessElasticStoredB
% 4. Create input data for NewFem (mesh,bc,material)
% 5. Create  LHSintergrator_StiffnessElastic 
% 6. DOF+ Boudnary
% 7. Force use integrator for assembly
% 8. LHS integrator must be composed by integrator_simple for assemnbly