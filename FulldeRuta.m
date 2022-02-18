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

%% Todo
% change name cartd per dNdx
% LHSintegrator_Stiffness : create geometry in init.
%try sparse vs accumarray in assemblyCmat, assembleMatrix (BMatrixComputer) and LHSintegrator 
% posI, posJ in LHSintergrator_StiffnessElasticStoredB
%LHSintergrator_StiffnessElastic (computeB) 
% 6. DOF+ Boudnary
% 4. Create input data for NewFem (mesh,bc,material)
% investigate why nnode is not correct in TopOptTests (reads 3, should be 2)
% make test3d_hexahedra go faster
% % time lost during squeeze (+33.000 calls) in computeElementalLHS

%%

% Inputa data
% DOF + bondary
% Stifness with integrator
% Forces * BC inside + integrator

%Top uses NewFem instaead of FEM

% 1. Vectorized test fix {{done}}
%   {{Done}}, committed and pushed. Pending approval of pull request.

% 2. Delete KgeneratorWithFullStoredB {{done}}
%   {{Done}}, committed and pushed. Pending approval of pull request.

% 3. Clean LHSintegratorStiffnessElasticStoredB {{done}} (?)
%   {{Done}} (?), committed and pushed. Pending approval of pull request.

% 4. Create input data for NewFem (mesh,bc,material)
%    {{Done}}-ish. It needs:
%         - mesh
%         - pdim
%         - ptype, scale
%         - problemID (fileName) for DOF compatibility
%         - dirichlet
%         - pointload

%% 5. Create  LHSintergrator_StiffnessElastic 
% 6. DOF+ Boudnary
% 7. Force use integrator for assembly
% 8. LHS integrator must be composed by integrator_simple for assemnbly