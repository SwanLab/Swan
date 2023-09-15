%% New FEM workflow
% (1) Create the mesh
mesh = UnitTriangleMesh(5,5);
mesh.plot()

% (2) Set boundary conditions
bc = createBoundaryConditions();

% (3) Create the material
material = createMaterial();

% (4) Create the problem
prob = createPhysicalProblem(mesh, material, bc);

% (5) Create the solver and solve
solver = createSolver();
solver.solve(problem)

% Post-processing
problem.print();
[u, stress, strain] = problem.getVariables();