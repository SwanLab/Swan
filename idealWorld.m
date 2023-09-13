% FEM
mesh = createMesh();
bc = createBoundaryConditions();
material = createMaterial();
prob = createPhysicalProblem(mesh, material, bc);
solver = createSolver();
solver.solve(problem)
problem.print();
[u, stress, strain] = problem.getVariables();

u.print()
u.plot();

% TopOpt
mesh = createMesh();
bc = createBoundaryConditions(mesh);
desVar = createDesignVariable(mesh);
material = createMaterial(desVar);
prob = createPhysicalProblem(mesh, material, bc);

compliance = createCompliance(u, desVar); % la u ve del prob PERO no l'ha de veure la compliance, passar-la al fer el compute
filter = createFilter(); % a la desVar? s'ha de pensar
volume = createVolume(desVar);
stress = createStress(u);
cost = createCost([compliance, stress]);

cnstr = createConstraint(volume);

optmProb = createOptimizationProblem(prob, cost, cnstr)

optmzr = createOptimizer()
optmzr.solve(optmProb)