%% New FEM workflow
clc; clear; close all
% (1) Create the mesh
    mesh = UnitTriangleMesh(5,5);

% (2) Set boundary conditions
% bc = createBoundaryConditions();
    dirichletNodes = abs(mesh.coord(:,1)-0) < 1e-12;
    rightSide  = max(mesh.coord(:,1));
    isInRight = abs(mesh.coord(:,1)-rightSide)< 1e-12;
    isInMiddleEdge = abs(mesh.coord(:,2)-0.5) < 0.1;
    forceNodes = isInRight & isInMiddleEdge;

    nodes = 1:mesh.nnodes;
    bcDir = [nodes(dirichletNodes)';nodes(dirichletNodes)'];
    nodesdir=size(nodes(dirichletNodes),2);
    bcDir(1:nodesdir,end+1) = 1;
    bcDir(nodesdir+1:end,end) = 2;
    bcDir(:,end+1)=0;
    bc.dirichlet = bcDir;
    bc.pointload(:,3) = -1;
    bc.pointload(:,2) = 2;
    bc.pointload(:,1) = nodes(forceNodes);


% (3) Create the material
% material = createMaterial();
    I = ones(mesh.nelem,1);
    s.ptype = 'ELASTIC';
    s.pdim  = '2D';
    s.nelem = mesh.nelem;
    s.mesh  = mesh;
    s.kappa = .9107*I;
    s.mu    = .3446*I;
    mat = Material.create(s);
    mat.compute(s);
    material = mat;

% (4) Create the problem
% prob = createPhysicalProblem(mesh, material, bc);
    s.mesh = mesh;
    s.type = 'ELASTIC';
    s.scale = 'MACRO';
    s.material = material;
    s.bc = bc;
%     problem = FEM.create(s);
    problem = NewElasticProblem(s);

% (5) Create the solver and solve
    ProblemSolver.solve(problem)

% Post-processing
% problem.print();
[u, stress, strain] = problem.getVariables();


%% random thoughts and stuff
% - time for boundary conditions
%       - pass coords of point load -> coord2dof matching
%       - logical functions for boundaries?

% - materials
%       - kappa + mu // E + nu
%       - C -> nstre x nstre x nelem x ngaus rly?
%       - how would we do multimaterial? (2 years ago...)
%           - generate 2 materials, each with a nstre x nstre matrix
%           - p0 function
%       - homogenizedvar...