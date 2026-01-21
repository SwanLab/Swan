clear;
close all;

file = 'punzon';
a.fileName = file;
s = FemDataContainer(a);
% mesh=createmesh(file);
% [young, poisson]=computeElasticProperties(mesh);
% createElasticProblem(file,mesh);
% bc=createBoundaryConditions(file,mesh);
% 
% 
% d = DirichletCondition(s.mesh,[],s.bc.dirichlet);
% s.boundaryConditions.dirichlet_dofs = d.dofs;
% s.boundaryConditions.dirichlet_vals = d.values;
% 
% t = TractionLoad(s.mesh,s.bc,'DIRAC');
% s.boundaryConditions.tractionFun = t;

fem = PhysicalProblem.create(s);
fem.solve();



%% Results
fem.uFun.plot()

fem.uFun.print('results_fem_dispPunzon', 'GiD') % print using GiD
fem.uFun.print('results_fem_dispPunzon', 'Paraview') % print using Paraview

fem.print('results_femPunzon', 'Paraview') % print using Paraview


%% Functions

function mesh= createMesh(file)
    a.fileName = file;
    s = FemDataContainer(a);
    m=s.mesh;
    con = m.connec;
    q = Quadrature.create(m,0);
    dv = m.computeDvolume(q);
    negElem=find(dv<=0);
    con(negElem,:) = [];
    sM.coord = m.coord;
    sM.connec = con;
    m2 = Mesh.create(sM);
    m2 = m2.computeCanonicalMesh();
    obj.mesh = m2;
end

function [young,poisson]=computeElasticProperties(mesh)
    E  = 1; %1; % canviar?
    nu = 1/3; % canviar?
    young   = ConstantFunction.create(E,mesh);
    poisson = ConstantFunction.create(nu,mesh);
end


function createElasticProblem(obj,file,mesh)
    s.mesh = mesh;
    s.scale = 'MACRO';
    s.material = obj.material;
    s.dim = '3D';
    s.boundaryConditions = createBoundaryConditions(file,mesh);
    s.interpolationType = 'LINEAR';
    s.solverType = 'REDUCED';
    s.solverMode = 'DISP';
    s.solverCase = 'DIRECT';
    fem = ElasticProblem(s);
    xD = obj.designVariable.obtainDomainFunction();
    obj.material.setDesignVariable(xD);
    C   = obj.material.obtainTensor();
    fem.updateMaterial(C);
    obj.physicalProblem = fem;
end

function bc = createBoundaryConditions(fileName,mesh)
    femReader = FemInputReaderGiD();
    s         = femReader.read(fileName);
    sPL       = obj.computeCondition(s.pointload);
    sDir      = obj.computeCondition(s.dirichlet);
    dirichletFun = [];
    for i = 1:numel(sDir)
        dir = DirichletCondition(mesh, sDir{i});
        dirichletFun = [dirichletFun, dir];
    end
    s.dirichletFun = dirichletFun;
    pointloadFun = [];
    for i = 1:numel(sPL)
        pl = PointLoad(mesh, sPL{i});
        pointloadFun = [pointloadFun, pl];
    end
    s.pointloadFun = pointloadFun;
    s.periodicFun  = [];
    s.mesh         = mesh;
    bc = BoundaryConditions(s);
end