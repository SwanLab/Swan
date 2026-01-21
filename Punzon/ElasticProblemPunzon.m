clear;
close all;


file = 'punzon';
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


d = DirichletCondition(s.mesh,[],s.bc.dirichlet);
s.boundaryConditions.dirichlet_dofs = d.dofs;
s.boundaryConditions.dirichlet_vals = d.values;

t = TractionLoad(s.mesh,s.bc,'DIRAC');
s.boundaryConditions.tractionFun = t;

fem = PhysicalProblem.create(s);
fem.solve();


%% Results
fem.uFun.plot()

fem.uFun.print('results_fem_dispPunzon', 'GiD') % print using GiD
fem.uFun.print('results_fem_dispPunzon', 'Paraview') % print using Paraview

fem.print('results_femPunzon', 'Paraview') % print using Paraview