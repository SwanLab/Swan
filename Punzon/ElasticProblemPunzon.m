clear;
close all;


file = 'punzon';
a.fileName = file;
s = FemDataContainer(a);
s.boundaryConditions.dirichlet_dofs

fem = PhysicalProblem.create(s);
fem.solve();


%% Results
fem.uFun.plot()

fem.uFun.print('results_fem_dispPunzon', 'GiD') % print using GiD
fem.uFun.print('results_fem_dispPunzon', 'Paraview') % print using Paraview

fem.print('results_femPunzon', 'Paraview') % print using Paraview