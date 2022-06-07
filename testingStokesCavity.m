%% Stokes: Triangle 2D
% file = 'test2d_stokes_triangle_steady';
% a.fileName = file;
% s = StokesDataContainer(a);
% s.material.nu = 1e3;
% fem = FEM.create(s);
% fem.computeVariables;
% % fem.printPressure(file);
% fem.printVelocity(file);

%% Stokes: Lid-driven cavity
file = 'test2d_stokes_cavity';
a.fileName = file;
s = StokesDataContainer(a);
s.material.nu = 100;
fem = FEM.create(s);
fem.computeVariables;
% fem.printPressure(file);
fem.printVelocity(file);

%% Elastic: Triangle
% fileII = 'test2d_triangle';
% b.fileName = fileII;
% d = FemDataContainer(b);
% femII = FEM.create(d);
% femII.solve();
% femII.print(fileII);