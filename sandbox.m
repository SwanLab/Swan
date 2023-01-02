file = 'test_micro_holeinclusion';
% file = 'test2d_triangle';
a.fileName = file;
s = FemDataContainer(a);
fem = FEM.create(s);
fem.computeChomog();
% fem.solve();
fem.print(file)