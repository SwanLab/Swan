a.fileName = 'test2d_simpleTriangleLinear';
s = FemDataContainer(a);
fem = FEM.create(s);
fem.solve();