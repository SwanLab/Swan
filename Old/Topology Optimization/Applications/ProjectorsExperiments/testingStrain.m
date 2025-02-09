fileII = 'test2d_triangle';
b.fileName = fileII;
d = FemDataContainer(b);
femII = FEM.create(d);
femII.solve();
femII.print(fileII);