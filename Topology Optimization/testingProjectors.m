%% Create sample FEM results
file = 'test2d_triangle';
a.fileName = file;
s = FemDataContainer(a);
fem = FEM.create(s);
fem.solve();

strain = fem.variables.strain;
u = fem.variables.d_u;

%% P1 to P0
aa.connec = s.mesh.connec;
aa.type   = s.mesh.type;
aa.fNodes = u;
fefunDisp = FeFunction(aa);
p1displac = fefunDisp.computeValueInCenterElement();

%% P0 to P1