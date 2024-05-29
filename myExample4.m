% My example 4
s = struct();
s.testName = 'Bridge_fine_mesh_Nesterov.m';
s.problemCase = 'bridge';
s.x1 = 3;
s.y1 = 1;
s.N   = 150;
s.M   = 50;
s.P    = 1;
s.DoF  = 2;
c = FEMInputWriter(s);
c.createTest();
