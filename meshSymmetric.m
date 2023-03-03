s.testName    = 'symmetricVerticalCantilever.m';
s.problemCase = 'cantilever';
s.x1          = 0.5;
s.y1          = 1;
s.N           = 80;
s.M           = 160;
s.P           = 15;
s.DoF         = 2;

Problem = FEMInputWriter(s);