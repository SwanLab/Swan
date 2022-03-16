function ComputingTopOpt

fileName = 'CantileverArnau';
% Data input
s.testName = [fileName,'.m'];
s.x1       = 2;
s.y1       = 1;
s.N        = 30;
s.M        = 30;
s.P        = -100;
s.DoF      = 2;

FEMWriter = FEMInputWriter(s);
FEMWriter.createTest;




s.testName = 'test_cantilever2';%'testJose';%'
t = TopOptComputer(s);
t.compute();
end