function ComputingTopOpt


% fileName = 'bridge_10_2';
fileName = 'CantileverArnau2';
% Data input
s.testName = [fileName,'.m'];
s.x1       = 2;
s.y1       = 1;
s.N        = 60;
s.M        = 30;
s.P        = -100;
s.DoF      = 2;
s.problemCase = 'cantilever';
% s.testName = [fileName,'.m'];
% s.x1       = 10;
% s.y1       = 2;
% s.N        = 100;
% s.M        = 20;
% s.P        = -100;
% s.DoF      = 2;
% s.problemCase = 'bridge';

FEMWriter = FEMInputWriter(s);
FEMWriter.createTest;
s.testName = 'test_cantilever2';%''testJose';
t = TopOptComputer(s);
t.compute();
end
