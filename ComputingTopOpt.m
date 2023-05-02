function ComputingTopOpt


% fileName = 'bridge_10_2';
fileName = 'CantileverArnau2';
% Data input
s.testName = [fileName,'.m'];
s.x1       = 2;
s.y1       = 1;
s.N        = 80;
s.M        = 40;
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
% s.testName = [fileName,'.m'];
% s.x1       = 2;
% s.y1       = 1;
% s.N        = 80;
% s.M        = 40;
% s.P        = -100;
% s.DoF      = 2;
% s.problemCase = 'arch';

% fileName = 'CantileverArnau3';
% Data input
% s.testName = [fileName,'.m'];
% s.x1       = 2;
% s.y1       = 1;
% s.z1       = 1;
% s.N        = 70;
% s.M        = 35;
% s.
% s.P        = -100;
% s.DoF      = 2;
% s.problemCase = 'cantilever3';

% FEMWriter = FEMInputWriter(s);
% FEMWriter.createTest;
s.testName = 'test_cantilever_IPM';%''testJose';
% s.testName = 'test_bridg_multipleLoads';
t = TopOptComputer(s);
t.compute();
end