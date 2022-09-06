function ComputingTopOpt

% fileName = 'Cantilever';
% % Data input
% s.testName = [fileName,'.m'];
% s.x1       = 2;
% s.y1       = 1;
% s.N        = 80;
% s.M        = 40;
% s.P        = -100;
% s.DoF      = 2;
% 
% FEMWriter = FEMInputWriter(s);
% FEMWriter.createTest;


s.testName = 'IsotropyTest2';%''testJose';
t = TopOptComputer(s);
t.compute();
Ch = t.computation.cost.shapeFunctions{1,1}.Chomog;
end