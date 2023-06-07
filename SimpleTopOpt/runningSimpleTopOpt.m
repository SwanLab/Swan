%% Simple topology optimization example

function runningSimpleTopOpt
clear
s.maxIter = 500;
s.TOL = 1e-12;
s.topOptProblem = createFullTopOptProblem();
solver = SimpleShapeOptimizationSolver(s);
solver.solve();
end


function t = createFullTopOptProblem()
% settings = Settings('Example1');
settings = Settings('Example1');
translator = SettingsTranslator();
translator.translate(settings);
fileName = translator.fileName;
settingsTopOpt = SettingsTopOptProblem(fileName);
t = TopOpt_Problem(settingsTopOpt);
end

