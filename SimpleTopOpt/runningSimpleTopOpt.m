%% Simple topology optimization example
% Note that the beta term
function runningSimpleTopOpt
s.maxIter = 100;
s.TOL = 1e-12;
s.topOptProblem = createFullTopOptProblem();
solver = SimpleShapeOptimizationSolver(s);
solver.solve();

end


function t = createFullTopOptProblem()
settings = Settings('Example1');
translator = SettingsTranslator();
translator.translate(settings);
fileName = translator.fileName;
settingsTopOpt = SettingsTopOptProblem(fileName);
t = TopOpt_Problem(settingsTopOpt);
end