%% Simple topology optimization example
% Note that the beta term

function runningSimpleTopOpt
tic
iter = 1;
% neededIter = zeros(102,1);
xNew = zeros(20200,1);
% solutions = zeros(20200,102);

% save('CounterPol\solutions.mat', 'solutions');
save('CounterPol\iter.mat', 'iter');
% save('CounterPol\neededIter.mat', 'neededIter');
save('CounterPol\xNew.mat', 'xNew');
clear
s.maxIter = 100;
s.TOL = 1e-12;
s.topOptProblem = createFullTopOptProblem();
solver = SimpleShapeOptimizationSolver(s);
solver.solve();
toc
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

