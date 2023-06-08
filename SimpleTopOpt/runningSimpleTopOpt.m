%% Simple topology optimization example

function runningSimpleTopOpt
clear
s.maxIter = 5;
s.TOL = 1e-12;
p.testName = 'LoopOfTests.m';
p.DoF = 2; % '2' for 2D or '3' for 3D
p = uploadMeshCases(p);
for i = 1:length(p.Nvec)
    p.iCase = i;
    p.N = p.Nvec(i);
    p.M = p.Mvec(i);
    Problem = FEMInputWriter(p);
    Problem.createTest();
    s.topOptProblem = createFullTopOptProblem();
    solver = SimpleShapeOptimizationSolver(s);
    solver.solve();
    dim = char(string(p.DoF));
    divx = char(string(p.N));
    divy = char(string(p.M));
    saveName = [p.problemCase,'_',dim,'D_',divx,'x',divy];
    save(['SimpleTopOpt/Results/',saveName,'.mat'],'solver');
    close  all
end
end


function t = createFullTopOptProblem()
settings = Settings('Example1');
translator = SettingsTranslator();
translator.translate(settings);
fileName = translator.fileName;
settingsTopOpt = SettingsTopOptProblem(fileName);
t = TopOpt_Problem(settingsTopOpt);
end

function p = uploadMeshCases(p)
p.problemCase = 'cantilever';
p.x1 = 2;
p.y1 = 1;
p.P  = -10; % load
p.Nvec = [20, 40, 80, 200]; % divisions along x direction
p.Mvec = [10, 20, 40, 100]; % divisions along y direction
end

