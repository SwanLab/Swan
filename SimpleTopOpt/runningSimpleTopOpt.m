%% Simple topology optimization example

function runningSimpleTopOpt
clear
s.maxIter = 100;
s.TOL = 1e-12;
p.testName = 'LoopOfTests.m';
p.DoF = 2; % '2' for 2D or '3' for 3D
p = uploadMeshCases(p);
timeData = zeros(size(p.Nvec,1),1);
for i = 1:length(p.Nvec)
    tic
    p.iCase = i;
    p.N = p.Nvec(i);
    p.M = p.Mvec(i);
    p.O = p.Ovec(i);
    Problem = FEMInputWriter(p);
    Problem.createTest();
    s.topOptProblem = createFullTopOptProblem();
    solver = SimpleShapeOptimizationSolver(s);
    solver.solve();
    dim = char(string(p.DoF));
    divx = char(string(p.N));
    divy = char(string(p.M));
    divz = char(string(p.O));
    timeData(i) = toc
    saveName = [p.problemCase,'_',dim,'D_',divx,'x',divy,'x',divz];    
    % save(['SimpleTopOpt/Results/',saveName,'.mat'],'solver');    
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
p.z1 = 0;
p.P  = -1; % load
p.Nvec = [200];
p.Mvec = [200];
p.Ovec = 0;
% p.Nvec = [100, 150, 200, 250, 300, 350, 400]; % divisions along x direction
% p.Mvec = [100, 150, 200, 250, 300, 350, 400]; % divisions along y direction
% p.Ovec = [0, 0, 0, 0, 0, 0, 0, 0, 0]; % divisions along z direction (if needed)
end

