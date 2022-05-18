function ComputingTopOpt

rho0Name = 'MicroUltraFine.mat';
jumpTo2ndPart = true;

if jumpTo2ndPart == false

%     fileName = 'jaCantilever';
%     % Data input
%     s.testName = [fileName,'.m'];
%     s.x1       = 2;
%     s.y1       = 1;
%     s.N        = 150;
%     s.M        = 75;
%     s.P        = -100;
%     s.DoF      = 2;
% 
%     FEMWriter = FEMInputWriter(s);
%     FEMWriter.createTest;




    fileName = 'test_anisotropy_cantilever';
    s = Settings(fileName);
    s.warningHoleBC = false;
    s.printIncrementalIter = false;
    s.printChangingFilter = false;
    s.printing = false;
    translator = SettingsTranslator();
    translator.translate(s);
    fileName = translator.fileName;
    settings  = SettingsTopOptProblem(fileName);
    topOptSolver = TopOpt_Problem(settings);
    while topOptSolver.incrementalScheme.hasNext()
        topOptSolver.incrementalScheme.next();
        topOptSolver.optimizer.solveProblem();
    end

    rho0 = topOptSolver.designVariable.value;
    save(rho0Name,'rho0');
else
    load(rho0Name);
    fileName = 'test_anisotropy_cantilever_rho0';

    s = Settings(fileName);
    s.warningHoleBC = false;
    s.printIncrementalIter = false;
    s.printChangingFilter = false;
    s.printing = false;
    translator = SettingsTranslator();
    translator.translate(s);
    fileName = translator.fileName;
    settings  = SettingsTopOptProblem(fileName);
    DesignVariable = convertCharsToStrings(settings.designVarSettings.type);
    if  DesignVariable == "Density"
        settings.designVarSettings.creatorSettings.type = 'Given';
        settings.designVarSettings.creatorSettings.rho0 = rho0;
    elseif DesignVariable == "LevelSet"
        settings.designVarSettings.initialCase = 'given';
        settings.designVarSettings.creatorSettings.value = rho0;
    end
    topOptSolver = TopOpt_Problem(settings);
    while topOptSolver.incrementalScheme.hasNext()
        topOptSolver.incrementalScheme.next();
        topOptSolver.optimizer.solveProblem();
    end
end

end