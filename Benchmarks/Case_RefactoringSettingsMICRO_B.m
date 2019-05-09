cParams.problemData.problemFileName = 'test2d_micro';
cParams.problemData.scale = 'MICRO';

cParams.designVarSettings.type        = 'LevelSet';
cParams.designVarSettings.initialCase = 'circleInclusion';

cParams.incrementalSchemeSettings.nSteps = 1;

cParams.costSettings.shapeFuncSettings = {
    struct('type','chomog_alphabeta','alpha',[1 0 0]','beta',[0 -1 0]')
    struct('type','perimeterConstraint','PerimeterTarget',5)
    };
cParams.costSettings.weights = [1 0.1];

cParams.constraintSettings.shapeFuncSettings = {struct('type','volumeConstraint')};

cParams.optimizerSettings.type       = 'AlternatingPrimalDual';
cParams.optimizerSettings.shallPrint = false;
cParams.optimizerSettings.settingsMonitor.showOptParams       = true;
cParams.optimizerSettings.settingsMonitor.refreshInterval     = 2;
cParams.optimizerSettings.settingsMonitor.shallDisplayDesignVar       = true;
cParams.optimizerSettings.settingsMonitor.shallShowBoundaryConditions = true;

cParams.optimizerSettings.uncOptimizerSettings.type       = 'SLERP';