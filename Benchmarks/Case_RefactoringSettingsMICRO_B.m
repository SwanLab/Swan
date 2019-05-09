problemData.problemFileName = 'test2d_micro';
problemData.scale = 'MICRO';

designVarSettings.type        = 'LevelSet';
designVarSettings.initialCase = 'circleInclusion';

incrementalSchemeSettings.nSteps = 1;

costSettings.shapeFuncSettings = {
    struct('type','chomog_alphabeta','alpha',[1 0 0]','beta',[0 -1 0]')
    struct('type','perimeterConstraint','PerimeterTarget',5)
    };
costSettings.weights = [1 0.1];

constraintSettings.shapeFuncSettings = struct('type','volumeConstraint');

optimizerSettings.type       = 'AlternatingPrimalDual';
optimizerSettings.shallPrint = false;
optimizerSettings.settingsMonitor.showOptParams       = true;
optimizerSettings.settingsMonitor.refreshInterval     = 2;
optimizerSettings.settingsMonitor.shallDisplayDesignVar       = true;
optimizerSettings.settingsMonitor.shallShowBoundaryConditions = true;

optimizerSettings.uncOptimizerSettings.type       = 'SLERP';