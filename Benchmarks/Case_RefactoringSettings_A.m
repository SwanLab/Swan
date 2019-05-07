problemData.problemFileName = 'Cantilever_triangle_coarse';
problemData.scale = 'MACRO';

designVarSettings             = SettingsDesignVariable();
designVarSettings.type        = 'LevelSet';
designVarSettings.initialCase = 'full';

incrementalSchemeSettings = SettingsIncrementalScheme();
incrementalSchemeSettings.nSteps = 1;

costSettings.shapeFuncSettings = {
    struct('type','compliance')%,...
%     struct('type','perimeter')
    };
costSettings.weights = [1 0.1];

constraintSettings.shapeFuncSettings = {struct('type','volume')};

optimizerSettings            = SettingsOptimizer();
optimizerSettings.name       = 'AlternatingPrimalDual';
optimizerSettings.shallPrint = true;
optimizerSettings.settingsMonitor.showOptParams       = true;
optimizerSettings.settingsMonitor.refreshInterval     = 2;
optimizerSettings.settingsMonitor.shallDisplayDesignVar       = true;
optimizerSettings.settingsMonitor.shallShowBoundaryConditions = true;

optimizerSettings.uncOptimizerSettings = SettingsOptimizerUnconstrained();
optimizerSettings.uncOptimizerSettings.type       = 'SLERP';