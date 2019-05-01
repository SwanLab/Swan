problemData.problemFileName = 'Cantilever_triangle_coarse';

designVarSettings             = SettingsDesignVariable();
designVarSettings.type        = 'LevelSet';
designVarSettings.initialCase = 'full';

incrementalSchemeSettings = SettingsIncrementalScheme();
incrementalSchemeSettings.nSteps = 1;

optimizerSettings            = SettingsOptimizer();
optimizerSettings.name       = 'AugmentedLagrangian';
optimizerSettings.shallPrint = true;
optimizerSettings.settingsMonitor.showOptParams       = true;
optimizerSettings.settingsMonitor.refreshInterval     = 2;
optimizerSettings.settingsMonitor.shallDisplayDesignVar       = true;
optimizerSettings.settingsMonitor.shallShowBoundaryConditions = true;

optimizerSettings.uncOptimizerSettings = SettingsOptimizerUnconstrained();
optimizerSettings.uncOptimizerSettings.type       = 'SLERP';