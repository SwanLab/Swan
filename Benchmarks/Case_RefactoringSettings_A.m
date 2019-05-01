fileName = 'Cantilever_triangle_coarse';

designVarSettings             = SettingsDesignVariable();
designVarSettings.type        = 'Density';
designVarSettings.initialCase = 'full';

incrementalSchemeSettings = SettingsIncrementalScheme();
incrementalSchemeSettings.nSteps = 1;

optimizerSettings            = SettingsOptimizer();
optimizerSettings.name       = 'IPOPT';
optimizerSettings.shallPrint = false;
optimizerSettings.settingsMonitor.showOptParams       = false;
optimizerSettings.settingsMonitor.refreshInterval     = 2;
optimizerSettings.settingsMonitor.shallDisplayDesignVar       = true;
optimizerSettings.settingsMonitor.shallShowBoundaryConditions = true;