classdef OptimizerFactory < handle
    
    methods (Access = public, Static)
        
        function optimizer = create(settings)
            
            optimizerName = settings.optimizer;
            designVar = settings.designVar;
            switch optimizerName
                case 'SLERP'
                    settings.unconstrainedOptimizer = Optimizer_SLERP(settings.uncOptimizerSettings);
                    settings.settings        = settings;
                    settings.designVariable  = settings.designVar;
                    settings.monitoring      = settings.monitoring;
                    settings.convergenceVars = settings.unconstrainedOptimizer.convergenceVars;
                    optimizer = Optimizer_AugLag(settings);
                case 'HAMILTON-JACOBI'
                    settings.unconstrainedOptimizer = Optimizer_HJ(settings.uncOptimizerSettings,designVar);
                    settings.settings        = settings;
                    settings.designVariable  = settings.designVar;
                    settings.monitoring      = settings.monitoring;
                    settings.convergenceVars = settings.unconstrainedOptimizer.convergenceVars;
                    optimizer = Optimizer_AugLag(settings);
                case 'PROJECTED GRADIENT'
                    settings.unconstrainedOptimizer = Optimizer_PG(settings.uncOptimizerSettings);
                    settings.settings        = settings;
                    settings.designVariable  = settings.designVar;
                    settings.monitoring      = settings.monitoring;
                    settings.convergenceVars = settings.unconstrainedOptimizer.convergenceVars;
                    optimizer = Optimizer_AugLag(settings);
                case 'MMA'
                    settings.settings        = settings;
                    settings.designVariable  = settings.designVar;
                    settings.monitoring      = settings.monitoring;
                    settings.convergenceVars = ConvergenceVariables(2);
                    optimizer = Optimizer_MMA(settings);
                case 'IPOPT'
                    settings.settings        = settings;
                    settings.designVariable  = settings.designVar;
                    settings.monitoring      = settings.monitoring;
                    settings.convergenceVars = ConvergenceVariables(1);
                    optimizer = Optimizer_IPOPT(settings);
                case 'PROJECTED SLERP'
                    settings.settings        = settings;
                    settings.designVariable  = settings.designVar;
                    settings.monitoring      = settings.monitoring;
                    settings.convergenceVars = ConvergenceVariables(3);
                    optimizer = Optimizer_Projected_Slerp(settings);
                otherwise
                    error('Invalid optimizer.')
            end
        end
        
    end
    
end