classdef OptimizerFactory < handle
    
    methods (Access = public, Static)
        
        function optimizer = create(cParams)
            
            optimizerName = cParams.optimizer;
            designVar = cParams.designVar;
            switch optimizerName
                case 'SLERP'
                    cParams.unconstrainedOptimizer = Optimizer_SLERP(cParams.uncOptimizerSettings);
                    cParams.convergenceVars = cParams.unconstrainedOptimizer.convergenceVars;
                    optimizer = Optimizer_AugLag(cParams);
                case 'HAMILTON-JACOBI'
                    cParams.unconstrainedOptimizer = Optimizer_HJ(cParams.uncOptimizerSettings,designVar);
                    cParams.convergenceVars = cParams.unconstrainedOptimizer.convergenceVars;
                    optimizer = Optimizer_AugLag(cParams);
                case 'PROJECTED GRADIENT'
                    cParams.unconstrainedOptimizer = Optimizer_PG(cParams.uncOptimizerSettings);
                    cParams.convergenceVars = cParams.unconstrainedOptimizer.convergenceVars;
                    optimizer = Optimizer_AugLag(cParams);
                case 'MMA'
                    cParams.convergenceVars = ConvergenceVariables(2);
                    optimizer = Optimizer_MMA(cParams);
                case 'IPOPT'
                    cParams.convergenceVars = ConvergenceVariables(1);
                    optimizer = Optimizer_IPOPT(cParams);
                case 'PROJECTED SLERP'
                    cParams.convergenceVars = ConvergenceVariables(3);
                    optimizer = Optimizer_Projected_Slerp(cParams);
                otherwise
                    error('Invalid optimizer.')
            end
        end
        
    end
    
end