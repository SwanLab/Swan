classdef OptimizerFactory < handle
    
    methods (Access = public, Static)
        
        function optimizer = create(optimizerName,settings,designVar,epsilon)
            switch optimizerName
                case 'SLERP'
                    unconstrainedOptimizer = Optimizer_SLERP(settings,epsilon);
                    optimizer = Optimizer_AugLag(settings,designVar,unconstrainedOptimizer);
                case 'HAMILTON-JACOBI'
                    unconstrainedOptimizer = Optimizer_HJ(settings,epsilon,designVar);
                    optimizer = Optimizer_AugLag(settings,designVar,unconstrainedOptimizer);
                case 'PROJECTED GRADIENT'
                    unconstrainedOptimizer = Optimizer_PG(settings,epsilon);
                    optimizer = Optimizer_AugLag(settings,designVar,unconstrainedOptimizer);
                case 'MMA'
                    optimizer = Optimizer_MMA(settings,designVar);
                case 'IPOPT'
                    optimizer = Optimizer_IPOPT(settings,designVar);
                case 'PROJECTED SLERP'
                    optimizer = Optimizer_Projected_Slerp(settings,designVar,epsilon);
                otherwise
                    error('Invalid optimizer.')
            end
        end
        
    end
    
end