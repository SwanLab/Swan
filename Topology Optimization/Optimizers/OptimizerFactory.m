classdef OptimizerFactory < handle
    
    methods (Access = public, Static)
        
        function optimizer = create(settings)
            
            optimizerName = settings.optimizer;
            designVar = settings.designVar;
            switch optimizerName
                case 'SLERP'
                    settings.optimizer_unconstr = Optimizer_SLERP(settings.uncOptimizerSettings);
                    optimizer = Optimizer_AugLag(settings);
                case 'HAMILTON-JACOBI'
                    settings.optimizer_unconstr = Optimizer_HJ(settings.uncOptimizerSettings,designVar);
                    optimizer = Optimizer_AugLag(settings);
                case 'PROJECTED GRADIENT'
                    settings.optimizer_unconstr = Optimizer_PG(settings.uncOptimizerSettings);
                    optimizer = Optimizer_AugLag(settings);
                case 'MMA'
                    optimizer = Optimizer_MMA(settings);
                case 'IPOPT'
                    optimizer = Optimizer_IPOPT(settings);
                case 'PROJECTED SLERP'
                    optimizer = Optimizer_Projected_Slerp(settings);
                otherwise
                    error('Invalid optimizer.')
            end
        end
        
    end
    
end