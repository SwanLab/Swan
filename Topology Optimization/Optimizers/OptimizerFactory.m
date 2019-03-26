classdef OptimizerFactory < handle
    
    methods (Access = public, Static)
        
        function optimizer = create(optimizer,settings,designVar,epsilon)
            mesh = designVar.mesh;
            switch optimizer
                case 'SLERP'
                    unconstrainedOptimizer = Optimizer_SLERP(settings,epsilon);
                    optimizer = Optimizer_AugLag(settings,mesh,unconstrainedOptimizer);
                case 'HAMILTON-JACOBI'
                    unconstrainedOptimizer = Optimizer_HJ(settings,epsilon,designVar);
                    optimizer = Optimizer_AugLag(settings,mesh,unconstrainedOptimizer);
                case 'PROJECTED GRADIENT'
                    unconstrainedOptimizer = Optimizer_PG(settings,epsilon);
                    optimizer = Optimizer_AugLag(settings,mesh,unconstrainedOptimizer);
                case 'MMA'
                    optimizer = Optimizer_MMA(settings,mesh);
                case 'IPOPT'
                    optimizer = Optimizer_IPOPT(settings,mesh);
                case 'PROJECTED SLERP'
                    optimizer = Optimizer_Projected_Slerp(settings,mesh,epsilon);
                otherwise
                    error('Invalid optimizer.')
            end
        end
        
    end
    
end