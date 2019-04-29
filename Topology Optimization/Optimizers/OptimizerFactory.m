classdef OptimizerFactory < handle
    
    methods (Access = public, Static)
        
        function op = create(cParams)
            switch cParams.optimizer
                case {'SLERP','HAMILTON-JACOBI','PROJECTED GRADIENT'}
                    op = Optimizer_AugLag(cParams);
                case 'MMA'
                    op = Optimizer_MMA(cParams);
                case 'IPOPT'
                    op = Optimizer_IPOPT(cParams);
                case 'PROJECTED SLERP'
                    op = Optimizer_Projected_Slerp(cParams);
                otherwise
                    error('Invalid optimizer.')
            end
        end
        
    end
    
end