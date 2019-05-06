classdef OptimizerFactory < handle
    
    methods (Access = public, Static)
        
        function op = create(cParams)
            switch cParams.name
                case 'AlternatingPrimalDual'
                    op = OptimizerAlternatingPrimalDual(cParams);
                case 'MMA'
                    op = Optimizer_MMA(cParams);
                case 'IPOPT'
                    op = Optimizer_IPOPT(cParams);
                case 'OptimizerDualNestedInPrimal'
                    op = OptimizerDualNestedInPrimal(cParams);
                otherwise
                    error('Invalid optimizer.')
            end
        end
        
    end
    
end