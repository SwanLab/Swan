classdef OptimizerFactory < handle
    
    methods (Access = public, Static)
        
        function op = create(cParams,optParams)
            switch cParams.type
                case 'AlternatingPrimalDual'
                    op = OptimizerAugmentedLagrangian(cParams,optParams);
                case 'MMA'
                    op = Optimizer_MMA(cParams);
                case 'IPOPT'
                    op = Optimizer_IPOPT(cParams);
                case 'DualNestedInPrimal'
                    op = OptimizerBisection(cParams);
                case 'fmincon'
                    op = Optimizer_fmincon(cParams);
                case 'NullSpace'
                    cParams.optimizerNames.aJmax = 1.7; % 1.7, 3
                    cParams.optimizerNames.aGmax = 0.05; % 0.05,  0.3
                    op = OptimizerNullSpace(cParams);
                otherwise
                    error('Invalid optimizer.')
            end
        end
        
    end
    
end