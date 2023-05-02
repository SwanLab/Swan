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
                    cParams.optimizerNames.aJmax = 2.5;
                    cParams.optimizerNames.aGmax = 0;
%                     cParams.optimizerNames.aJmax = optParams.aJmax;
%                     cParams.optimizerNames.aGmax = optParams.aGmax;
                    op = OptimizerNullSpace(cParams);
                otherwise
                    error('Invalid optimizer.')
            end
        end
        
    end
    
end