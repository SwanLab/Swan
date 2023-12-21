classdef OptimizerFactory < handle
    
    methods (Access = public, Static)
        
        function op = create(cParams)
            switch cParams.type
                case 'AlternatingPrimalDual'
                    op = OptimizerAugmentedLagrangian(cParams);
                case 'MMA'
                    op = Optimizer_MMA(cParams);
                case 'IPOPT'
                    op = Optimizer_IPOPT(cParams);
                case 'DualNestedInPrimal'
                    op = OptimizerBisection(cParams);
                case 'fmincon'
                    op = Optimizer_fmincon(cParams);
                case 'NullSpace'
                    %aJmax --> constrain, weight to decrease cost
                    %aGmax --> constrain, weight constaint value
                    %cParams.optimizerNames.aJmax = 8; % 4 
                    %aJ previous: 10e-04
                    cParams.optimizerNames.aJmax = 2.5; % earlier calibrated 2.5
                    cParams.optimizerNames.aGmax = 0.15; % earlier calibrated 0.15
                    %cParams.optimizerNames.aGmax = 0.0;
                    op = OptimizerNullSpace(cParams);
                otherwise
                    error('Invalid optimizer.')
            end
        end
        
    end
    
end