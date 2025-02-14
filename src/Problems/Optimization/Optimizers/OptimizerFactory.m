classdef OptimizerFactory < handle
    
    methods (Access = public, Static)
        
        function op = create(cParams)
            switch cParams.type
                case 'AlternatingPrimalDual'
                    op = OptimizerAugmentedLagrangian(cParams);
                case 'MMA'
                    op = OptimizerMMA(cParams);
                case 'IPOPT'
                    op = Optimizer_IPOPT(cParams);
                case 'DualNestedInPrimal'
                    op = OptimizerBisection(cParams);
                case 'fmincon'
                    op = Optimizerfmincon(cParams);
                case 'NullSpace'
%                     cParams.optimizerNames.aJmax = 2;
%                     cParams.optimizerNames.aGmax = 0.05;
                    op = OptimizerNullSpace(cParams);
                case 'IPM'
                    op = OptimizerInteriorPoint(cParams);
                otherwise
                    error('Invalid optimizer.')
            end
        end
        
    end
    
end