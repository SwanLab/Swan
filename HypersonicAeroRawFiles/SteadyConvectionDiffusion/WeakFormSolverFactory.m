classdef WeakFormSolverFactory < handle
    
    methods (Access = public, Static)  
        function op = create(cParams)
            switch cParams.stab
                case 1
                    op = OptimizerAugmentedLagrangian(cParams);
                case 2
                    op = OptimizerMMA(cParams);
                case 3
                    op = Optimizer_IPOPT(cParams);
                case 4
                    op = OptimizerBisection(cParams);
                case 5
                    op = Optimizer_fmincon(cParams);
            end
        end
    end
end