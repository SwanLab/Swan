classdef WeakFormSolverFactory < handle
    
    methods (Access = public, Static)  
        function wf = create(cParams)
            switch cParams.stab
                case 1
                    wf = ConvDifGalerkinSystem(cParams);
                case 2
                    wf = ConvDifSUSystem(cParams);
                case 3
                    wf = Optimizer_IPOPT(cParams);
                case 4
                    wf = OptimizerBisection(cParams);
                case 5
                    wf = Optimizer_fmincon(cParams);
            end
        end
    end
end