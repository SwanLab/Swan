classdef WeakFormSolverFactory < handle
    
    methods (Access = public, Static)  
        function wf = create(cParams)
            switch cParams.stab
                case 1
                    wf = ConvDifGalerkinSystem(cParams);
                case 2
                    wf = ConvDifSUSystem(cParams);
                case 3
                    wf = ConvDifSUPGSystem(cParams);
            end
        end
    end
end