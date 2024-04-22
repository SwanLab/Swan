classdef WeakFormSolverFactory < handle
    
    methods (Access = public, Static)  
        function wf = create(cParams)
            switch cParams.stab
                case 'Galerkin'
                    wf = ConvDifGalerkinSystem(cParams);
                case 'Upwind'
                    wf = ConvDifSUSystem(cParams);
                case 'SUPG'
                    wf = ConvDifSUPGSystem(cParams);
            end
        end
    end
end