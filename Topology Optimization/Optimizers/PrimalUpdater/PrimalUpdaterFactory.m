classdef PrimalUpdaterFactory < handle
    
   
    methods (Access = public, Static)
    
        function p = create(cParams)
            switch cParams.optimizerNames.primal
                case 'SLERP'
                    p = SLERP(cParams);
                case 'PROJECTED GRADIENT'
                    p = ProjectedGradient(cParams);
                case 'HAMILTON-JACOBI'
                    p = HamiltonJacobi(cParams);
            end
        end
    
    end
    
end