classdef PrimalUpdaterFactory < handle
    
   
    methods (Access = public, Static)
    
        function p = create(cParams)
            switch cParams.primal
                case 'SLERP'
                    s.mesh    = cParams.designVariable.fun.mesh;
                    s.filter  = cParams.filter;
                    p = SLERP(s);
                case 'PROJECTED GRADIENT'
                    p = ProjectedGradient(cParams);
                case 'HAMILTON-JACOBI'
                    p = HamiltonJacobi(cParams);
            end
        end
    
    end
    
end