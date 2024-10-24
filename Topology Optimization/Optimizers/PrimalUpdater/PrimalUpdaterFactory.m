classdef PrimalUpdaterFactory < handle
    
   
    methods (Access = public, Static)
    
        function p = create(cParams)
            switch cParams.primal
                case 'SLERP'
                    lsF = cParams.designVariable.obtainFunctionInCell();
                    s.mesh = lsF{1}.mesh;
                    p = SLERP(s);
                case 'PROJECTED GRADIENT'
                    p = ProjectedGradient(cParams);
                case 'HAMILTON-JACOBI'
                    p = HamiltonJacobi(cParams);
            end
        end
    
    end
    
end