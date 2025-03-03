classdef Optimizer_UnconstrainedFactory < handle
    
    methods (Access = public, Static)
        
        function op = create(cParams)
            switch cParams.type
                case 'SLERP'
                    op = Optimizer_SLERP(cParams);
                case 'HAMILTON-JACOBI'
                    op = Optimizer_HJ(cParams);
                case 'PROJECTED GRADIENT'
                    op = Optimizer_PG(cParams);
                otherwise
                    error('Invalid Unconstrained optimizer.')
            end
        end
        
    end
    
end