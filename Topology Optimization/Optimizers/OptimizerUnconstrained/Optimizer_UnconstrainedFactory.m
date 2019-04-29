classdef Optimizer_UnconstrainedFactory < handle
    
    methods (Access = public, Static)
        
        function op = create(cParams)
            switch cParams.optimizer
                case 'SLERP'
                    op = Optimizer_SLERP(cParams.uncOptimizerSettings);
                case 'HAMILTON-JACOBI'
                    op = Optimizer_HJ(cParams.uncOptimizerSettings,cParams.designVar);
                case 'PROJECTED GRADIENT'
                    op = Optimizer_PG(cParams.uncOptimizerSettings);     
                otherwise
                    error('Invalid Unconstrained optimizer.')
            end
        end
        
    end
    
end