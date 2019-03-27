classdef OptimizerToDesignVariableTranslator <  handle
    
    methods (Access = public, Static)
        
        function n = translate(optimizer)
            switch optimizer
                case {'SLERP', 'PROJECTED SLERP', 'HAMILTON-JACOBI'}
                    n = 'LevelSet';
                case {'PROJECTED GRADIENT', 'MMA', 'IPOPT'}
                    n = 'Density';
            end
            
        end
        
    end
    
end

