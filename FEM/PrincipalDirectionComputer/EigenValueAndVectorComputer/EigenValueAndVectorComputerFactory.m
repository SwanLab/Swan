classdef EigenValueAndVectorComputerFactory < handle
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            switch cParams.type
                case 'SYMBOLIC'
                    obj = EigenValueAndVectorComputerSymbolic(cParams);
                case 'PRECOMPUTED'
                    obj = EigenValueAndVectorComputerPrecomputed(cParams);
            end
        end
        
    end
    
end