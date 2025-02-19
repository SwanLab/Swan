classdef HomogenizedVarComputerFactory < handle
    
    methods (Access = public, Static)
        
        function h = create(cParams)
            switch cParams.type
                case 'ByVademecum'
                    h = HomogenizedVarComputerFromVademecum(cParams);
                case 'ByInterpolation'
                    h = HomogenizedVarComputerFromInterpolation(cParams);
            end
            
        end
        
    end
    
end