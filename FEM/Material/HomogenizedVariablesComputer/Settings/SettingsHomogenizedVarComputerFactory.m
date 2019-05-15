classdef SettingsHomogenizedVarComputerFactory < handle
   
    methods (Access = public, Static)
        
        function obj = create(cParams)
            switch cParams.type
                case 'ByInterpolation'
                    obj = SettingsHomogenizedVarComputerFromInterpolation(cParams);
                case 'ByVademecum'
                    obj = SettingsHomogenizedVarComputerFromVademecum(cParams);
                otherwise
                    error('Invalid HomogenizerVarComputer type.');
            end
        end
        
    end
    
end