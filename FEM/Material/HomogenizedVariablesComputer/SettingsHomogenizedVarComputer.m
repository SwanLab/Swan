classdef SettingsHomogenizedVarComputer < AbstractSettings
    
    properties (Access = public)
        type
        designVariable
        dim
    end
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            factory = SettingsHomogenizedVarComputerFactory();
            obj     = factory.create(cParams);
        end
        
    end
    
end