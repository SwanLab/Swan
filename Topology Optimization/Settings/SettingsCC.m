classdef SettingsCC < AbstractSettings
    
    properties (Access = protected, Abstract)
        defaultParamsName
    end
    
    properties (Access = public)
        settings
        shapeFuncSettings
        nShapeFuncs        
        designVar
        homogenizedVarComputer
        targetParameters
    end
    
    methods (Access = public)
        
        function obj = SettingsCost(varargin)
            if nargin == 1
                obj.loadParams(varargin{1});
            end
        end
        
    end
    
    methods (Access = public, Static)
        
        function s = create(cParams,settings)

        end
        
    end
    
end