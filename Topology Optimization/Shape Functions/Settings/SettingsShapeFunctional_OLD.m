classdef SettingsShapeFunctional_OLD < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsShapeFunctional.json'
    end
    
    properties (Access = public)
        type
        filterParams
        filename
        scale
        homogVarComputer
        designVariable
        targetParameters
    end
    
    methods (Access = public)
        
        function obj = SettingsShapeFunctional_OLD(varargin)
            if nargin == 1
                obj.loadParams(varargin{1});
            end
        end
        
    end
    
    methods (Access = public, Static)
        
        function s = create(cParams,settings)
            f = SettingsShapeFunctionalFactory_OLD();
            s = f.create(cParams,settings);
        end
        
    end
    
end