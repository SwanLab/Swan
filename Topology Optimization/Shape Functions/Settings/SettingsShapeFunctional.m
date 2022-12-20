classdef SettingsShapeFunctional < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsShapeFunctional.json'
    end
    
    properties (Access = public)
        type
        filterParams
        femSettings
        homogVarComputer
        designVariable
        targetParameters
        mesh
        domainNotOptimizable
        shNumber
    end
    
    methods (Access = public)
        
        function obj = SettingsShapeFunctional(varargin)
            if nargin == 1
                obj.loadParams(varargin{1});
            end
        end
        
    end
    
    methods (Access = public, Static)
        
        function s = create(cParams)
            f = SettingsShapeFunctionalFactory();
            s = f.create(cParams);
        end
        
    end
    
end