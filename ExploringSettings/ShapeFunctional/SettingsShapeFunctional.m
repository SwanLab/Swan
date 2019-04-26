classdef SettingsShapeFunctional < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsShapeFunctional'
    end
    
    properties (Access = public)
        filterParams 
        filename 
        domainType 
       homogVarComputer
    end
    
     methods (Access = public)
        
        function obj = SettingsShapeFunctional(varargin)
            if nargin == 1
                    obj.loadParams(varargin{1});
            end
        end        
     end
end