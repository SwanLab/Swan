classdef SettingsShapeFunctional < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsShapeFunctional'
    end
    
    properties
        filterParams 
        filename 
        domainType 
        materialInterpolationParams
    end
    
     methods (Access = public)
        
        function obj = SettingsShapeFunctional(varargin)
            switch nargin
                case 1
                    obj.loadParams(varargin{1});
                case 2
                    obj.filterParams = varargin{1};
                    obj.materialInterpolationParams = varargin{2};

            end
        end        
     end
end