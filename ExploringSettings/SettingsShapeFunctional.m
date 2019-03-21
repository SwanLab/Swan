classdef SettingsShapeFunctional < AbstractSettings
    
    properties
        filterParams = SettingsFilter()
        filename 
        ptype = 'MACRO'
    end
    
     methods (Access = public)
        
        function obj = SettingsShapeFunctional(varargin)
            if nargin == 1
                obj.loadParams(varargin{1})
            end
        end
        
    end
    
end