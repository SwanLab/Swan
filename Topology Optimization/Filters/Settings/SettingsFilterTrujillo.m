classdef SettingsFilter < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsFilter'
    end
    
    properties (Access = public)
        filename
        mesh
        scale
    end
    
    methods (Access = public)
        
        function obj = SettingsFilter(varargin)
            switch nargin
                case 1
                    obj.loadParams(varargin{1});
            end
        end
        
    end
    
end