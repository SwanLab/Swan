classdef SettingsFilter < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsFilter'
    end
    
    properties (Access = public)
        filterType 
        domainType
        designVar
    end
    
    methods (Access = public)
        function obj = SettingsFilter(varargin)
            if nargin == 1
                    obj.loadParams(varargin{1});
            end
        end
    end
    
end