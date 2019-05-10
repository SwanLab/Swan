classdef SettingsFilter < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsFilter.json'
    end
    
    properties (Access = public)
        filterType 
        domainType
        designVar
        quadratureOrder
    end
    
    methods (Access = public)
        function obj = SettingsFilter(varargin)
            if nargin == 1
                    obj.loadParams(varargin{1});
            end
        end
    end
    
end