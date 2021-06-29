classdef SettingsVideoManager < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsVideoManager.json'
    end
    
    properties (Access = public)
        caseFileName
        shallPrint
        designVarType
        pdim
    end
    
    methods (Access = public)
        
        function obj = SettingsVideoManager(varargin)
            if nargin == 1
                obj.loadParams(varargin{1});
            end
        end
        
    end
    
end