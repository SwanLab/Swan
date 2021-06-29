classdef SettingsIncrementalScheme < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsIncrementalScheme.json'
    end
    
    properties (Access = public)
        nSteps
        targetParamsSettings
        shallPrintIncremental
        mesh
    end
    
    methods (Access = public)
        
        function obj = SettingsIncrementalScheme(varargin)
            if nargin == 1
                obj.loadParams(varargin{1});
            end
            obj.init();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.createTargetParamsSettings();
        end
        
        function createTargetParamsSettings(obj)
            s = obj.targetParamsSettings;
            s.mesh = obj.mesh;
            s.nSteps = obj.nSteps;
            obj.targetParamsSettings = SettingsTargetParamsManager(s);
        end
        
        function setTargetParamsMesh(obj,m)
            obj.targetParamsSettings.mesh = m;
        end
        
    end
    
    methods
        
        function set.mesh(obj,m)
            obj.mesh = m;
            obj.setTargetParamsMesh(m);
        end
        
    end
    
end