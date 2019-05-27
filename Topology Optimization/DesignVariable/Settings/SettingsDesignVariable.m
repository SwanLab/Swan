classdef SettingsDesignVariable < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsDesignVariable.json'
    end
    
    properties (Access = public)
        value
        mesh
        type
        initialCase
        levelSetCreatorSettings
        scalarProductSettings
        femData
    end
    
    methods (Access = public)
        
        function obj = SettingsDesignVariable(varargin)
            if nargin == 1
                obj.loadParams(varargin{1});
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.initLevelSetCreator();
        end
        
        function initLevelSetCreator(obj)
            s = obj.levelSetCreatorSettings;
            s.type = obj.initialCase;
            obj.levelSetCreatorSettings = SettingsLevelSetCreator().create(s);
        end
        
        function initScalarProductSettings(obj)
            obj.scalarProductSettings.femSettings.fileName = obj.femData.fileName;
            obj.scalarProductSettings.scale = obj.femData.scale;
        end
        
    end
    
    methods
        
        function set.femData(obj,pD)
            obj.femData = pD;
            obj.initScalarProductSettings();
        end
        
        function set.mesh(obj,m)
            obj.mesh = m;
            obj.init();
        end
        
    end
    
end