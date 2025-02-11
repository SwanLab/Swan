classdef SettingsDesignVariable < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsDesignVariable.json'
    end
    
    properties (Access = public)
        fun
        mesh
        type
        initialCase
        creatorSettings
        scalarProductSettings
        femData
        isFixed
    end
    
    methods (Access = public)
        
        function obj = SettingsDesignVariable(varargin)
            if nargin == 1
                obj.loadParams(varargin{1});
            end
            obj.init();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            switch obj.type
                case 'Density'
                    if isfield(obj.creatorSettings,'rho0') 
                        obj.creatorSettings.type = 'Given';
                    else
                        obj.creatorSettings.type = 'FromLevelSet';
                    end
                case 'LevelSet'
                    obj.initLevelSetCreator();
            end
        end
        
        function initLevelSetCreator(obj)
            s = obj.creatorSettings;
            s.type = obj.initialCase;
            obj.creatorSettings = SettingsLevelSetCreator().create(s);
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
        end
        
    end
    
end