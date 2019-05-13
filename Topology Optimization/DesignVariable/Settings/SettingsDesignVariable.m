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
            obj.initMesh();
            obj.initValue();
            obj.initLevelSetCreator();
        end
        
        function initMesh(obj)
            if ischar(obj.mesh)
                meshFile = obj.mesh;
                obj.mesh = Mesh_GiD(meshFile);
            end
        end
        
        function initValue(obj)
            obj.value = ones(size(obj.mesh.coord,1),1);
        end
        
        function initLevelSetCreator(obj)
            s.ndim  = obj.mesh.ndim;
            s.coord = obj.mesh.coord;
            s.type = obj.initialCase;
            obj.levelSetCreatorSettings = SettingsLevelSetCreator().create(s);
        end
        
    end
    
end