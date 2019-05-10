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
            if ischar(obj.mesh)
                obj.mesh = Mesh_GiD(obj.mesh);
            end
            obj.value = ones(size(obj.mesh.coord,1),1);
        end
        
    end
    
end