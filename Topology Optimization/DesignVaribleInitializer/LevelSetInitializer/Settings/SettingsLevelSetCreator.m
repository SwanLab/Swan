classdef SettingsLevelSetCreator < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsLevelSetCreator.json'
    end
    
    properties (Access = public)
        type
        ndim
        coord
        value
        geomParams
    end
    
    methods (Access = public)
        
        function obj = SettingsLevelSetCreator(varargin)
            if nargin == 1
                obj.loadParams(varargin{1})
            end
            obj.init();
        end
        
        function obj = create(obj,s)
           factory = SettingsLevelSetCreatorFactory();
           obj     = factory.create(s);
        end
        
    end
    
     methods (Access = private)
        
        function init(obj)
            obj.initCoord();
            obj.initValue();
        end
        
        function initCoord(obj)
            if ischar(obj.coord)
                meshFile = obj.coord;
                mesh = Mesh_GiD(meshFile);
                obj.coord = mesh.coord;
            end
        end
        
        function initValue(obj)
            obj.value = ones(size(obj.coord,1),1);
        end

    end
    
end