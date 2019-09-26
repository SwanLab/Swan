classdef PlottingLevelSetVigergauz < handle
    
    properties (Access = private)
        mesh
        levelSet
    end
    
    methods (Access = public)
        
        function obj = PlottingLevelSetVigergauz()
            obj.init();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.createMesh();
            obj.createLevelSet();
        end
        
        function createLevelSet(obj)
            cParams.type = 'LevelSet';
            cParams.mesh = obj.mesh;
            cParams.levelSetCreatorSettings = SettingsLevelSetCreator;
            cParams.scalarProductSettings = [];
            obj.levelSet = LevelSet(cParams);
        end
        
        function createMesh(obj)
            eval('RVE_Square_Triangle_FineFine');
            coord = coord(:,2:end);
            mesh = Mesh();
            obj.mesh = mesh.create(coord,connec);          
        end
        
    end
    
    
end