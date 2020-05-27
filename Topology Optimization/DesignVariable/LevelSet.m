classdef LevelSet < DesignVariable
    
    
    properties (Access = private)
        levelSetCreatorSettings
        unfittedMesh
    end    
    
    methods (Access = public)
        
        function obj = LevelSet(cParams)
            obj.nVariables = 1;            
            obj.init(cParams);
            obj.levelSetCreatorSettings = cParams.levelSetCreatorSettings;
            obj.createValue();
            obj.createUnfittedMesh();
        end
        
        function update(obj,value)
            obj.value = value;
            obj.updateUnfittedMesh();
        end
        
        function m = getUnfittedMesh(obj)
            m = obj.unfittedMesh;
        end
        
    end
    
    methods (Access = private)
        
        function createValue(obj)
            s = obj.levelSetCreatorSettings;
            s.ndim  = obj.mesh.ndim;
            s.coord = obj.mesh.coord;
            lsCreator  = LevelSetCreator.create(s);
            obj.value  = lsCreator.getValue();
        end        
        
        function createUnfittedMesh(obj)
            cParams = SettingsMeshUnfitted('INTERIOR',obj.mesh);
            obj.unfittedMesh = UnfittedMesh(cParams);
            obj.updateUnfittedMesh();
        end
        
        function updateUnfittedMesh(obj)
            obj.unfittedMesh.compute(obj.value);
        end
        
    end
    
end

