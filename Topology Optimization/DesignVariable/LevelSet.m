classdef LevelSet < DesignVariable
    
    properties (Access = private)
        levelSetCreatorSettings
    end    
    
    methods (Access = public)
        
        function obj = LevelSet(cParams)
            obj.init(cParams);
            obj.levelSetCreatorSettings = cParams.levelSetCreatorSettings;
            obj.createValue();
            obj.createUnfittedMesh();
        end
        
        function update(obj,value)
            obj.value = value;
            obj.updateUnfittedMesh();
        end
        
    end
    
    methods (Access = private)
        
        function createValue(obj)
            s = obj.levelSetCreatorSettings;
            lsCreator  = LevelSetCreator.create(s);
            obj.value  = lsCreator.getValue();
        end        
        
        function createUnfittedMesh(obj)
            unfittedSettings = SettingsMeshUnfitted('INTERIOR',obj.meshGiD);
            obj.mesh = Mesh_Unfitted(unfittedSettings);
            obj.updateUnfittedMesh();
        end
        
        function updateUnfittedMesh(obj)
            obj.mesh.computeMesh(obj.value);
        end
        
    end
    
end

