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
        
    end
    
    methods (Access = private)
        
        function createValue(obj)
            s = obj.levelSetCreatorSettings;
            lsCreator  = LevelSetCreator.create(s);
            obj.value  = lsCreator.getValue();
        end        
        
        function createUnfittedMesh(obj)
            unfittedSettings = SettingsMeshUnfitted('INTERIOR',obj.mesh);
            obj.unfittedMesh = Mesh_Unfitted(unfittedSettings);
            obj.updateUnfittedMesh();
        end
        
        function updateUnfittedMesh(obj)
            obj.unfittedMesh.computeMesh(obj.value);
        end
        
    end
    
end

