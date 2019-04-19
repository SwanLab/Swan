classdef LevelSet < DesignVariable
    
    methods (Access = public)
        
        function obj = LevelSet(cParams)
            obj.init(cParams);
            obj.createUnfittedMesh(cParams.mesh);
        end
        
        function update(obj,value)
            obj.value = value;
            obj.updateUnfittedMesh();
        end
        
    end
    
    methods (Access = private)
        
        function createUnfittedMesh(obj,mesh)
            unfittedSettings = SettingsMeshUnfitted('INTERIOR',mesh);
            obj.mesh = Mesh_Unfitted(unfittedSettings);
            obj.updateUnfittedMesh();
        end
        
        function updateUnfittedMesh(obj)
            obj.mesh.computeMesh(obj.value);
        end
        
    end
    
end

