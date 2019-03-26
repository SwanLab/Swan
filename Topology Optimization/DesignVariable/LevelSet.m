classdef LevelSet < DesignVariable
    
    methods (Access = public)
        
        function obj = LevelSet(cParams)
            obj.createUnfittedMesh(cParams.mesh);
            obj.updateUnfittedMesh(cParams.value);
        end
        
    end
    
    methods (Access = private)
        
        function createUnfittedMesh(obj,mesh)
            unfittedSettings = SettingsMeshUnfitted('INTERIOR',mesh);
            obj.mesh = Mesh_Unfitted(unfittedSettings);
        end
        
        function updateUnfittedMesh(obj,phi)
           obj.value = phi;
           obj.mesh.computeMesh(phi);
        end
        
    end
    
end

