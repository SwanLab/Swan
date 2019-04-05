classdef LevelSet < DesignVariable
    
    properties (GetAccess = public, SetAccess = protected)
        type = 'LevelSet';
    end
    
    methods (Access = public)
        
        function obj = LevelSet(cParams)
            obj.value = cParams.value;
            obj.meshGiD = cParams.mesh;
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

