classdef Filter_LevelSet_3D_Interior < Filter_LevelSet_3D & Filter_LevelSet_Interior
    methods (Access = public)
        function createUnfittedMesh(obj)
            obj.unfitted_mesh = Mesh_Unfitted_3D_Interior(obj.mesh.duplicate,obj.diffReacProb.geometry.interpolation);
        end
    end
end

