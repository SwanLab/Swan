classdef Filter_LevelSet_2D_Boundary < Filter_LevelSet_2D & Filter_LevelSet_Boundary
    methods (Access = public)
        function createUnfittedMesh(obj)
            obj.unfitted_mesh = Mesh_Unfitted_2D_Boundary(obj.mesh.duplicate,obj.diffReacProb.geometry.interpolation);
        end
    end
end

