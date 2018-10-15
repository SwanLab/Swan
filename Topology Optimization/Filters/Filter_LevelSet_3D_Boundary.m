classdef Filter_LevelSet_3D_Boundary < Filter_LevelSet_3D & Filter_LevelSet_Boundary
    methods (Access = public)
        function createUnfittedMesh(obj)
            obj.unfitted_mesh = Mesh_Unfitted_3D_Boundary(obj.mesh.duplicate,obj.diffReacProb.geometry.interpolation);
        end
        
        function setInterpolation_Unfitted(obj)
            obj.interpolation_unfitted = Triangle_Linear(obj.unfitted_mesh);
        end
    end
    
    methods (Static, Access = public)
        function quadrature = getQuadrature_Unfitted
            quadrature = Quadrature_Triangle;
        end
    end
end

