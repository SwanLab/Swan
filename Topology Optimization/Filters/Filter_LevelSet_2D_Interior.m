classdef Filter_LevelSet_2D_Interior < Filter_LevelSet_2D & Filter_LevelSet_Interior
    methods (Access = public)
        function createUnfittedMesh(obj)
            obj.unfitted_mesh = Mesh_Unfitted_2D_Interior(obj.mesh.duplicate,obj.diffReacProb.geometry.interpolation);
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

