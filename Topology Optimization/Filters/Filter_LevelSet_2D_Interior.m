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
        function djacob = mapping(points,dvolu)
            v1 = [diff(points([1 2],:)) 0];
            v2 = [diff(points([1 3],:)) 0];
            A = 0.5*norm(cross(v1,v2));
            djacob = A/dvolu;
        end
        
        function quadrature = getQuadrature_Unfitted
            quadrature = Quadrature_Triangle;
        end
    end
end

