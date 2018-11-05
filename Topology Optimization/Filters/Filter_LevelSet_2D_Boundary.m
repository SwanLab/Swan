classdef Filter_LevelSet_2D_Boundary < Filter_LevelSet_2D & Filter_LevelSet_Boundary
    methods (Access = public)
        function createUnfittedMesh(obj)
            obj.unfitted_mesh = Mesh_Unfitted_2D_Boundary(obj.mesh.duplicate,obj.diffReacProb.geometry.interpolation);
        end
        
        function setInterpolation_Unfitted(obj)
            obj.interpolation_unfitted = Line_Linear;
        end
    end
    
    methods (Static, Access = public)
        function djacob = mapping(points,dvolu)
            v = diff(points);
            L = norm(v);
            djacob = L/dvolu;
        end
        
        function quadrature = getQuadrature_Unfitted
            quadrature = Quadrature_Line;
        end
    end
end

