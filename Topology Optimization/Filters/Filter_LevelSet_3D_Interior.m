classdef Filter_LevelSet_3D_Interior < Filter_LevelSet_3D & Filter_LevelSet_Interior
    methods (Access = public)
        function createUnfittedMesh(obj)
            obj.unfitted_mesh = Mesh_Unfitted_3D_Interior(obj.mesh.duplicate,obj.diffReacProb.geometry.interpolation);
        end
        
        function setInterpolation_Unfitted(obj)
            obj.interpolation_unfitted = Tetrahedra_Linear(obj.unfitted_mesh);
        end
    end
    
    methods (Static, Access = public)
        function djacob = mapping(points,dvolu)
            v1 = diff(points([1 2],:));
            v2 = diff(points([1 3],:));
            v3 = diff(points([1 4],:));
            V = (1/6)*det([v1;v2;v3]);
            djacob = V/dvolu;
        end
        
        function quadrature = getQuadrature_Unfitted
            quadrature = Quadrature_Tetrahedra;
        end
    end
end

