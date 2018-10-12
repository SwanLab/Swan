classdef Filter_LevelSet_2D < Filter_LevelSet
    methods
        function obj = Filter_LevelSet_2D
            obj.max_subcells = 6;
            obj.nnodes_subelem = 3;
            obj.ndim = 2;
        end
        
        function getQuadrature_Unfitted(obj)
            obj.quadrature_unfitted = Quadrature_Triangle;
        end
        
        function getInterpolation_Unfitted(obj)
            obj.interpolation_unfitted = Triangle_Linear(obj.unfitted_mesh);
        end
        
        function createUnfittedMesh_Interior(obj)
            obj.unfitted_mesh = Mesh_Unfitted_2D_Interior(obj.mesh.duplicate,obj.diffReacProb.geometry.interpolation);
        end
        
        function createUnfittedMesh_Boundary(obj)
            obj.unfitted_mesh = Mesh_Unfitted_2D_Boundary(obj.mesh.duplicate,obj.diffReacProb.geometry.interpolation);
        end

        function [interp_facet,quadrature_facet] = createFacet(obj)
            quadrature_facet = Quadrature.set('LINE');
            interp_facet = Line_Linear;
            quadrature_facet.computeQuadrature(obj.quadrature_fitted.order);
            interp_facet.computeShapeDeriv(quadrature_facet.posgp);
        end
    end
    
    methods (Static)        
        function djacob = mapping(points,dvolu)
            v = diff(points);
            L = norm(v);
            djacob = L/dvolu;
        end
    end
end

