classdef Integrator < handle
    properties (GetAccess = public, SetAccess = private)
        mesh
    end
    
    properties (Access = private)
        output_size
        shape_size
    end
    
    methods (Static, Access = public)
        function obj = create(mesh)
            factory = IntegratorFactory;
            obj = factory.create(mesh);
        end
    end
    
    methods (Access = public)
        function A = integrateMesh(obj,F,mesh)
            if exist('mesh_unfitted','var')
                obj.saveMesh(mesh);
            end
            A = obj.computeIntegral(F);
        end
    end
    
    methods (Access = private)
        function A = computeIntegral(obj,F1,mesh)
            shapeValues = obj.integrate(F1);
            shapeValues = obj.assembleShapeValues(shapeValues,mesh);
            A = obj.rearrangeOutputRHS(shapeValues,mesh);
        end
        
        function shapeValues = assembleShapeValues(obj,shapeValues_CutCells,mesh)
            interpolation = Interpolation.create(mesh,'LINEAR');
            shapeValues = zeros(size(mesh.connec));
            
            for i_subcell = 1:size(shapeValues_CutCells,2)
                shapeValues(:,i_subcell) = shapeValues(:,i_subcell)+accumarray(obj.mesh_unfitted.cell_containing_subcell,shapeValues_CutCells(:,i_subcell),[interpolation.nelem,1],@sum,0);
            end
        end
    end
    
    methods (Access = protected)
        function saveMesh(obj,mesh)
            obj.mesh = mesh;
        end
    end
    
    methods (Static, Access = protected)
        function M2 = rearrangeOutputRHS(shapeValues,mesh)
            interpolation = Interpolation.create(mesh,'LINEAR');
            
            M2 = zeros(interpolation.npnod,1);
            for inode = 1:interpolation.nnode
                M2 = M2 + accumarray(mesh.connec(:,inode),shapeValues(:,inode),[interpolation.npnod,1],@sum,0);
            end
        end
        
        function quadrature = computeQuadrature(geometryType)
            quadrature = Quadrature.set(geometryType);
            quadrature.computeQuadrature('LINEAR');
        end
        
        function posgp = computePosGP(subcell_coord,interpolation,quadrature)
            interpolation.computeShapeDeriv(quadrature.posgp);
            posgp = zeros(quadrature.ngaus,size(subcell_coord,3),size(subcell_coord,1));
            for igaus = 1:quadrature.ngaus
                for idime = 1:size(subcell_coord,3)
                    posgp(igaus,idime,:) = subcell_coord(:,:,idime)*interpolation.shape(:,igaus);
                end
            end
        end
    end
end

