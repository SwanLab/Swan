classdef Integrator_Interior < Integrator_Composed
    %     properties (Access = private)
    %         shapeValues_FullCells
    %     end
    
    methods (Access = public)
        function obj = Integrator_Interior(mesh)
            obj@Integrator_Composed(mesh);
        end
        
        function A = integrate(obj,F)
            A_list = integrate@Integrator_Composed(obj,F);
            if ~isempty(A_list)
                A = summatory(A_list);
            end
        end
    end

    methods (Static, Access = private)
        function A = summatory(A_list)
            A = zeros(size(A_list{1}));
            for imesh = 1:obj.n_meshes
                A = A + A_list{imesh};
            end
        end
    end
    
%     methods (Access = private)
%         function shapeValues_FullCells = integrateFullCells(obj,F1)
%             % !! F1 should be evalutated in Integrator and integration at Interior Full
%             % Cells should be allowed !!
%             
%             interpolation = Interpolation.create(obj.mesh_background,'LINEAR');
%             quadrature = obj.computeQuadrature(obj.mesh_background.geometryType);
%             interpolation.computeShapeDeriv(quadrature.posgp);
%             geometry = Geometry(obj.mesh_background,'LINEAR');
%             geometry.computeGeometry(quadrature,interpolation);
%             
%             shapeValues_FullCells = zeros(size(obj.mesh_background.connec));
%             for igauss = 1:quadrature.ngaus
%                 shapeValues_FullCells = shapeValues_FullCells + interpolation.shape(:,igauss)'.*geometry.dvolu(:,igauss);
%                 %                 shapeValues_FullCells = shapeValues_FullCells + interpolation.shape(:,igauss)'.*F1.*geometry.dvolu(:,igauss);
%             end
%         end
%     end
end

