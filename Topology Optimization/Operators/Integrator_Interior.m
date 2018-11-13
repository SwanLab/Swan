classdef Integrator_Interior < Integrator
    %     properties (Access = private)
    %         shapeValues_FullCells
    %     end
    
    methods (Access = public)
        function A = computeIntegral(obj,F1)
            shapeValues_CutCells = obj.integrateCutCells(F1);
            shapeValues_FullCells = obj.integrateFullCells(F1); % !! shapeValues_FullCells could saved in a property instead of being re-computed all the time !!
            shapeValues_All = obj.assembleShapeValues(shapeValues_CutCells,shapeValues_FullCells);
            A = obj.rearrangeOutputRHS(shapeValues_All);
        end
    end
    
    methods (Access = private)
        function shapeValues_FullCells = integrateFullCells(obj,F1)
            % !! F1 should be evalutated in Integrator and integration at Interior Full
            % Cells should be allowed !!
            
            interpolation = Interpolation.create(obj.mesh_background,'LINEAR');
            quadrature = obj.computeQuadrature(obj.mesh_background.geometryType);
            interpolation.computeShapeDeriv(quadrature.posgp);
            geometry = Geometry(obj.mesh_background,'LINEAR');
            geometry.computeGeometry(quadrature,interpolation);
            
            shapeValues_FullCells = zeros(size(obj.mesh_background.connec));
            for igauss = 1:quadrature.ngaus
                shapeValues_FullCells = shapeValues_FullCells + interpolation.shape(:,igauss)'.*geometry.dvolu(:,igauss);
                %                 shapeValues_FullCells = shapeValues_FullCells + interpolation.shape(:,igauss)'.*F1.*geometry.dvolu(:,igauss);
            end
        end
        
        function shapeValues_AllCells = assembleShapeValues(obj,shapeValues_CutCells,shapeValues_FullCells)
            interpolation = Interpolation.create(obj.mesh_background,'LINEAR');
            shapeValues_AllCells = zeros(size(obj.mesh_background.connec));
            shapeValues_AllCells(obj.mesh_unfitted.background_full_cells,:) = shapeValues_FullCells(obj.mesh_unfitted.background_full_cells,:);
            
            for i_subcell = 1:size(shapeValues_CutCells,2)
                shapeValues_AllCells(:,i_subcell) = shapeValues_AllCells(:,i_subcell)+accumarray(obj.mesh_unfitted.cell_containing_subcell,shapeValues_CutCells(:,i_subcell),[interpolation.nelem,1],@sum,0);
            end
        end
    end
end

