classdef Integrator_Interior < IntegratorUnfitted  
    
    methods (Access = public)
       
    function A = computeIntegral(obj,F1)
            if obj.isLeveSetCuttingMesh()
                shapeValues_CutCells = obj.integrateCutCells(F1);
                shapeValues_FullCells = obj.integrateFullCells(F1);
                shapeValues_All = obj.assembleShapeValues(shapeValues_CutCells,shapeValues_FullCells);
            else
                shapeValues_All = obj.integrateFullCells(F1);
            end
            
            A = obj.rearrangeOutputRHS(shapeValues_All);
        end        
        
    end
    
    methods (Access = private)
        
        function shapeValues_FullCells = integrateFullCells(obj,F1)            
            interpolation = Interpolation.create(obj.meshBackground,'LINEAR');
            quadrature = obj.computeQuadrature(obj.meshBackground.geometryType);
            interpolation.computeShapeDeriv(quadrature.posgp);
            geometry = Geometry(obj.meshBackground,'LINEAR');
            geometry.computeGeometry(quadrature,interpolation);
            
            shapeValues_FullCells = zeros(size(obj.meshBackground.connec));
            for igauss = 1:quadrature.ngaus
                shapeValues_FullCells = shapeValues_FullCells + interpolation.shape(:,igauss)'.*geometry.dvolu(:,igauss);
            end
        end
        
        function shapeValues_AllCells = assembleShapeValues(obj,shapeValues_CutCells,shapeValues_FullCells)
            interpolation = Interpolation.create(obj.meshBackground,'LINEAR');
            shapeValues_AllCells = zeros(size(obj.meshBackground.connec));
            shapeValues_AllCells(obj.meshUnfitted.backgroundFullCells,:) = shapeValues_FullCells(obj.meshUnfitted.backgroundFullCells,:);
            
            for i_subcell = 1:size(shapeValues_CutCells,2)
                shapeValues_AllCells(:,i_subcell) = shapeValues_AllCells(:,i_subcell)+accumarray(obj.meshUnfitted.cellContainingSubcell,shapeValues_CutCells(:,i_subcell),[interpolation.nelem,1],@sum,0);
            end
        end
        
    end
    
end

