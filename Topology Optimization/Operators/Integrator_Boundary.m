classdef Integrator_Boundary < Integrator
    
    methods (Access = public)
        
    end
    
    methods (Access = ?Integrator)
        function A = computeIntegral(obj,F1)
            shapeValues = obj.integrateCutCells(F1);
            shapeValues = obj.assembleShapeValues(shapeValues);
            A = obj.rearrangeOutputRHS(shapeValues);
        end
        
        function shapeValues_AllCells = assembleShapeValues(obj,shapeValues_CutCells)
            interpolation = Interpolation.create(obj.meshBackground,'LINEAR');
            shapeValues_AllCells = zeros(size(obj.meshBackground.connec));
            
            for i_subcell = 1:size(shapeValues_CutCells,2)
                shapeValues_AllCells(:,i_subcell) = shapeValues_AllCells(:,i_subcell)+accumarray(obj.mesh_unfitted.cell_containing_subcell,shapeValues_CutCells(:,i_subcell),[interpolation.nelem,1],@sum,0);
            end
        end
    end
end

