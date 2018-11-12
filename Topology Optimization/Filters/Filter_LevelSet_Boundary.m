classdef Filter_LevelSet_Boundary < Filter_LevelSet
    properties (Access = public)
        domainType = 'BOUNDARY';
    end
    
    methods (Access = public)
        function preProcess(obj)
            preProcess@Filter_LevelSet(obj);
        end
        
        function M2 = computeRHS(obj,F1)           
            shapeValues = obj.integrateFoverMesh(F1);
            shapeValues = obj.assembleShapeValues(shapeValues);
            M2 = obj.rearrangeOutputRHS(shapeValues);
        end
    end
    
    methods (Access = private)
        function shapeValues_AllCells = assembleShapeValues(obj,shapeValues_CutCells)
            shapeValues_AllCells = zeros(size(obj.mesh.connec,1),size(obj.mesh.connec,2));

            for i_subcell = 1:size(shapeValues_CutCells,2)
                shapeValues_AllCells(:,i_subcell) = shapeValues_AllCells(:,i_subcell)+accumarray(obj.unfitted_mesh.cell_containing_subcell,shapeValues_CutCells(:,i_subcell),[obj.nelem,1],@sum,0);
            end
        end
    end
end