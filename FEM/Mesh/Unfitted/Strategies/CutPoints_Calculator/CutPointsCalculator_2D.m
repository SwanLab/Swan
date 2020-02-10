classdef CutPointsCalculator_2D < CutPointsCalculator_Abstract
    
    methods (Access = protected)
        
        function computeCutPoints_Iso(obj)
            pos_nodes = obj.backgroundGeomInterpolation.pos_nodes;
            
            ls = obj.levelSet_background;
            
            ls1 = permute(ls(obj.meshBackground.connec(obj.backgroundCutCells,:)),[2 3 1]);
            ls2 = circshift(ls1,[-1 0 0]);
            
            P1 = repmat(pos_nodes,[1 1 size(obj.backgroundCutCells)]);
            P2 = circshift(P1,[-1 0 0]);
            P = P1 + ls1.*(P2-P1)./(ls1-ls2);
            
            obj.cutPointsIso = P;
            obj.activeCutPoints = sign(ls1.*ls2) <= 0;
        end
        
        function computeCutPoints_Global(obj)
            index1 = permute(obj.meshBackground.connec(obj.backgroundCutCells,:),[2 3 1]);
            index2 = circshift(index1,[-1 0 0]);
            
            ls = obj.levelSet_background;
           
            ls1 = ls(index1);
            ls2 = ls(index2);
            
            coord1 = obj.meshBackground.coord(:,1);
            coord2 = obj.meshBackground.coord(:,2);
            
            P1 = [coord1(index1) coord2(index1)];
            P2 = [coord1(index2) coord2(index2)];
            P = P1 + ls1.*(P2-P1)./(ls1-ls2);
            
            obj.cutPointsGlobal = P;
            obj.activeCutPoints = sign(ls1.*ls2) <= 0;       
        end
        
    end
    
end

