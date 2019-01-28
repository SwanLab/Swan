classdef CutPointsCalculator_1D < CutPointsCalculator_Abstract
   
    methods (Access = protected)
        
        function computeCutPoints_Iso(obj)
            pos_nodes = obj.backgroundGeomInterpolation.pos_nodes;
            
            gamma1 = permute(obj.levelSet_background(obj.meshBackground.connec(obj.backgroundCutCells,:)),[2 3 1]);
            gamma2 = circshift(gamma1,[-1 0 0]);
            
            P1 = repmat(pos_nodes,[1 1 size(obj.backgroundCutCells)]);
            P2 = circshift(P1,[-1 0 0]);
            P = P1 + gamma1.*(P2-P1)./(gamma1-gamma2);
            
            obj.cutPointsIso = P;
            obj.activeCutPoints = repmat([true;false],[1 1 size(obj.backgroundCutCells)]);
        end
        
        function computeCutPoints_Global(obj)
            index1 = permute(obj.meshBackground.connec(obj.backgroundCutCells,:),[2 3 1]);
            index2 = circshift(index1,[1 0 0]);
            
            gamma1 = obj.levelSet_background(index1);
            gamma2 = obj.levelSet_background(index2);
            
            coord1 = obj.meshBackground.coord(:,1);
            
            P1 = coord1(index1);
            P2 = coord1(index2);
            P = P1 + gamma1.*(P2-P1)./(gamma1-gamma2);
            
            obj.cutPointsGlobal = P;
            obj.activeCutPoints = sign(gamma1.*gamma2)<=0;
        end
        
    end
    
end

