classdef CutPointsCalculator_3D < CutPointsCalculator_Abstract
    
    methods (Access = protected)
        
        function computeCutPoints_Iso(obj)
            pos_nodes = obj.backgroundGeomInterpolation.pos_nodes;
            
            iteration1 = obj.backgroundGeomInterpolation.iteration(1,:);
            iteration2 = obj.backgroundGeomInterpolation.iteration(2,:);
            
            gamma1 = permute(obj.levelSet_background(obj.meshBackground.connec(obj.backgroundCutCells,iteration1)),[2 3 1]);
            gamma2 = permute(obj.levelSet_background(obj.meshBackground.connec(obj.backgroundCutCells,iteration2)),[2 3 1]);
            
            P1 = repmat(pos_nodes(iteration1,:),[1 1 size(obj.backgroundCutCells)]);
            P2 = repmat(pos_nodes(iteration2,:),[1 1 size(obj.backgroundCutCells)]);
            P = P1 + gamma1.*(P2-P1)./(gamma1-gamma2);
            
            obj.cutPointsIso = P;
            obj.activeCutPoints = sign(gamma1.*gamma2)<=0;
        end
        
        function computeCutPoints_Global(obj)
            iteration1 = obj.backgroundGeomInterpolation.iteration(1,:);
            iteration2 = obj.backgroundGeomInterpolation.iteration(2,:);
            
            index1 = permute(obj.meshBackground.connec(obj.backgroundCutCells,iteration1),[2 3 1]);
            index2 = permute(obj.meshBackground.connec(obj.backgroundCutCells,iteration2),[2 3 1]);
            
            gamma1 = obj.levelSet_background(index1);
            gamma2 = obj.levelSet_background(index2);
            
            coord1 = obj.meshBackground.coord(:,1);
            coord2 = obj.meshBackground.coord(:,2);
            coord3 = obj.meshBackground.coord(:,3);
            
            P1 = [coord1(index1) coord2(index1) coord3(index1)];
            P2 = [coord1(index2) coord2(index2) coord3(index2)];
            P = P1+gamma1.*(P2-P1)./(gamma1-gamma2);
            
            obj.cutPointsGlobal = P;
            obj.activeCutPoints = sign(gamma1.*gamma2)<=0;
        end
        
    end
    
end

