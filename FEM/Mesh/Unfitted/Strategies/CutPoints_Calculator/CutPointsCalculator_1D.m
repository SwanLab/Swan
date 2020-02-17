classdef CutPointsCalculator_1D < CutPointsCalculator_Abstract
   
    methods (Access = protected)
        
        function computeCutPoints_Iso(obj)
            obj.computeCutPointsInIsoCoordinates();
            obj.computeActiveCutPoints();
        end
        
        function computeCutPoints_Global(obj)
            [ls1,ls2] = obj.computeLevelSetInNodes2();
            [x1,x2] = obj.computeNodesIsoPosition2();            
            xCut = x1 + ls1.*(x2-x1)./(ls1-ls2);           
            obj.cutPointsGlobal = xCut;
            obj.activeCutPoints = sign(ls1.*ls2)<=0;
        end
        
    end
    
    methods (Access = private)
        
        function computeCutPointsInIsoCoordinates(obj)
            [ls1,ls2] = obj.computeLevelSetInNodes();
            [x1,x2]   = obj.computeNodesIsoPosition();                        
            xCut = x1 + ls1.*(x2-x1)./(ls1-ls2);            
            obj.cutPointsIso = xCut;            
        end
        
        function computeActiveCutPoints(obj)
            nCutCells = length(obj.backgroundCutCells);            
            obj.activeCutPoints = repmat([true;false],[1 1 nCutCells]);            
        end
        
        function [ls1,ls2] = computeLevelSetInNodes(obj)
            connec        = obj.meshBackground.connec;
            cutCells      = obj.backgroundCutCells;
            nodesCutCells = connec(cutCells,:);  
            ls            = obj.levelSet_background;
            lsCutNodes    = ls(nodesCutCells);
            ls1 = permute(lsCutNodes,[2 3 1]);
            ls2 = circshift(ls1,[-1 0 0]);                        
        end
        
        function [x1,x2] = computeNodesIsoPosition(obj)
            nCutCells = length(obj.backgroundCutCells);            
            xNodes = obj.backgroundGeomInterpolation.pos_nodes;                        
            x1 = repmat(xNodes,[1 1 nCutCells]);
            x2 = circshift(x1,[-1 0 0]);                        
        end
        
        function [ls1,ls2] = computeLevelSetInNodes2(obj)
            connec        = obj.meshBackground.connec;
            cutCells      = obj.backgroundCutCells;
            nodesCutCells = connec(cutCells,:);              
            index1 = permute(nodesCutCells,[2 3 1]);
            index2 = circshift(index1,[1 0 0]);            
            ls = obj.levelSet_background;
            ls1 = ls(index1);
            ls2 = ls(index2);                     
        end
        
        function [x1,x2] = computeNodesIsoPosition2(obj)
            connec        = obj.meshBackground.connec;
            cutCells      = obj.backgroundCutCells;
            nodesCutCells = connec(cutCells,:);              
            index1 = permute(nodesCutCells,[2 3 1]);
            index2 = circshift(index1,[1 0 0]);              
            coord = obj.meshBackground.coord(:,1);            
            x1 = coord(index1);
            x2 = coord(index2);            
        end
        
    end
    
end

