classdef CutPointsCalculator < handle
    
    properties(GetAccess = public, SetAccess = private)
        cutPointsIso
        cutPointsGlobal
        activeCutPoints
    end
        
    properties (Access = private)
        meshBackground
        levelSet_background
        backgroundCutCells
        backgroundGeomInterpolation
        index1
        index2
    end
    
    methods (Access = public)
        
        function init(obj,mesh)
            obj.meshBackground = mesh.meshBackground;
            obj.levelSet_background = mesh.levelSet_background;
            obj.backgroundCutCells = mesh.backgroundCutCells;
            obj.backgroundGeomInterpolation = mesh.backgroundGeomInterpolation;
        end
        
        function computeCutPoints(obj)
            obj.computeCutPoints_Iso();
            obj.computeCutPoints_Global();
        end
        
        function cutPoints = getThisCellCutPoints(obj,i)
            cutPoints.ISO = obj.getThisCellActiveCutPointsIso(i);
            cutPoints.GLOBAL = obj.getThisCellActiveCutPointsGlobal(i);
        end
        
    end
    
    methods (Access = private)
        
        function cP = getThisCellActiveCutPointsIso(obj,i)
            cP = obj.cutPointsIso(obj.activeCutPoints(:,:,i),:,i);
        end
        
        function cP = getThisCellActiveCutPointsGlobal(obj,i)
            cP = obj.cutPointsGlobal(obj.activeCutPoints(:,:,i),:,i);
        end
        
        function cP = getThisCellActiveCutPoints(obj,cutPoints,i)
            cP = cutPoints(obj.activeCutPoints(:,:,i),:,i);
        end
                
        function [x1,x2] = computeNodesIsoPosition(obj)
            iteration1 = obj.backgroundGeomInterpolation.iteration(1,:);
            iteration2 = obj.backgroundGeomInterpolation.iteration(2,:);
            xNodes = obj.backgroundGeomInterpolation.pos_nodes;
            xNodes1 = xNodes(iteration1,:);
            xNodes2 = xNodes(iteration2,:);
            nCutCells = length(obj.backgroundCutCells);
            x1 = repmat(xNodes1,[1 1 nCutCells]);
            x2 = repmat(xNodes2,[1 1 nCutCells]);
        end
        
        function computeIndex(obj)
            iteration1 = obj.backgroundGeomInterpolation.iteration(1,:);
            iteration2 = obj.backgroundGeomInterpolation.iteration(2,:);
            connec        = obj.meshBackground.connec;
            cutCells      = obj.backgroundCutCells;
            nodesCutCells1 = connec(cutCells,iteration1);
            nodesCutCells2 = connec(cutCells,iteration2);
            obj.index1 = permute(nodesCutCells1,[2 3 1]);
            obj.index2 = permute(nodesCutCells2,[2 3 1]);
        end
        
        function [x1,x2] = computeNodesIsoPosition2(obj)
            ndim  = obj.meshBackground.ndim;
            nPossibleCut = length(obj.backgroundGeomInterpolation.iteration(1,:));
            nCutCells = length(obj.backgroundCutCells);
            x1 = zeros(nPossibleCut,ndim,nCutCells);
            x2 = zeros(nPossibleCut,ndim,nCutCells);
            for idim = 1:ndim
                xNode = obj.meshBackground.coord(:,idim);
                xNode1 = xNode(obj.index1);
                xNode2 = xNode(obj.index2);
                x1(:,idim,:) = xNode1;
                x2(:,idim,:) = xNode2;
            end
        end
        
        function [ls1,ls2] = computeLevelSetInNodes2(obj)
            ls = obj.levelSet_background;
            ls1 = ls(obj.index1);
            ls2 = ls(obj.index2);
        end
        
        function [ls1,ls2] = computeLevelSetInNodes(obj)
            iteration1 = obj.backgroundGeomInterpolation.iteration(1,:);
            iteration2 = obj.backgroundGeomInterpolation.iteration(2,:);
            ls = obj.levelSet_background;
            
            cutCells = obj.backgroundCutCells;
            connec   = obj.meshBackground.connec(cutCells,:);
            nodesCutCells1 = connec(:,iteration1);
            nodesCutCells2 = connec(:,iteration2);
            lsNodes1 = ls(nodesCutCells1);
            lsNodes2 = ls(nodesCutCells2);
            ls1 = permute(lsNodes1,[2 3 1]);
            ls2 = permute(lsNodes2,[2 3 1]);
        end
        
        function computeCutPoints_Iso(obj)
            [ls1,ls2] = obj.computeLevelSetInNodes();
            [x1,x2]   = obj.computeNodesIsoPosition();
            xCut = x1 + ls1.*(x2-x1)./(ls1-ls2);
            obj.cutPointsIso = xCut;
            obj.activeCutPoints = sign(ls1.*ls2)<=0;
        end
        
        function computeCutPoints_Global(obj)
            obj.computeIndex();
            [ls1,ls2] = obj.computeLevelSetInNodes2();
            [ls1,ls2] = obj.computeLevelSetInNodes();
            [x1,x2]   = obj.computeNodesIsoPosition2();
            %ndime = size(x1,2);
            %for idime = 1:ndime
            %    x1v = squeeze(x1(:,idime,:));
            %    x2v = squeeze(x2(:,idime,:));
            %    xCut(:,idime,:) = x1v +ls1.*(x2v - x1v)./(ls1-ls2);
            %end
            xCut = x1+ls1.*(x2-x1)./(ls1-ls2);
            obj.cutPointsGlobal = xCut;
           % obj.activeCutPoints = sign(ls1.*ls2)<=0;
        end
        
        
        
        
    end
    
end

