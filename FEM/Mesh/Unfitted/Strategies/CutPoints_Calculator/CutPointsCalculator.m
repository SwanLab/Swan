classdef CutPointsCalculator < handle
    
    properties(GetAccess = public, SetAccess = private)
        cutPointsIso
        cutPointsGlobal
    end
        
    properties (Access = private)
        backgroundMesh
        backgroundCutCells
        backgroundGeomInterpolation
        levelSet  
        activeCutPoints
        index1
        index2
        lsNode1
        lsNode2
    end
    
    methods (Access = public)
        
        function init(obj,cParams)
            obj.backgroundMesh      = cParams.backgroundMesh;
            obj.levelSet            = cParams.levelSet_background;
            obj.backgroundCutCells  = cParams.backgroundCutCells;
            obj.backgroundGeomInterpolation =  Interpolation.create(obj.backgroundMesh.type,'LINEAR');
        end
        
        function computeCutPoints(obj)
            obj.computeIndex();
            obj.computeLevelSetInNodes();
            obj.computeActiveCutPoints();
            obj.computeCutPointsIso();
            obj.computeCutPointsGlobal();
        end
        
        function cutPoints = getThisCellCutPoints(obj,i)
            cutPoints.ISO = obj.getThisCellActiveCutPointsIso(i);
            cutPoints.GLOBAL = obj.getThisCellActiveCutPointsGlobal(i);
        end
        
    end
    
    methods (Access = private)
        
        function cP = getThisCellActiveCutPointsIso(obj,i)
            cutNode = obj.activeCutPoints(:,:,i);
            cP = obj.cutPointsIso(cutNode,:,i);
        end
        
        function cP = getThisCellActiveCutPointsGlobal(obj,i)
            cutNode = obj.activeCutPoints(:,:,i);
            cP = obj.cutPointsGlobal(cutNode,:,i);
        end
                

        function computeIndex(obj)
            node1 = obj.backgroundGeomInterpolation.iteration(1,:);
            node2 = obj.backgroundGeomInterpolation.iteration(2,:);
            connec        = obj.backgroundMesh.connec;
            cutCells      = obj.backgroundCutCells;
            cutNode1 = connec(cutCells,node1);
            cutNode2 = connec(cutCells,node2);
            obj.index1 = permute(cutNode1,[2 3 1]);
            obj.index2 = permute(cutNode2,[2 3 1]);
        end
        
      function [ls1,ls2] = computeLevelSetInNodes(obj)
            ls = obj.levelSet;
            ls1 = zeros(size(obj.index1));
            ls2 = zeros(size(obj.index2));
            ls1(:,:,:) = ls(obj.index1);
            ls2(:,:,:) = ls(obj.index2);
            obj.lsNode1 = ls1;
            obj.lsNode2 = ls2;
      end
            
        function computeActiveCutPoints(obj)
            ls1 = obj.lsNode1;
            ls2 = obj.lsNode2;
            obj.activeCutPoints = sign(ls1.*ls2)<=0;
        end
        
        function computeCutPointsIso(obj)
            [x1,x2] = obj.computeNodesIsoCoords();
            obj.cutPointsIso = obj.computeCutPoint(x1,x2);
        end
        
        function computeCutPointsGlobal(obj)
            [x1,x2] = obj.computeNodesGlobalCoords();
            obj.cutPointsGlobal = obj.computeCutPoint(x1,x2);
        end
        
        function xCut = computeCutPoint(obj,x1,x2)
            ls1 = obj.lsNode1;
            ls2 = obj.lsNode2;
            xCut = x1+ls1.*(x2-x1)./(ls1-ls2);
        end
        
        function [x1,x2] = computeNodesIsoCoords(obj)
            node1 = obj.backgroundGeomInterpolation.iteration(1,:);
            node2 = obj.backgroundGeomInterpolation.iteration(2,:);
            xNodes = obj.backgroundGeomInterpolation.pos_nodes;
            xNodes1 = xNodes(node1,:);
            xNodes2 = xNodes(node2,:);
            nCutCells = length(obj.backgroundCutCells);
            x1 = repmat(xNodes1,[1 1 nCutCells]);
            x2 = repmat(xNodes2,[1 1 nCutCells]);
        end
        
        function [x1,x2] = computeNodesGlobalCoords(obj)
            ndim  = obj.backgroundMesh.ndim;
            nPossibleCut = length(obj.backgroundGeomInterpolation.iteration(1,:));
            nCutCells = length(obj.backgroundCutCells);
            x1 = zeros(nPossibleCut,ndim,nCutCells);
            x2 = zeros(nPossibleCut,ndim,nCutCells);
            for idim = 1:ndim
                xNode = obj.backgroundMesh.coord(:,idim);
                xNode1 = xNode(obj.index1);
                xNode2 = xNode(obj.index2);
                x1(:,idim,:) = xNode1;
                x2(:,idim,:) = xNode2;
            end
        end
        
    end
    
end

