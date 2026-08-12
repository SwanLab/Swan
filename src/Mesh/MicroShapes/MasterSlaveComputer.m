classdef MasterSlaveComputer < handle
    
    properties (Access = private)
        coord
        nodes
        div
        latticeVectors
        tol = 1e-10
    end
    
    properties (Access = public)
        masterSlaveIndex
    end
    
    methods (Access = public)
        
        function obj = MasterSlaveComputer(cParams)
            obj.init(cParams);
            obj.computeMasterSlaveNodes();
        end
        
        function computeMasterSlaveNodes(obj)
            nVectors = size(obj.latticeVectors, 1);
            
            
            if nVectors == 2
                obj.computeMasterSlaveParallelogram();
            else
                obj.computeMasterSlaveGeneric();
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.coord = cParams.coord;
            obj.nodes = cParams.nodes;
            obj.div = cParams.div;
            obj.latticeVectors = cParams.latticeVectors;
        end
        
        function computeMasterSlaveParallelogram(obj)
            nDiv1 = obj.div(1);
            nDiv2 = obj.div(2); 
            base = obj.nodes.vert;  
            startL1 = base + 1;
            startL2 = startL1 + (nDiv1 - 1); 
            startL3 = startL2 + (nDiv2 - 1);
            startL4 = startL3 + (nDiv1 - 1);
            nPairs = (nDiv1 - 1) + (nDiv2 - 1);
            obj.masterSlaveIndex = zeros(nPairs, 2);
            idx = 1;
            for i = 1:nDiv1-1
                master = startL1 + (i - 1);
                slave  = startL3 + (nDiv1 - 1 - i); 
                obj.masterSlaveIndex(idx, :) = [master, slave];
                idx = idx + 1;
            end
            for i = 1:nDiv2-1
                master = startL2 + (i - 1);
                slave  = startL4 + (nDiv2 - 1 - i);  
                obj.masterSlaveIndex(idx, :) = [master, slave];
                idx = idx + 1;
            end
        end
        
        function computeMasterSlaveGeneric(obj)
            boundaryNodes = 1:obj.nodes.bound;
            obj.masterSlaveIndex = obj.createPeriodicPairs(boundaryNodes);
        end
        function boundaryNodes = findBoundaryNodes(obj)
            boundaryNodes = 1:obj.nodes.bound;
        end
        
        function pairs = createPeriodicPairs(obj, boundaryNodes)
            v1 = obj.latticeVectors(1, :);
            v2 = obj.latticeVectors(2, :);
            
            nBoundary = length(boundaryNodes);
            pairs = zeros(nBoundary, 2);
            pairCount = 0;
            
            for i = 1:nBoundary
                nodeIdx = boundaryNodes(i);
                coordNode = obj.coord(nodeIdx, :);
                target1 = coordNode + v1;
                j1 = obj.findClosestNode(target1, boundaryNodes);
                target2 = coordNode + v2;
                j2 = obj.findClosestNode(target2, boundaryNodes);
                
                if ~isempty(j1) && j1 > nodeIdx
                    pairCount = pairCount + 1;
                    pairs(pairCount, :) = [nodeIdx, j1];
                elseif ~isempty(j2) && j2 > nodeIdx
                    pairCount = pairCount + 1;
                    pairs(pairCount, :) = [nodeIdx, j2];
                end
            end
            pairs = pairs(1:pairCount, :);
        end
        
        function idx = findClosestNode(obj, target, candidateNodes)
            minDist = inf;
            idx = [];
            
            for i = 1:length(candidateNodes)
                nodeIdx = candidateNodes(i);
                dist = norm(obj.coord(nodeIdx, :) - target);
                if dist < obj.tol
                    idx = nodeIdx;
                    return;
                elseif dist < minDist
                    minDist = dist;
                    idx = nodeIdx;
                end
            end
            
            if minDist > 1e-6
                idx = [];
            end
        end
        
    end
    
end