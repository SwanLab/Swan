% classdef MasterSlaveComputer < handle
% 
%     properties (Access = private)
%         tol = 10e-6
%         coord
%         div
%         nodes
%         normalVec
%         coord2A
%         ortoNodes
%         pairs
%     end
% 
%     properties (Access = public)
%         masterSlaveIndex
%     end
% 
%     methods (Access = public)
% 
%         function obj = MasterSlaveComputer(cParams)
%             obj.init(cParams);
%         end
% 
%         function computeMasterSlaveNodes(obj)
%             obj.computeSidesNormalVector();
%             obj.findVectorCandidates();
%             obj.findOrtogonalVectors();
%             obj.computeSidesPariety();
%             obj.createMasterAndSlaves();
%         end
% 
%     end
% 
%     methods (Access = private)
% 
%         function init(obj,cParams)
%             obj.coord = cParams.coord;
%             obj.div = cParams.div;
%             obj.nodes = cParams.nodes;
%         end
% 
%         function computeSidesNormalVector(obj)
%             obj.normalVec = zeros(obj.nodes.vert,2);
%             for iVert = 1:obj.nodes.vert
%                 [vertexA,vertexB] = obj.obtainVertices(iVert);
%                 vecU = MasterSlaveComputer.computeUnitaryVector(vertexA,vertexB);
%                 obj.computeNormalVec(iVert,vecU);
%             end
%         end
% 
%         function findVectorCandidates(obj)
%             bound = obj.coord(1:obj.nodes.bound,:);
%             vert = obj.coord(1:obj.nodes.vert,:);
%             nIntNodes = sum(2*(obj.div-1));
%             intNodes = bound(obj.nodes.vert+1:end,:);
%             obj.coord2A = zeros(nIntNodes,2,obj.nodes.vert);
%             for iVert = 1:obj.nodes.vert
%                 vertexA = vert(iVert,:);
%                 for iVec = 1:nIntNodes
%                     vertexC = intNodes(iVec,:);
%                     vecU = MasterSlaveComputer.computeUnitaryVector(vertexA,vertexC);
%                     obj.coord2A(iVec,:,iVert) = obj.coord2A(iVec,:,iVert)+vecU;
%                 end
%             end
%         end
% 
%         function findOrtogonalVectors(obj)
%             nIntNodes = sum(2*(obj.div-1));
%             obj.ortoNodes = zeros(max(obj.div)-1,obj.nodes.vert);
%             for iVert = 1:obj.nodes.vert
%                 iPos = 1;
%                 vectorA = obj.normalVec(iVert,:);
%                 for iVec = 1:nIntNodes
%                     vectorB = obj.coord2A(iVec,:,iVert);
%                     if abs(dot(vectorA,vectorB)) < obj.tol
%                         nodeNumber = obj.nodes.vert+iVec;
%                         obj.ortoNodes(iPos,iVert) = obj.ortoNodes(iPos,iVert)+nodeNumber;
%                         iPos = iPos+1;
%                     end
%                 end
%                 obj.turnNodes(iVert);
%             end
%         end
% 
%         function computeSidesPariety(obj)
%             obj.pairs = zeros(obj.nodes.vert/2,2);
%             irow = 1;
%             for iSrch = 1:obj.nodes.vert/2
%                 jSrch = iSrch+1;
%                 vectorA = obj.normalVec(iSrch,:);
%                 found = 0;
%                 while found == 0
%                     vectorB = obj.normalVec(jSrch,:);
%                     cosAngle = abs(dot(vectorA,vectorB));
%                     if abs(cosAngle-1) < obj.tol
%                         obj.pairs(irow,1) = obj.pairs(irow,1)+iSrch;
%                         obj.pairs(irow,2) = obj.pairs(irow,2)+jSrch;
%                         irow = irow+1;
%                         found = 1;
%                     else
%                         jSrch = jSrch+1;
%                     end
%                 end
%             end
%         end
% 
%         function createMasterAndSlaves(obj)
%             obj.masterSlaveIndex = zeros((obj.nodes.bound-obj.nodes.vert)/2,2);
%             cont = 1;
%             for iPair = 1:obj.nodes.vert/2
%                 lineA = obj.pairs(iPair,1);
%                 lineB = obj.pairs(iPair,2);
%                 for iDiv = 1:obj.div(iPair)-1
%                     obj.masterSlaveIndex(cont,1) = obj.masterSlaveIndex(cont,1)+obj.ortoNodes(iDiv,lineA);
%                     obj.masterSlaveIndex(cont,2) = obj.masterSlaveIndex(cont,2)+obj.ortoNodes(iDiv,lineB);
%                     cont = cont+1;
%                 end
%             end
%         end
% 
%         function [vertexA,vertexB] = obtainVertices(obj,iVert)
%             vertexA = obj.coord(iVert,:);
%             if iVert == obj.nodes.vert
%                 vertexB = obj.coord(1,:);
%             else
%                 vertexB = obj.coord(iVert+1,:);
%             end
%         end
% 
%         function computeNormalVec(obj,iVert,vecU)
%             obj.normalVec(iVert,1) = obj.normalVec(iVert,1)+vecU(1,2);
%             obj.normalVec(iVert,2) = obj.normalVec(iVert,2)-vecU(1,1);
%         end
% 
% %         function checkOrtogonality(obj,vectorA,vectorB,iVec,iVert,iPos)
% % 
% %         end
% 
%         function turnNodes(obj,iVert)
%             if iVert > obj.nodes.vert/2
%                 obj.ortoNodes(:,iVert) = sort(obj.ortoNodes(:,iVert),'descend');
%             end
%         end
% 
%     end
% 
%     methods (Static)
% 
%         function v_norm = computeUnitaryVector(A,B)
%             v = B-A;
%             m = norm(v);
%             v_norm = v/m;
%         end
% 
%     end
% 
% end
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