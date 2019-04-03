classdef MasterSlaveRelator < handle
    
    properties (Access = private)
        x
        y
        nodesInXmin
        nodesInXmax
        nodesInYmin
        nodesInYmax
        cellDescriptor
        msRelation
        allNodes
    end
    
    methods (Access = public)
        
        function obj = MasterSlaveRelator(coord)
            obj.init(coord);
            obj.computeNodesInFaces(coord)
            obj.computeMasterSlaveRelation();
        end
        
        function r = getRelation(obj)
            r = obj.msRelation;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,coord)
            obj.x = coord(:,1);
            obj.y = coord(:,2);
            obj.allNodes(:,1) = 1:size(obj.x,1);
        end
                
        function computeNodesInFaces(obj,coord)
            cD = CellNodesDescriptor(coord);                        
            corners = cD.cornerNodes;
            obj.nodesInXmin = setdiff(cD.nodesInXmin,corners);
            obj.nodesInXmax = setdiff(cD.nodesInXmax,corners);
            obj.nodesInYmin = setdiff(cD.nodesInYmin,corners);
            obj.nodesInYmax = setdiff(cD.nodesInYmax,corners);
        end
               
        function computeMasterSlaveRelation(obj)
            [masterFaceX,slaveFaceX] = obj.computeMasterSlaveNodesInFaceX();
            [masterFaceY,slaveFaceY] = obj.computeMasterSlaveNodesInFaceY();
            obj.msRelation = [masterFaceX slaveFaceX;
                masterFaceY slaveFaceY];
        end
        
        function [master,slave] = computeMasterSlaveNodesInFaceX(obj)
            [master,slave] = obj.computeMasterSlaveNode(obj.nodesInXmin,obj.nodesInXmax,obj.y);
        end
        
        function [master,slave] = computeMasterSlaveNodesInFaceY(obj)
            [master,slave] = obj.computeMasterSlaveNode(obj.nodesInYmin,obj.nodesInYmax,obj.x);
        end
        
        function [master,slave] = computeMasterSlaveNode(obj,isLowerFace,isUpperFace,dir2compare)
            nLF = obj.allNodes(isLowerFace,1);
            nUF = obj.allNodes(isUpperFace,1);
            master = nLF;
            slave  = obj.obtainClosestNodeInOppositeFace(nLF,nUF,dir2compare);
        end
        
        function closestNodesUF = obtainClosestNodeInOppositeFace(obj,nodesLF,nodesUF,pos)
            nNodes         = length(nodesLF);
            closestNodesUF = zeros(nNodes,1);
            for inode = 1:nNodes
                nodeLF = nodesLF(inode);
                closestNodesUF(inode) = obj.obtainClosestNode(nodeLF,nodesUF,pos);
            end
        end
        
    end
    
    methods (Access = private, Static)
        
        function closestNodeB = obtainClosestNode(nodeA,nodeB,pos)
            xA = pos(nodeA,1);
            xB = pos(nodeB,1);
            distAB = abs(xA - xB);
            [distMin,isClosest] = min(distAB);
            if distMin > 1e-10
                error('non slave node')
            end
            closestNodeB = nodeB(isClosest);
        end
        
    end
    
    
end