classdef MasterSlaveRelator < handle
    
    properties (Access = private)
        x
        y
        z
        nodesInXmin
        nodesInXmax
        nodesInYmin
        nodesInYmax
        nodesInZmin
        nodesInZmax
        cellDescriptor
        msRelation
        allNodes
        cornerNodes 
        commonNodesEdges
    end
    
    methods (Access = public)
        
        function obj = MasterSlaveRelator(coord)
            switch size(coord,2)
                case 2
                    obj.init2(coord);
                    obj.computeNodesInFaces(coord);
                    obj.computeMasterSlaveRelation();
                case 3
                    obj.init3(coord);
                    obj.computeNodesInFaces3D(coord);
                    obj.computeMasterSlaveRelation3D();
%                     obj.avoidRepeatedEdges();
                    obj.avoidEdges();
            end
        end
        
        function r = getRelation(obj)
            r = obj.msRelation;
        end

        function e = getRepeatedEdges(obj)
            e = obj.commonNodesEdges;
        end
        
    end
    
    methods (Access = private)
        
        function init2(obj,coord)
            obj.x = coord(:,1);
            obj.y = coord(:,2);
            obj.allNodes(:,1) = 1:size(obj.x,1);
        end

        function init3(obj,coord)
            obj.x = coord(:,1);
            obj.y = coord(:,2);
            obj.z = coord(:,3);
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

        function computeNodesInFaces3D(obj,coord)
            cD = CellNodesDescriptor(coord);
            corners = cD.cornerNodes;
            obj.nodesInXmin = setdiff(cD.nodesInXmin,corners);
            obj.nodesInXmax = setdiff(cD.nodesInXmax,corners);
            obj.nodesInYmin = setdiff(cD.nodesInYmin,corners);
            obj.nodesInYmax = setdiff(cD.nodesInYmax,corners);
            obj.nodesInZmin = setdiff(cD.nodesInZmin,corners);
            obj.nodesInZmax = setdiff(cD.nodesInZmax,corners);
        end
        
        function computeMasterSlaveRelation(obj)
            [masterFaceX,slaveFaceX] = obj.computeMasterSlaveNodesInFaceX();
            [masterFaceY,slaveFaceY] = obj.computeMasterSlaveNodesInFaceY();
            obj.msRelation = [masterFaceX slaveFaceX;
                masterFaceY slaveFaceY];
        end

        function computeMasterSlaveRelation3D(obj)
            % not checked
            [masterFaceX,slaveFaceX] = obj.computeMasterSlaveNodesInFaceX3D();
            [masterFaceY,slaveFaceY] = obj.computeMasterSlaveNodesInFaceY3D();
            [masterFaceZ,slaveFaceZ] = obj.computeMasterSlaveNodesInFaceZ3D();
            obj.msRelation = [masterFaceX slaveFaceX;
                masterFaceY slaveFaceY;
                masterFaceZ slaveFaceZ];
        end
        
        function [master,slave] = computeMasterSlaveNodesInFaceX(obj)
            [master,slave] = obj.computeMasterSlaveNode(obj.nodesInXmin,obj.nodesInXmax,obj.y);
        end
        
        function [master,slave] = computeMasterSlaveNodesInFaceY(obj)
            [master,slave] = obj.computeMasterSlaveNode(obj.nodesInYmin,obj.nodesInYmax,obj.x);
        end

        function [master,slave] = computeMasterSlaveNodesInFaceX3D(obj)
            [master,slave] = obj.computeMasterSlaveNode3D(obj.nodesInXmin,obj.nodesInXmax,obj.y,obj.z);
        end
        
        function [master,slave] = computeMasterSlaveNodesInFaceY3D(obj)
            [master,slave] = obj.computeMasterSlaveNode3D(obj.nodesInYmin,obj.nodesInYmax,obj.x,obj.z);
        end
        
        function [master,slave] = computeMasterSlaveNodesInFaceZ3D(obj)
            [master,slave] = obj.computeMasterSlaveNode3D(obj.nodesInZmin,obj.nodesInZmax,obj.x,obj.y);
        end
        
        function [master,slave] = computeMasterSlaveNode(obj,isLowerFace,isUpperFace,dir2compare)
            nLF = obj.allNodes(isLowerFace,1);
            nUF = obj.allNodes(isUpperFace,1);
            master = nLF;
            slave  = obj.obtainClosestNodeInOppositeFace(nLF,nUF,dir2compare);
        end
        
        function [master,slave] = computeMasterSlaveNode3D(obj,isLowerFace,isUpperFace,dir2compare1,dir2compare2)
            nLF = obj.allNodes(isLowerFace,1);
            nUF = obj.allNodes(isUpperFace,1);
            master = nLF;
            slave  = obj.obtainClosestNodeInOppositeFace3D(nLF,nUF,dir2compare1,dir2compare2);
        end

        function closestNodesUF = obtainClosestNodeInOppositeFace(obj,nodesLF,nodesUF,pos)
            nNodes         = length(nodesLF);
            closestNodesUF = zeros(nNodes,1);
            for inode = 1:nNodes
                nodeLF = nodesLF(inode);
                closestNodesUF(inode) = obj.obtainClosestNode(nodeLF,nodesUF,pos);
            end
        end
        
        function closestNodesUF = obtainClosestNodeInOppositeFace3D(obj,nodesLF,nodesUF,pos1,pos2)
            nNodes         = length(nodesLF);
            closestNodesUF = zeros(nNodes,1);
            for inode = 1:nNodes
                nodeLF = nodesLF(inode);
                closestNodesUF(inode) = obj.obtainClosestNode3D(nodeLF,nodesUF,pos1,pos2);
            end
        end

        function avoidRepeatedEdges(obj)
            MS = obj.msRelation;
            uniqueCol1 = unique(MS(:,1));
            uniqueCol2 = unique(MS(:,2));
            commonNodes = intersect(uniqueCol1, uniqueCol2);
            rowsRemoved = ismember(MS(:,1), commonNodes) | ismember(MS(:,2), commonNodes);
            MS(rowsRemoved, :) = [];
            obj.msRelation = MS;
            obj.commonNodesEdges = commonNodes;
        end

        function avoidEdges(obj)
            MS = obj.msRelation;
            uniqueValues = unique(MS(:));
            repeatedValues = uniqueValues(histc(MS(:), uniqueValues) > 1);
            rowsRemoved = ismember(MS(:,1), repeatedValues) | ismember(MS(:,2), repeatedValues);
            MS(rowsRemoved, :) = [];
            obj.msRelation = MS;
            obj.commonNodesEdges = repeatedValues;
        end

    end
    
    methods (Access = private, Static)
        
        function closestNodeB = obtainClosestNode(nodeA,nodeB,pos)
            xA = pos(nodeA,1);
            xB = pos(nodeB,1);
            distAB = abs(xA - xB);
            [distMin,isClosest] = min(distAB);
            if distMin > 1e-4
                error('non slave node')
            end
            closestNodeB = nodeB(isClosest);
        end

        function closestNodeB = obtainClosestNode3D(nodeA, nodeB, pos1,pos2)
            xA = pos1(nodeA);
            yA = pos2(nodeA);
            
            xB = pos1(nodeB);
            yB = pos2(nodeB);
            
            distAB = sqrt((xA - xB).^2 + (yA - yB).^2);
            [distMin, isClosest] = min(distAB);
            
            if distMin > 1e-4
                error('non slave node')
            end
            closestNodeB = nodeB(isClosest);
        end
        
    end
    
end