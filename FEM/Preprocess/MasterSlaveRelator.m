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
%                     obj.get_cornerNodes(coord);
%                     obj.get_MasterSlave(coord);
%                     obj.get_MasterSlave2(coord); 
            end
        end
        
        function r = getRelation(obj)
            r = obj.msRelation;
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

%         function get_cornerNodes(obj,coord)
%             v = [];
%             for i = 1:length(coord) %   find ????
%                 nodeGID = coord(i,:);
%                 if (isequal(nodeGID, [0,0,0]) || isequal(nodeGID, [0,0,1]) || ...
%                         isequal(nodeGID, [0,1,0]) || isequal(nodeGID, [1,0,0]) || ...
%                         isequal(nodeGID, [0,1,1]) || isequal(nodeGID, [1,0,1]) || ...
%                         isequal(nodeGID, [1,1,0]) || isequal(nodeGID, [1,1,1]))
%                     v = [v; nodeGID];
%                 end
%             end
%             obj.cornerNodes = v;
%         end 
% 
%         function get_MasterSlave(obj,gidcoord)
%             MS = [ ];
%             v = obj.cornerNodes;
%             for i = 1:length(gidcoord)
%                 if ismember(gidcoord(i,:),v,'rows')
%                 else % if node is not a vertice
%                     if any(gidcoord==0)
%                         xCoord = gidcoord(i,1);
%                         yCoord = gidcoord(i,2);
%                         zCoord = gidcoord(i,3);
%                         if xCoord==0
%                             for j = 1:length(gidcoord)
%                                 nodeCompared = gidcoord(j,1:3);
%                                 nodeSuposed = [xCoord+1,yCoord,zCoord];
%                                 if (isequal(nodeCompared,nodeSuposed))
%                                     MS = [MS; [i j]];
%                                 end
%                             end
%                         end
%                         if yCoord==0
%                             for j = 1:length(gidcoord)
%                                 nodeCompared = gidcoord(j,1:3);
%                                 nodeSuposed = [xCoord,yCoord+1,zCoord];
%                                 if (isequal(nodeCompared,nodeSuposed))
%                                     MS = [MS; [i j]];
%                                 end
%                             end
%                         end
%                         if zCoord==0
%                             for j = 1:length(gidcoord)
%                                 nodeCompared = gidcoord(j,1:3);
%                                 nodeSuposed = [xCoord,yCoord,zCoord+1];
%                                 if (isequal(nodeCompared,nodeSuposed))
%                                     MS = [MS; [i j]];
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%             obj.msRelation = MS;
%         end
% 
%         function get_MasterSlave2(obj,gidcoord)
%             MS = [ ];
%             v = obj.cornerNodes;
%             n = length(gidcoord);
%             [~, idx] = ismember(gidcoord,v,'rows');
%             noncornerNodes = find(~idx);
%             
%             if any(gidcoord(noncornerNodes,:) == 0)
%                 for i = noncornerNodes'
%                     xCoord = gidcoord(i,1);
%                     yCoord = gidcoord(i,2);
%                     zCoord = gidcoord(i,3);
%                     for j = 1:n
%                         if j == i
%                             continue
%                         end
%                         nodeCompared = gidcoord(j,1:3);
%                         if xCoord == 0 && isequal(nodeCompared, [xCoord+1,yCoord,zCoord])
%                             MS = [MS; [i j]];
%                         elseif yCoord == 0 && isequal(nodeCompared, [xCoord,yCoord+1,zCoord])
%                             MS = [MS; [i j]];
%                         elseif zCoord == 0 && isequal(nodeCompared, [xCoord,yCoord,zCoord+1])
%                             MS = [MS; [i j]];
%                         end
%                     end
%                 end
%             end
%             
%             obj.msRelation = MS;
%         end

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

        function closestNodeB = obtainClosestNode3D(nodeA, nodeB, pos1,pos2)
            xA = pos1(nodeA);
            yA = pos2(nodeA);
            
            xB = pos1(nodeB);
            yB = pos2(nodeB);
            
            distAB = sqrt((xA - xB).^2 + (yA - yB).^2);
            [distMin, isClosest] = min(distAB);
            
            if distMin > 1e-10
                error('non slave node')
            end
            closestNodeB = nodeB(isClosest);
        end
        
    end
    
end