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
        vertices
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
%                     obj.get_vertices(coord);
%                     obj.get_MasterSlave(coord)
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
            [masterY,slaveY] = obj.computeMasterSlaveNode(obj.nodesInXmin,obj.nodesInXmax,obj.y);
            [masterZ,slaveZ] = obj.computeMasterSlaveNode(obj.nodesInXmin,obj.nodesInXmax,obj.z);
            master = [masterY; masterZ];
            slave  = [slaveY; slaveZ];
        end
        
        function [master,slave] = computeMasterSlaveNodesInFaceY3D(obj)
            [master,slave] = obj.computeMasterSlaveNode(obj.nodesInYmin,obj.nodesInYmax,obj.x);
        end
        
        function [master,slave] = computeMasterSlaveNodesInFaceZ3D(obj)
            [master,slave] = obj.computeMasterSlaveNode(obj.nodesInZmin,obj.nodesInZmax,obj.x);
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

        function init3(obj,coord)
            obj.x = coord(:,1);
            obj.y = coord(:,2);
            obj.z = coord(:,3);
            obj.allNodes(:,1) = 1:size(obj.x,1);
        end

        function get_vertices(obj,coord)
            v = [];
            for i = 1:length(coord) %   find ????
                nodeGID = coord(i,:);
                if (isequal(nodeGID, [0,0,0]) || isequal(nodeGID, [0,0,1]) || ...
                        isequal(nodeGID, [0,1,0]) || isequal(nodeGID, [1,0,0]) || ...
                        isequal(nodeGID, [0,1,1]) || isequal(nodeGID, [1,0,1]) || ...
                        isequal(nodeGID, [1,1,0]) || isequal(nodeGID, [1,1,1]))
                    v = [v; nodeGID];
                end
            end
            obj.vertices = v;
        end 

        function get_MasterSlave(obj,gidcoord)
            MS = [ ];
            v = obj.vertices;
            for i = 1:length(gidcoord)
                if ismember(gidcoord(i,:),v,'rows')
                else % if node is not a vertice
                    if any(gidcoord==0)
                        xCoord = gidcoord(i,1);
                        yCoord = gidcoord(i,2);
                        zCoord = gidcoord(i,3);
                        if xCoord==0
                            for j = 1:length(gidcoord)
                                nodeCompared = gidcoord(j,1:3);
                                nodeSuposed = [xCoord+1,yCoord,zCoord];
                                if (isequal(nodeCompared,nodeSuposed))
                                    MS = [MS; [i j]];
                                end
                            end
                        end
                        if yCoord==0
                            for j = 1:length(gidcoord)
                                nodeCompared = gidcoord(j,1:3);
                                nodeSuposed = [xCoord,yCoord+1,zCoord];
                                if (isequal(nodeCompared,nodeSuposed))
                                    MS = [MS; [i j]];
                                end
                            end
                        end
                        if zCoord==0
                            for j = 1:length(gidcoord)
                                nodeCompared = gidcoord(j,1:3);
                                nodeSuposed = [xCoord,yCoord,zCoord+1];
                                if (isequal(nodeCompared,nodeSuposed))
                                    MS = [MS; [i j]];
                                end
                            end
                        end
                    end
                end
            end
            obj.msRelation = MS;
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