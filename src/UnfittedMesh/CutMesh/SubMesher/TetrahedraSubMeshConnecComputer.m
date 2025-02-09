classdef TetrahedraSubMeshConnecComputer < handle
    
    properties (Access = public)
        nSubCellNodes
        nSubCellsByElem
    end
    
    properties (Access = private)
        connecCase
    end
    
    properties (Access = private)
        nElemInCase
        cellNodes
    end
    
    methods (Access = public)
        
        function obj = TetrahedraSubMeshConnecComputer(cParams)
            obj.init(cParams)
        end
        
        function compute(obj,nodes,icase)
            
            switch mode(size(nodes,2))
                case 7
                    obj.nSubCellsByElem = 4;
                case 8
                    obj.nSubCellsByElem = 6;
            end
            
            obj.nElemInCase = size(nodes,1);
            obj.cellNodes   = nodes;
            nodesC = zeros(obj.nSubCellsByElem,obj.nSubCellNodes,obj.nElemInCase);
            isActive = true(obj.nElemInCase,1);
            switch mode(size(nodes,2))
                case 7
                     [nodesT,nodesP] = obj.computeNodesTAndNodesP(icase);
                    for inode = 1:obj.nSubCellNodes
                        nodesC(1,inode,isActive) = nodesT(inode);
                        nodesC(2,inode,isActive) = nodesP(1,inode);
                        nodesC(3,inode,isActive) = nodesP(2,inode);
                        nodesC(4,inode,isActive) = nodesP(3,inode);
                    end
                case 8
                    nodesCC = obj.computeNodesFourCutNodes(icase);
                    for inode = 1:obj.nSubCellNodes
                        nodesC(1,inode,isActive) = nodesCC(1,inode);
                        nodesC(2,inode,isActive) = nodesCC(2,inode);
                        nodesC(3,inode,isActive) = nodesCC(3,inode);
                        nodesC(4,inode,isActive) = nodesCC(4,inode);
                        nodesC(5,inode,isActive) = nodesCC(5,inode);
                        nodesC(6,inode,isActive) = nodesCC(6,inode);
                    end 
                    
                    
            end
            obj.connecCase = nodesC;
        end
        
        function nodePartition = partition(obj,nodes)
            connec = permute(obj.connecCase,[3 1 2]);
            nodePartition = obj.initNodePartition();
            nNodesSubCell = size(nodes,2);
            for isubcell = 1:obj.nSubCellsByElem
                for inode = 1:obj.nSubCellNodes
                    node = connec(:,isubcell,inode);
                    ind  = obj.computeAbsoluteMatrixIndex(node,nNodesSubCell);
                    nodePartition(inode,isubcell,:) = nodes(ind);
                end
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.nSubCellNodes = 4;
        end
        
        function nodes = computeNodesFourCutNodes(obj,icase)
            
            switch icase
                case 3
                    nodes1 = obj.prismaTriangulation([1 5 6],[3 7 8]);
                    nodes2 = obj.prismaTriangulation([2 5 7],[4 6 8]);
                    nodes = [nodes1;nodes2];
                case 2
                    nodes1 = obj.prismaTriangulation([2 7 8],[1 5 6]);
                    nodes2 = obj.prismaTriangulation([5 3 7],[6 4 8]);
                    nodes = [nodes1;nodes2];
                case 1
                    nodes1 = obj.prismaTriangulation([4 7 8],[1 5 6]);
                    nodes2 = obj.prismaTriangulation([5 2 7],[6 3 8]);
                    nodes = [nodes1;nodes2];
            end
            
            
        end
        
        function [nodesT,nodesP] = computeNodesTAndNodesP(obj,icase)
            switch icase
                case 4
                    nodesT = [5 6 7 4];
                    nodesP = obj.prismaTriangulation([5 6 7],[1 2 3]);
                    %nodesP = [1 5 6 7;1 2 7 6;1 2 3 7];
                case 3
                    nodesP = [5 6 7 1 ;2 7 6 1;2 4 7 1 ];
                    nodesP = obj.prismaTriangulation([5 7 6],[1 4 2]);
                    nodesT = [5 7 6 3];
                    
                    
                    %nodesP = [1 5 7 6;1 4 6 7;1 4 2 6];
                    %nodesP  = [5 6 7 2;1 5 7 2;4 1 7 2];
                case 2
                    nodesT = [5 6 7 2];
%                    nodesP = [1 5 6 7;1 4 7 6;6 1 4 3];
                    nodesP = obj.prismaTriangulation([5 6 7],[1 3 4]);
                    
                case 1
                    nodesT = [6 7 1 5];
%                    nodesP = [5 6 7 2;4 2 7 6;2 3 6 4];
                    nodesP = obj.prismaTriangulation([5 7 6],[2 4 3]);
                    
            end
            
            
        end
        
       function nodes = prismaTriangulation2(obj,nodesA,nodesB)
             nodes =  [nodesA(1) nodesA(2) nodesA(3) nodesB(1);
                       nodesA(1) nodesA(2) nodesA(3) nodesB(2);
                       nodesA(1) nodesA(2) nodesA(3) nodesB(3)];
       end
        
        function nodes = prismaTriangulation(obj,nodesA,nodesB)
             nodes =  [nodesB(1) nodesA(1) nodesA(2) nodesA(3);
                       nodesB(1) nodesB(2) nodesA(3) nodesA(2);
                       nodesB(1) nodesB(2) nodesB(3) nodesA(3)];
        end
        
        
        function c = initNodePartition(obj)
            c = zeros(obj.nSubCellNodes,obj.nSubCellsByElem,obj.nElemInCase);
        end
        
        
        function ind = computeAbsoluteMatrixIndex(obj,colums,nNodesSubCell)
            nElem = obj.nElemInCase;
            rows = (1:nElem)';
            ind = sub2ind([nElem,nNodesSubCell],rows,colums);
        end
        
    end
    
end