classdef TethaedraSubMeshConnecComputer < handle
    
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
        
        function obj = TethaedraSubMeshConnecComputer(cParams)
            obj.init(cParams)            
        end
        
        function compute(obj,nodes,icase)
            obj.nElemInCase = size(nodes,1);
            obj.cellNodes   = nodes;
            nodes = zeros(obj.nSubCellsByElem,obj.nSubCellNodes,obj.nElemInCase);
            isActive = true(obj.nElemInCase,1);
            [nodesT,nodesP] = obj.computeNodesTAndNodesP(icase);
            for inode = 1:obj.nSubCellNodes
                nodes(1,inode,isActive) = nodesT(inode);
                nodes(2,inode,isActive) = nodesP(1,inode);
                nodes(3,inode,isActive) = nodesP(2,inode);
                nodes(4,inode,isActive) = nodesP(3,inode);
            end
            obj.connecCase = nodes;
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
            obj.nSubCellsByElem = 4;
        end
        
        function [nodesT,nodesP] = computeNodesTAndNodesP(obj,icase)
            icase
            switch icase
                case 4
                    nodesT = [5 6 7 4];
                    nodesP = [1 5 6 7;1 2 7 6;1 2 3 7];                    
                case 3
                    nodesT = [5 7 6 3];
                    nodesP = [5 6 7 1 ;2 7 6 1;2 4 7 1 ];                                        
                    %nodesP = [1 5 7 6;1 4 6 7;1 4 2 6];                     
                    %nodesP  = [5 6 7 2;1 5 7 2;4 1 7 2];
                case 2
                    nodesT = [5 6 7 2];
                    nodesP = [1 5 6 7;1 4 7 6;6 1 4 3];                      
                case 1
                    nodesT = [6 7 1 5];
                    nodesP = [5 6 7 2;4 2 7 6;2 3 6 4];                    
            end
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