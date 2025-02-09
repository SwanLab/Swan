classdef AllSubCellsConnecComputer < handle
    
    properties (GetAccess = public, SetAccess = private)
        allSubCellsConnec
        xNodesInSubCells
    end
    
    properties (Access = private)
        xNodesInSubCellsByElem
        cellMesher
    end
    
    properties (Access = private)
        allNodesInElem
        xAllNodesInElem
        subCellCases
        nElem
        nCases
        nodesInSubCells
        subMeshConnecParams
    end
    
    methods (Access = public)
        
        function obj = AllSubCellsConnecComputer(cParams)
            obj.init(cParams)
        end
        
        function compute(obj)
            obj.permuteXallNodes();
            obj.initNodesInSubCells();
            obj.initXnodesInSubCellsByElem();
            for icase = 1:obj.nCases
                subCells = obj.subCellCases(:,icase);
                if sum(subCells) > 0
                    obj.updateCellMesherPartitioner(subCells,icase);
                    obj.computePartitionConnecSubCell(subCells);
                    obj.computePartitionCoordSubCell(subCells);
                end
            end
            obj.concatenateAllNodesConnec();
            obj.concatenateAllNodesCoordinates();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.allNodesInElem  = cParams.allNodesInElem;
            obj.subCellCases    = cParams.subCellCases;
            obj.xAllNodesInElem = cParams.xAllNodesInElem;
            obj.nElem  = size(obj.subCellCases,1);
            obj.nCases = size(obj.subCellCases,2);
            obj.subMeshConnecParams = cParams.subMeshConnecParams;
            obj.createSubCellMesher();
        end
        
        function createSubCellMesher(obj)
            s = obj.subMeshConnecParams;
            switch size(obj.xAllNodesInElem,2)
                case 5
                    obj.cellMesher = TriangleSubMeshConnecComputer(s);
                case {7,8}
                    obj.cellMesher = TetrahedraSubMeshConnecComputer(s);
            end
        end
        
        function permuteXallNodes(obj)
            xNodes = obj.xAllNodesInElem;
            xNodes = permute(xNodes,[3 2 1]);
            obj.xAllNodesInElem = xNodes;
        end
        
        function initNodesInSubCells(obj)
            switch mode(size(obj.allNodesInElem,2))
                case 5
                    nSubCellsByElem = 3;
                case 7
                    nSubCellsByElem = 4;
                case 8
                    nSubCellsByElem = 6;
            end
            nSubCellNodes   = obj.cellMesher.nSubCellNodes;
            nodes = zeros(nSubCellNodes,nSubCellsByElem,obj.nElem);
            obj.nodesInSubCells = nodes;
        end
        
        function initXnodesInSubCellsByElem(obj)
            switch mode(size(obj.allNodesInElem,2))
                case 5
                    nSubCellsByElem = 3;
                case 7
                    nSubCellsByElem = 4;
                case 8
                    nSubCellsByElem = 6;
            end
            nSubCellNodes   = obj.cellMesher.nSubCellNodes;
            nDim = size(obj.xAllNodesInElem,3);
            nodes = zeros(nSubCellNodes,nSubCellsByElem,obj.nElem,nDim);
            obj.xNodesInSubCellsByElem = nodes;
        end
        
        function updateCellMesherPartitioner(obj,subCells,icase)
            nodesCase = obj.allNodesInElem(subCells,:);
            obj.cellMesher.compute(nodesCase,icase);
        end
        
        function computePartitionConnecSubCell(obj,subCells)
            nodes = obj.allNodesInElem(subCells,:);
            dividedNodes = obj.cellMesher.partition(nodes);
            obj.nodesInSubCells(:,:,subCells) = dividedNodes;
        end
        
        function computePartitionCoordSubCell(obj,subCells)
            xNodes = obj.xAllNodesInElem;
            nDim = size(xNodes,3);
            for idim = 1:nDim
                xN = xNodes(subCells,:,idim);
                xP = obj.cellMesher.partition(xN);
                obj.xNodesInSubCellsByElem(:,:,subCells,idim) = xP;
            end
        end
        
        function concatenateAllNodesCoordinates(obj)
            nDim = size(obj.xAllNodesInElem,3);
            nSubCellsByElem = obj.cellMesher.nSubCellsByElem;
            nSubCellNodes   = obj.cellMesher.nSubCellNodes;
            nSubCells       = obj.nElem*nSubCellsByElem;
            obj.xNodesInSubCells = zeros(nSubCellNodes,nSubCells,nDim);
            for idim = 1:nDim
                x = obj.xNodesInSubCellsByElem(:,:,:,idim);
                xAll = obj.concatenateSubCells(x);
                obj.xNodesInSubCells(:,:,idim) = xAll;
            end
        end
        
        function concatenateAllNodesConnec(obj)
            nodes = obj.nodesInSubCells;
            allNodes = obj.concatenateSubCells(nodes);
            allNodes = transpose(allNodes);
            obj.allSubCellsConnec = allNodes;
        end
        
        function allNodes = concatenateSubCells(obj,nodesByElem)
            nSubCellsByElem = obj.cellMesher.nSubCellsByElem;
            nSubCellNodes   = obj.cellMesher.nSubCellNodes;
            nSubCells       = obj.nElem*nSubCellsByElem;
            allNodes = reshape(nodesByElem,nSubCellNodes,nSubCells);
        end
        
    end
    
end