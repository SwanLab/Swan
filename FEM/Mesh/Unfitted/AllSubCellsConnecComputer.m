classdef AllSubCellsConnecComputer < handle
    
    properties (GetAccess = public, SetAccess = private)
        allSubCellsConnec
        xNodesInSubCells
    end
    
    properties (Access = private)
        nodesInSubCells
        xNodesInSubCellsByElem
        cellMesher
        
    end
    
    properties (Access = private)
        allNodesInElem
        xAllNodesInElem        
        subCellCases
        nElem
        nCases
        
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
                obj.updateCellMesherPartitioner(subCells,icase);
                obj.computePartitionConnecSubCell(subCells);
                obj.computePartitionCoordSubCell(subCells);
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
            obj.cellMesher = TriangleSubMeshConnecComputer(s);
        end
        
        function permuteXallNodes(obj)
            xNodes = obj.xAllNodesInElem;
            xNodes = permute(xNodes,[3 2 1]);
            obj.xAllNodesInElem = xNodes;
        end
        
        function initNodesInSubCells(obj)
            nSubCellsByElem = obj.cellMesher.nSubCellsByElem;
            nSubCellNodes   = obj.cellMesher.nSubCellNodes;
            nodes = zeros(nSubCellsByElem,nSubCellNodes,obj.nElem);
            obj.nodesInSubCells = nodes;
        end
        
        function initXnodesInSubCellsByElem(obj)
            nSubCellsByElem = obj.cellMesher.nSubCellsByElem;
            nSubCellNodes   = obj.cellMesher.nSubCellNodes;
            nDim = size(obj.xAllNodesInElem,3);
            nodes = zeros(nSubCellsByElem,nSubCellNodes,obj.nElem,nDim);
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
            for idim = 1:nDim
                x = obj.xNodesInSubCellsByElem(:,:,:,idim);
                xAll = obj.concatenateSubCells(x);
                obj.xNodesInSubCells(:,:,idim) = xAll;
            end
        end
        
        function concatenateAllNodesConnec(obj)
            nodes = obj.nodesInSubCells;
            allNodes = obj.concatenateSubCells(nodes);
            obj.allSubCellsConnec = allNodes;
        end
        
        function allNodes = concatenateSubCells(obj,nodesByElem)
            nSubCellsByElem = obj.cellMesher.nSubCellsByElem;
            nSubCellNodes   = obj.cellMesher.nSubCellNodes;
            nSubCells       = obj.nElem*nSubCellsByElem;
            nodesByElem = permute(nodesByElem,[1 3 2]);
            allNodes = reshape(nodesByElem,nSubCells,nSubCellNodes);
        end
        
    end
    
end