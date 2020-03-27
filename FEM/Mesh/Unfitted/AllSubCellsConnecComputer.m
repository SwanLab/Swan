classdef AllSubCellsConnecComputer < handle
    
    properties (GetAccess = public, SetAccess = private)
        allSubCellsConnec                
    end
    
    properties (Access = private)
        allNodesInElem        
        subCellCases                
        nElem
        nCases               
        cellMesher 
        
        subMeshConnecParams
    end
    
    methods (Access = public)
        
        function obj = AllSubCellsConnecComputer(cParams)
            obj.init(cParams)
        end
         
        function compute(obj)                      
            nodesInSubCells = obj.initNodesInSubCells();
            for icase = 1:obj.nCases
                subCells  = obj.subCellCases(:,icase);               
                nodesCase = obj.allNodesInElem(subCells,:);                
                connecSubCell = obj.cellMesher.compute(nodesCase,icase);
                nodesInSubCells(:,:,subCells) = connecSubCell;                
            end
            obj.allSubCellsConnec = obj.computeAllSubCells(nodesInSubCells);            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)            
            obj.allNodesInElem = cParams.allNodesInElem;
            obj.subCellCases   = cParams.subCellCases;               
            obj.nElem  = size(obj.subCellCases,1);
            obj.nCases = size(obj.subCellCases,2);  
            obj.subMeshConnecParams = cParams.subMeshConnecParams;
            obj.createSubCellMesher();                        
        end
        
        function createSubCellMesher(obj)
            s = obj.subMeshConnecParams;
            obj.cellMesher = TriangleSubMeshConnecComputer(s);
        end
        
        function nodes = initNodesInSubCells(obj)
            nSubCellsByElem = obj.cellMesher.nSubCellsByElem;
            nSubCellNodes   = obj.cellMesher.nSubCellNodes;
            nodes = zeros(nSubCellsByElem,nSubCellNodes,obj.nElem);            
        end
    
        function nodes = computeAllSubCells(obj,nodes)
            nSubCellsByElem = obj.cellMesher.nSubCellsByElem;            
            nSubCellNodes   = obj.cellMesher.nSubCellNodes;            
            nSubCells       = obj.nElem*nSubCellsByElem;            
            nodes = permute(nodes,[1 3 2]);
            nodes = reshape(nodes,nSubCells,nSubCellNodes);
        end
        
    end
    
end