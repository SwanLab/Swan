classdef AllSubCellsConnecComputer < handle
    
    properties (GetAccess = public, SetAccess = private)
        allSubCellsConnec 
        xNodesInSubCells
        yNodesInSubCells
    end
    
    properties (Access = private)
        allNodesInElem        
        subCellCases       
        xAllNodesInElem
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
            xNodesInSubCells = obj.initNodesInSubCells();
            yNodesInSubCells = obj.initNodesInSubCells();
            for icase = 1:obj.nCases
                subCells  = obj.subCellCases(:,icase);                                           

                nodesCase = obj.allNodesInElem(subCells,:);
                obj.cellMesher.compute(nodesCase,icase);
                
                nodesInSubCells(:,:,subCells) = obj.cellMesher.partition(nodesCase);
                
                %nodesInSubCells(:,:,subCells) = obj.computeCase(subCells,icase,obj.allNodesInElem);              
                xNodesInSubCells(:,:,subCells) = obj.cellMesher.partition(obj.xAllNodesInElem(subCells,:,1)); 
                yNodesInSubCells(:,:,subCells) = obj.cellMesher.partition(obj.xAllNodesInElem(subCells,:,2));              


            end
            obj.allSubCellsConnec = obj.computeAllSubCells(nodesInSubCells);      
            obj.xNodesInSubCells = obj.computeAllSubCells(xNodesInSubCells);
            obj.yNodesInSubCells = obj.computeAllSubCells(yNodesInSubCells);            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)            
            obj.allNodesInElem = cParams.allNodesInElem;
            obj.subCellCases   = cParams.subCellCases;            
            obj.xAllNodesInElem = cParams.xAllNodesInElem;
            obj.nElem  = size(obj.subCellCases,1);
            obj.nCases = size(obj.subCellCases,2);  
            obj.subMeshConnecParams = cParams.subMeshConnecParams;
            obj.createSubCellMesher();                        
        end
        
        function connecSubCell = computeCase(obj,subCells,icase,allNodesInElem)
            nodesCase = allNodesInElem(subCells,:);
            connecSubCell = obj.cellMesher.compute(nodesCase,icase);            
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