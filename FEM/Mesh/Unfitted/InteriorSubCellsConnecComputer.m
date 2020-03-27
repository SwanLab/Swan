classdef InteriorSubCellsConnecComputer < handle
    
    properties (GetAccess = public, SetAccess = private)
        connec
    end
    
    properties (Access = private)
        isEdgeCutInElem
        levelSet
        allNodesInElem
        nSubCellsByElem
        
        subCellCases                
        isSubCellInterior        
        allSubCellsConnec
        
        allSubCellsConnecParams
    end
    
    methods (Access = public)
        
        function obj = InteriorSubCellsConnecComputer(cParams)
            obj.init(cParams);
            obj.compute();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.allNodesInElem    = cParams.allNodesInElem;
            obj.isEdgeCutInElem   = cParams.isEdgeCutInElem;
            obj.levelSet          = cParams.levelSet;
            obj.nSubCellsByElem   = 3;
            obj.allSubCellsConnecParams = cParams.allSubCellsConnecParams;
        end
        
        function compute(obj)
            obj.computeSubCellCases();
            obj.computeAllSubCellsConnec();
            obj.computeIsSubCellsInterior();
            obj.computeConnec();
        end
        
        function computeSubCellCases(obj)            
           s.isEdgeCutInElem = obj.isEdgeCutInElem;
           subCells = SubCellsCasesComputer(s);
           subCells.compute();
           obj.subCellCases = subCells.subCellCases;
        end
        
        function computeAllSubCellsConnec(obj)
            s = obj.allSubCellsConnecParams;
            s.allNodesInElem            = obj.allNodesInElem;
            s.subCellCases              = obj.subCellCases;                     
            a = AllSubCellsConnecComputer(s);
            a.compute();            
            obj.allSubCellsConnec = a.allSubCellsConnec;            
        end
           
        function computeIsSubCellsInterior(obj)
            s.subCellCases    = obj.subCellCases;
            s.levelSet        = obj.levelSet;
            s.allNodesInElem  = obj.allNodesInElem;
            s.nSubCellsByElem = obj.nSubCellsByElem; 
            subCellInt = IsSubCellInteriorComputer(s);
            subCellInt.compute();
            obj.isSubCellInterior = subCellInt.isSubCellInterior;
        end
        
        function computeConnec(obj)
            allConnec  = obj.allSubCellsConnec;
            isInterior = obj.isSubCellInterior(:);
            obj.connec = allConnec(isInterior,:);           
        end
                
    end
    
end
