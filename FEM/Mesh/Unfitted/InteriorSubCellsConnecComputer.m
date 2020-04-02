classdef InteriorSubCellsConnecComputer < handle
    
    properties (GetAccess = public, SetAccess = private)
        connec
        xCoordsIso
    end
    
    properties (Access = private)
        subCellCases                
        isSubCellInterior        
        allSubCellsConnec        
        xNodesInSubCells        
    end
    
    properties (Access = private)
        nSubCellsByElem
        allSubCellsConnecParams
        isSubCellInteriorParams
        subCellsCasesParams
    end
    
    methods (Access = public)
        
        function obj = InteriorSubCellsConnecComputer(cParams)
            obj.init(cParams);
            obj.compute();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)        
            obj.allSubCellsConnecParams = cParams.allSubCellsConnecParams;
            obj.isSubCellInteriorParams = cParams.isSubCellInteriorParams;
            obj.subCellsCasesParams     = cParams.subCellsCasesParams;
            obj.nSubCellsByElem   = 3;            
        end
        
        function compute(obj)
            obj.computeSubCellCases();
            obj.computeAllSubCellsConnec();
            obj.computeIsSubCellsInterior();
            obj.computeConnec();
            obj.computeXcoordsIso();
        end
        
        function computeSubCellCases(obj)            
           s = obj.subCellsCasesParams;
           subCells = SubCellsCasesComputer(s);
           subCells.compute();
           obj.subCellCases = subCells.subCellCases;
        end
        
        function computeAllSubCellsConnec(obj)
            s = obj.allSubCellsConnecParams;
            s.subCellCases              = obj.subCellCases;                     
            a = AllSubCellsConnecComputer(s);
            a.compute();            
            obj.allSubCellsConnec = a.allSubCellsConnec;  
            obj.xNodesInSubCells  = a.xNodesInSubCells;
        end
           
        function computeIsSubCellsInterior(obj)
            s = obj.isSubCellInteriorParams;
            s.subCellCases    = obj.subCellCases;
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
                
        function computeXcoordsIso(obj)
            isInterior = obj.isSubCellInterior(:);            
            x = obj.xNodesInSubCells;
            nDim = size(x,3);
            for idim = 1:nDim
                obj.xCoordsIso(:,idim,:) = x(isInterior,:,idim);
            end
        end
                
    end
    
end
