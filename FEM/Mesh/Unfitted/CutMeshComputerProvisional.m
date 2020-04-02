classdef CutMeshComputerProvisional < handle
    
    properties (Access = public)
        connec
        coord
        xCoordsIso
    end
    
    properties (Access = private)
        cutEdgesComputer  
        cutPointsInElemComputer            
    end
    
    properties (Access = private)        
        cutEdgesParams
        cutCoordParams
        interiorSubCellsParams
        cutEdgesComputerParams
    end
    
    methods (Access = public)
        
        function  obj = CutMeshComputerProvisional(cParams)
            obj.init(cParams);
            obj.compute()
        end
        
    end    
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.cutEdgesParams         = cParams.cutEdgesParams;
            obj.cutCoordParams         = cParams.cutCoordParams;
            obj.cutEdgesComputerParams = cParams.cutEdgesComputerParams;
            obj.interiorSubCellsParams = cParams.interiorSubCellsParams;
        end
        
        function compute(obj)
            obj.computeCutEdges();
            obj.computeCoordinates();  
            obj.computeCutPointsInElemComputer();
            obj.computeConnec();
        end
        
        function computeCutEdges(obj)
            s = obj.cutEdgesParams;
            c = CutEdgesComputer(s);
            c.compute();
            obj.cutEdgesComputer = c;    
        end
        
        function computeCoordinates(obj)   
            s = obj.cutCoordParams;
            s.xCutEdgePoint    = obj.cutEdgesComputer.xCutEdgePoint;
            s.isEdgeCut        = obj.cutEdgesComputer.isEdgeCut;
            coordComputer = CutCoordinatesComputer(s);
            coordComputer.compute();
            obj.coord = coordComputer.coord;
        end        
        
        function computeCutPointsInElemComputer(obj)
            s = obj.cutEdgesComputerParams;
            s.isEdgeCut        = obj.cutEdgesComputer.isEdgeCut;
            s.allNodesInElemCoordParams.xCutEdgePoint = obj.cutEdgesComputer.xCutEdgePoint;
            c = CutPointsInElemComputer(s);
            c.compute();
            obj.cutPointsInElemComputer = c;
        end
         
        function computeConnec(obj)
            c = obj.cutPointsInElemComputer;
            sS.bestSubCellCaseSelector.coord = obj.coord;
            sA.subMeshConnecParams = sS;
            sA.xAllNodesInElem = c.xAllNodesInElem;
            sA.allNodesInElem  = c.allNodesInElem;
            sC.isEdgeCutInElem = c.isEdgeCutInElem;
            s = obj.interiorSubCellsParams;          
            sI = s.isSubCellInteriorParams;
            sI.allNodesInElem = c.allNodesInElem;
            s.allSubCellsConnecParams = sA;
            s.subCellsCasesParams = sC; 
            s.isSubCellInteriorParams = sI;
            subCell = InteriorSubCellsConnecComputer(s);
            obj.connec = subCell.connec;
            obj.xCoordsIso = subCell.xCoordsIso;
        end               
        
    end
    
end