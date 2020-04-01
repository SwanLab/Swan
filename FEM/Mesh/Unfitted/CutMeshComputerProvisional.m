classdef CutMeshComputerProvisional < handle
    
    properties (Access = public)
        connec
        coord
    end
    
    properties (Access = private)
        backgroundConnec
        levelSet
        
        cutEdgesComputer        
        allNodesInElem  
        isEdgeCutInElem
        xAllNodesInElem
        
        cutEdgesParams
        cutCoordParams
        allNodesinElemParams
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
            obj.backgroundConnec       = cParams.connec;
            obj.cutEdgesParams         = cParams.cutEdgesParams;
            obj.cutCoordParams         = cParams.cutCoordParams;
            obj.cutEdgesComputerParams = cParams.cutEdgesComputerParams;
            obj.levelSet = cParams.levelSet;
        end
        
        function compute(obj)
            obj.computeCutEdges();
            obj.computeCoordinates();  
            obj.computeCutEdgesComputer();
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
        
        function computeCutEdgesComputer(obj)
            s = obj.cutEdgesComputerParams;
            s.backgroundConnec = obj.backgroundConnec;
            s.isEdgeCut        = obj.cutEdgesComputer.isEdgeCut;
            s.xCutEdgePoint    = obj.cutEdgesComputer.xCutEdgePoint;
            c = CutPointsInElemComputer(s);
            c.compute();
            obj.allNodesInElem = c.allNodesInElem;  
            obj.isEdgeCutInElem = c.isEdgeCutInElem;
            obj.xAllNodesInElem = c.xAllNodesInElem;
        end
         
        function computeConnec(obj)
            sS.bestSubCellCaseSelector.coord = obj.coord;
            sA.subMeshConnecParams = sS;
            sA.xAllNodesInElem = obj.xAllNodesInElem;
            s.allSubCellsConnecParams = sA;
            s.allNodesInElem   = obj.allNodesInElem;
            s.isEdgeCutInElem  = obj.isEdgeCutInElem;            
            s.levelSet         = obj.levelSet;            
            subCell = InteriorSubCellsConnecComputer(s);
            obj.connec = subCell.connec;
        end               
        
    end
    
end