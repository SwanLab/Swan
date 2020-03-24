classdef CutMeshComputerProvisional < handle
    
    properties (Access = public)
        connec
        coord
    end
    
    properties (Access = private)
        backgroundConnec
        backgroundCoord
        levelSet
        
        cutEdgesComputer        
        allNodesInElem    
        cutEdgesParams
    end
    
    methods (Access = public)
        
        function  obj = CutMeshComputerProvisional(cParams)
            obj.init(cParams);
            obj.compute()
        end
        
    end    
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.backgroundConnec   = cParams.connec;
            obj.backgroundCoord    = cParams.coord;
            obj.cutEdgesParams     = cParams.cutEdgesParams;
            obj.levelSet = cParams.levelSet;
        end
        
        function compute(obj)
            obj.computeCutEdges();    
            obj.computeAllNodesInElem();
            obj.computeCoordinates();                                      
            obj.computeConnec();
        end
        
        function computeCutEdges(obj)
            s = obj.cutEdgesParams;
            s.levelSet = obj.levelSet;
            c = CutEdgesComputer(s);
            c.compute();
            obj.cutEdgesComputer = c;    
        end
        
        function computeAllNodesInElem(obj)
            s.finalNodeNumber  = size(obj.backgroundCoord,1);
            s.firstCutEdge     = obj.cutEdgesComputer.firstCutEdge;
            s.backgroundConnec = obj.backgroundConnec;
            aComputer = AllNodesInElemComputer(s);
            aComputer.compute();
            obj.allNodesInElem = aComputer.allNodesInElem;
        end
        
        function computeCoordinates(obj)           
            s.nodesInCutEdges  = obj.cutEdgesComputer.nodesInCutEdges;
            s.levelSet         = obj.levelSet;
            s.backgroundCoord  = obj.backgroundCoord;
            coordComputer = CutCoordinatesComputer(s);
            coordComputer.compute();
            obj.coord = coordComputer.coord;
        end        
        
        function computeConnec(obj)
            s.firstCutEdge     = obj.cutEdgesComputer.firstCutEdge;
            s.allNodesInElem   = obj.allNodesInElem;
            s.nElem            = size(obj.backgroundConnec,1);
            s.isEdgeCutInElem  = obj.cutEdgesComputer.isEdgeCutInElem;
            s.coord            = obj.coord;
            s.levelSet         = obj.levelSet;            
            subCell = InteriorSubCellsConnecComputer(s);
            obj.connec = subCell.connec;
        end               
        
    end
    
end