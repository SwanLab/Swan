classdef CutMeshComputerProvisional < handle
    
    properties (Access = public)
        connec
        coord
    end
    
    properties (Access = private)
        backgroundConnec
        backgroundCoord
        levelSet
        
        edgesComputer
        cutEdgesComputer        
   
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
            obj.levelSet = cParams.levelSet;
        end
        
        function compute(obj)
            obj.computeEdges();
            obj.computeCutEdges();         
            obj.computeCoordinates();                                      
            obj.computeConnec();
        end
        
        function computeEdges(obj)
            s.nodesByElem = obj.backgroundConnec;
            e = EdgesConnectivitiesComputer(s);
            e.compute();
            obj.edgesComputer = e;
        end
        
        function computeCutEdges(obj)
            s.levelSet = obj.levelSet;
            s.edgesComputer = obj.edgesComputer;
            c = CutEdgesComputer(s);
            c.compute();
            obj.cutEdgesComputer = c;    
        end
        
        function computeCoordinates(obj)           
            s.nodesInCutEdges     = obj.cutEdgesComputer.nodesInCutEdges;
            s.levelSet     = obj.levelSet;
            s.backgroundCoord = obj.backgroundCoord;
            coordComputer = CutCoordinatesComputer(s);
            coordComputer.compute();
            obj.coord = coordComputer.coord;
        end        
        
        function computeConnec(obj)
            s.cutNodePerElemen = obj.computeCutNodePerElem();
            s.nCutEdgeByElem = obj.cutEdgesComputer.nCutEdgeByElem;
            s.backgroundConnec = obj.backgroundConnec;
            s.nElem = size(obj.backgroundConnec,1);
            s.elemCases = obj.cutEdgesComputer.elemCases;
            s.coord = obj.coord;
            s.levelSet = obj.levelSet;
            subCell = SubCellConnecComputer(s);
            obj.connec = subCell.connec;
        end        
        
        function cutNodePerElemen = computeCutNodePerElem(obj)
            firstCutEdgePerElem = obj.cutEdgesComputer.firstCutEdge;
            finalNode = size(obj.backgroundCoord,1);                                    
            cutNodePerElemen = firstCutEdgePerElem + finalNode;
        end        
        
        
    end
    
end