classdef DiagonalCoordComputer < TotalCoordinatesCalculator
    
    properties (Access = public)
        c
        theta
        nodes
        vertCoord
        boundCoord
        div
        totalCoord
    end
    
    properties (Access = private)
        intNode
        diagDiv
        div_aux
        O
    end
    
    methods (Access = public)
        
        function obj = DiagonalCoordComputer(cParams)
            obj.init(cParams);
            obj.initBoundary();
            obj.placeCentralNode();
            obj.computeTotalCoord();
        end
        
    end
    
    methods (Access = private)
        
        function placeCentralNode(obj)
            obj.intNode = obj.nodes.bound+1;
            vA = obj.vertCoord(4,:)-obj.vertCoord(1,:);
            pA = obj.vertCoord(1,:);
            centerVec = vA/2;
            obj.O = pA+centerVec;
            obj.totalCoord(obj.intNode,:) = obj.O;
        end
        
        function computeTotalCoord(obj)
            obj.diagDiv = max(obj.div);
            obj.div_aux = obj.div-1;
            obj.intNode = obj.intNode+1;
            for iDiv = 1:obj.diagDiv-1
                newVert = obj.computeNewVertices(iDiv);
                obj.obtainMasterCoord(newVert);
                obj.obtainSlaveCoord(newVert);
                obj.div_aux = obj.div_aux-1;
            end
        end
        
        function newVert = computeNewVertices(obj,iDiv)
            nsides = obj.nodes.vert;
            for iDiag = 1:nsides
                diagA = obj.O-obj.vertCoord(iDiag,:);
                vecDiv = iDiv*diagA/obj.diagDiv;
                pos = obj.vertCoord(iDiag,:)+vecDiv;
                obj.totalCoord(obj.intNode,:) = pos;
                obj.intNode = obj.intNode+1;
            end
            newVert = obj.totalCoord(obj.intNode-nsides:obj.intNode-1,:);
        end
        
        function obtainMasterCoord(obj,newVert)
            nsides = obj.nodes.vert;
            for iMaster = 1:nsides/2
                vertA = newVert(iMaster,:);
                vertB = newVert(iMaster+1,:);
                if (obj.c(1)~=obj.c(2)) || (obj.c(2)~=obj.c(3))
                    bool = 0;
                    while bool == 0
                        if norm((vertB-vertA)/obj.div_aux(iMaster)) > obj.c(iMaster)/obj.div(iMaster)
                            obj.div_aux(iMaster) = obj.div_aux(iMaster)+1;
                        else
                            bool = 1;
                        end
                    end
                end
                for intDiv = 1:obj.div_aux(iMaster)-1
                    sideVec = intDiv*(vertB-vertA)/obj.div_aux(iMaster);
                    sidePos = vertA+sideVec;
                    obj.totalCoord(obj.intNode,:) = sidePos;
                    obj.intNode = obj.intNode+1;
                end
            end
        end
        
        function obtainSlaveCoord(obj,newVert)
            nsides = obj.nodes.vert;
            iMaster = nsides/2;
            for iSlave = 1:nsides/2
                if iSlave == nsides/2
                    vertA = newVert(end,:);
                    vertB = newVert(1,:);
                else
                    vertA = newVert(iMaster+iSlave,:);
                    vertB = newVert(iMaster+iSlave+1,:);
                end
                for intDiv = 1:obj.div_aux(iSlave)-1
                    sideVec = intDiv*(vertB-vertA)/obj.div_aux(iSlave);
                    sidePos = vertA+sideVec;
                    obj.totalCoord(obj.intNode,:) = sidePos;
                    obj.intNode = obj.intNode+1;
                end
            end
        end
        
    end
     
end