classdef XcoordIsoComputer < handle
    
    properties (Access = private)
        iFull
        iCut
        nElem
        localSubCells
        allSubCells
        localSubMesh
    end
    
    properties (Access = private)
        globalToLocal
        fullCells
        cutCells
        localMesh
        xIsoCutCoord
    end
    
    methods (Access = public)
        
        function obj = XcoordIsoComputer(cParams)
            obj.init(cParams);
        end
        
        function xCoords = compute(obj)
            obj.computeNelem();
            obj.computeFullCellIndex();
            obj.computeCutCellIndex();
            obj.computeAllSubCells();
            obj.computeLocalSubCells();
            obj.computeLocalSubMesh();            
            xCoords = obj.computeXcoords();            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.fullCells     = cParams.fullCells;
            obj.cutCells      = cParams.cutCells;
            obj.globalToLocal = cParams.globalToLocal;
            obj.localMesh     = cParams.localMesh;
            obj.xIsoCutCoord  = cParams.xIsoCutCoord;
        end
        
        function computeLocalSubMesh(obj)
            s.coord  = obj.localMesh.coord;
            s.connec = obj.localMesh.connec(obj.localSubCells,:);
            obj.localSubMesh = Mesh().create(s);                        
        end
        
        function xCoords = computeXcoords(obj)
            xIsoAll = obj.computeXisoAll();            
            xCoords = obj.localSubMesh.computeXgauss(xIsoAll);
        end                
        
        function xIsoAll = computeXisoAll(obj)   
            xIsoFull = obj.computeXisoFull();            
            xIsoCut  = obj.xIsoCutCoord;
            nDim  = size(xIsoCut,1);
            nNode = size(xIsoCut,2);
            xIsoAll = zeros(nDim,nNode,obj.nElem);
            xIsoAll(:,:,obj.iFull) = xIsoFull;
            xIsoAll(:,:,obj.iCut)  = xIsoCut;
        end
        
        function xIsoFull = computeXisoFull(obj)
            xNodalAllIso = obj.localSubMesh.coordElem;            
            xIsoFull = xNodalAllIso(:,:,obj.iFull);                        
        end
        
        function computeNelem(obj)
            nFull = size(obj.fullCells,1);
            nCut  = size(obj.cutCells,1);
            obj.nElem = nFull + nCut;
        end
        
        function computeFullCellIndex(obj)
            nFull = size(obj.fullCells,1);
            iF = false(obj.nElem,1);
            iF(1:nFull,1) = true;            
            obj.iFull = iF;           
        end
        
        function computeCutCellIndex(obj)
            nFull = size(obj.fullCells,1);
            nCut  = size(obj.cutCells,1);
            iC  = false(obj.nElem,1);            
            iC(nFull + (1:nCut),1) = true;            
            obj.iCut  = iC;            
        end

        function computeLocalSubCells(obj)
            aSubCells = obj.allSubCells();
            obj.localSubCells = obj.globalToLocal(aSubCells);
        end
        
        function computeAllSubCells(obj)
            aC = zeros(obj.nElem,1);
            aC(obj.iFull,1) = obj.fullCells;
            aC(obj.iCut,1)  = obj.cutCells;
            obj.allSubCells = aC;            
        end
        
    end
    
    
end