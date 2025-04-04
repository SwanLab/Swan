classdef XcoordIsoComputer < handle
    
    properties (Access = private)
        iFull
        iCut
        nElem
        localSubCells
        allSubCells
        localSubMesh
        xIsoFull
        xIsoCut
    end
    
    properties (Access = private)
        globalToLocal
        fullCells
        cutCells
        localMesh
        xIsoSubCut
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
            obj.computeXisoFull(); 
            obj.computeXisoSubCut();
            xCoords = obj.computeXcoords();
        end
        
        function x = computeXSubCut(obj)
            obj.computeNelem();
            obj.computeFullCellIndex();
            obj.computeCutCellIndex();
            obj.computeAllSubCells();
            obj.computeLocalSubCells();
            obj.computeXisoSubCut();
            x = obj.xIsoCut;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.fullCells     = cParams.fullCells;
            obj.cutCells      = cParams.cutCells;
            obj.globalToLocal = cParams.globalToLocal;
            obj.localMesh     = cParams.localMesh;
            obj.xIsoSubCut    = cParams.xIsoCutCoord;
        end
        
        function xCoords = computeXcoords(obj)
            nDim  = size(obj.xIsoCut,1);
            nNode = size(obj.xIsoCut,2);
            xCoords = zeros(nDim,nNode,obj.nElem);
            xCoords(:,:,obj.iFull) = obj.xIsoFull;
            xCoords(:,:,obj.iCut)  = obj.xIsoCut;
        end
        
        function computeXisoFull(obj)
            s.coord  = obj.localMesh.coord;
            s.connec = obj.localMesh.connec(obj.localSubCells(obj.iFull),:);
            m = Mesh.create(s);
            xNodalAllIso = m.coordElem;
            obj.xIsoFull = xNodalAllIso;
        end
        
        function computeXisoSubCut(obj)
            xIsoTri  = obj.xIsoSubCut;
            s.coord  = obj.localMesh.coord;
            s.connec = obj.localMesh.connec(obj.localSubCells(obj.iCut),:);
            m = Mesh.create(s);
            obj.xIsoCut = m.computeXgauss(xIsoTri);
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
            aSubCells = obj.allSubCells;
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