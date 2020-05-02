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
            xCoords = obj.computeXcoords();            
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
            xIsoFull = obj.computeXisoFull();            
            xIsoCut  = obj.computeXisoSubCut();
            nDim  = size(xIsoCut,1);
            nNode = size(xIsoCut,2);
            xCoords = zeros(nDim,nNode,obj.nElem);
            xCoords(:,:,obj.iFull) = xIsoFull;
            xCoords(:,:,obj.iCut)  = xIsoCut;
        end
        
        function xIsoFull = computeXisoFull(obj)
            s.coord  = obj.localMesh.coord;
            s.connec = obj.localMesh.connec(obj.localSubCells(obj.iFull),:);
            m = Mesh().create(s);               
            xNodalAllIso = m.coordElem;
            xIsoFull = xNodalAllIso(:,:,:);                        
        end
        
        function xIsoQuadCut = computeXisoSubCut(obj)
            xIsoTri  = obj.xIsoSubCut;
            s.coord  = obj.localMesh.coord;
            s.connec = obj.localMesh.connec(obj.localSubCells(obj.iCut),:);
            m = Mesh().create(s);  
            xIsoQuadCut = m.computeXgauss(xIsoTri);            
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