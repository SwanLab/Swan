classdef CutMesh < handle
    
    properties (Access = private)
        cutMeshOfSubCellLocal
        cutMeshOfSubCellGlobal
        cellContainingSubcell        
        mesh        
        
        backgroundCutMesh
        subcellIsoCoords
    end
    
    properties (Access = private)
        backgroundMesh        
        levelSet
        type
        cutCells
        cM
    end
    
    methods (Access = public)
        
        function c = computeInteriorMesh(obj)
            obj.mesh = obj.cM.computeMesh();
            obj.subcellIsoCoords      = obj.cM.xCoordsIso;
            obj.cellContainingSubcell = obj.cM.cellContainingSubcell;
            obj.computeCutMeshOfSubCellGlobal();
            obj.computeCutMeshOfSubCellLocal('INTERIOR');               
            c = obj.compute();
        end
        
        function c = computeBoundaryMesh(obj)
            obj.mesh = obj.cM.computeBoundaryMesh();
            obj.subcellIsoCoords      = obj.cM.obtainBoundaryXcutIso();
            obj.cellContainingSubcell = obj.cM.obtainBoundaryCellContainingSubCell();
            obj.computeCutMeshOfSubCellGlobal();
            obj.computeCutMeshOfSubCellLocal('BOUNDARY');               
            c = obj.compute();
        end
        
        function c = compute3D(obj)
            obj.mesh                  = obj.cM.computeMesh();
            obj.subcellIsoCoords      = obj.cM.xCoordsIso;
            obj.cellContainingSubcell = obj.cM.cellContainingSubcell;
            obj.computeCutMeshOfSubCellGlobal();
            obj.computeCutMeshOfSubCellLocal(obj.cM.type);               
            c = obj.compute();            
        end
        
        function c = compute(obj)         
            c.cutMeshOfSubCellLocal  = obj.cutMeshOfSubCellLocal;
            c.cutMeshOfSubCellGlobal = obj.cutMeshOfSubCellGlobal;
            c.cellContainingSubcell  = obj.cellContainingSubcell;
            c.mesh                   = obj.mesh;
        end
        
        function obj = CutMesh(cParams)
            obj.init(cParams)
            obj.computeBackgroundCutMesh();
            
            bMtype = obj.backgroundMesh.type;
            isTriangle = isequal(bMtype,'TRIANGLE');
            isQuad     = isequal(bMtype,'QUAD');
            isLine     = isequal(bMtype,'LINE');
            
            if isTriangle
                s.backgroundMesh = obj.backgroundCutMesh;
                s.cutCells       = obj.cutCells;
                s.levelSet       = obj.levelSet;
                obj.cM = CutMeshComputerProvisional(s);
                obj.cM.compute();

            elseif isQuad
                s.backgroundMesh = obj.backgroundCutMesh;
                s.cutCells       = obj.cutCells;
                s.levelSet       = obj.levelSet;
                s.lastNode       = max(obj.backgroundMesh.connec(:));
                obj.cM = CutMeshProvisionalQuadrilater(s);
                obj.cM.compute();
                
            elseif isLine
                s.backgroundMesh = obj.backgroundCutMesh;
                s.cutCells       = obj.cutCells;
                s.levelSet       = obj.levelSet;
                obj.cM = CutMeshProvisionalLine(s);
                obj.cM.compute();
            else
                obj.cM = CutMeshProvisionalOthers(cParams);
                obj.cM.compute();
            end
  
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.levelSet       = cParams.levelSet;
            obj.backgroundMesh = cParams.backgroundMesh;
            obj.cutCells       = cParams.cutCells;
        end
        
        function computeBackgroundCutMesh(obj)
            cells = obj.cutCells;
            m = obj.computeSubMesh(cells);
            obj.backgroundCutMesh = m;
        end
        
        function m = computeSubMesh(obj,cells)
            s.coord  = obj.backgroundMesh.coord;
            s.connec = obj.backgroundMesh.connec(cells,:);
            s.kFace  = obj.backgroundMesh.kFace;
            m = Mesh(s);            
        end
                
        function computeCutMeshOfSubCellGlobal(obj)
            cells = obj.cellContainingSubcell;
            m = obj.computeSubMesh(cells);
            obj.cutMeshOfSubCellGlobal = m;
        end
        
        function computeCutMeshOfSubCellLocal(obj,type)
            coord = obj.subcellIsoCoords;
            nElem = size(coord,3);
            nNode = size(coord,2);
            nDim  = size(coord,1);
            s.coord = reshape(coord,nDim,[])';
            s.connec = reshape(1:nElem*nNode,nNode,nElem)';
            if isequal(type,'INTERIOR')
                kFace = 0;
            else
                kFace = -1;
            end
            s.kFace = kFace;
            m = Mesh(s);
            obj.cutMeshOfSubCellLocal = m;            
        end
        
    end
    
    
end