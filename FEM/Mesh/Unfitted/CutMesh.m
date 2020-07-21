classdef CutMesh < handle
    
    properties (GetAccess = public, SetAccess = private)
        cutMeshOfSubCellLocal
        cutMeshOfSubCellGlobal
        cellContainingSubcell        
        mesh
    end
    
    properties (Access = private)
        backgroundCutMesh
        subcellIsoCoords
    end
    
    properties (Access = private)
        backgroundMesh        
        levelSet
        type
        cutElems
    end
    
    methods (Access = public)
        
        function obj = CutMesh(cParams)
            
            obj.init(cParams)
            obj.computeBackgroundCutMesh();
            
            bMtype = obj.backgroundMesh.type;
            isTriangle = isequal(bMtype,'TRIANGLE');
            isQuad     = isequal(bMtype,'QUAD');
            isLine     = isequal(bMtype,'LINE');
            
            isInterior = isequal(obj.type,'INTERIOR');
            isBoundary = isequal(obj.type,'BOUNDARY');
            
            if isTriangle
                
                s.backgroundMesh = obj.backgroundCutMesh;
                s.cutElems       = obj.cutElems;
                s.levelSet       = obj.levelSet;
                cM = CutMeshComputerProvisional(s);
                cM.compute();
                
                if isInterior
                    
                    obj.mesh = cM.computeMesh();                    
                    obj.subcellIsoCoords      = cM.xCoordsIso;
                    obj.cellContainingSubcell = cM.cellContainingSubcell;
                    
                elseif isBoundary
                    
                    obj.mesh = cM.computeBoundaryMesh();                    
                    obj.subcellIsoCoords      = cM.obtainXcutIso();
                    obj.cellContainingSubcell = obj.cutElems;
                end
                
            elseif isQuad
                
                s.backgroundMesh = obj.backgroundCutMesh;
                s.cutElems       = obj.cutElems;
                s.levelSet       = obj.levelSet;
                s.lastNode       = max(obj.backgroundMesh.connec(:));
                cM = CutMeshProvisionalQuadrilater(s);
                cM.compute();
                
                if isInterior       
                    obj.mesh = cM.computeMesh();                                        
                    obj.subcellIsoCoords      = cM.xCoordsIso;
                    obj.cellContainingSubcell = cM.cellContainingSubcell;
                elseif isBoundary                    
                    obj.mesh = cM.computeBoundaryMesh();
                    obj.subcellIsoCoords      = cM.obtainXcutIsoBoundary();
                    obj.cellContainingSubcell = cM.obtainBoundaryCellContainingSubCell();                    
                end
                
                
            elseif isLine && isInterior
                s.backgroundMesh = obj.backgroundCutMesh;
                s.cutElems       = obj.cutElems;
                s.levelSet       = obj.levelSet;
                cM = CutMeshProvisionalLine(s);
                cM.compute();
                
                obj.mesh                  = cM.computeMesh();
                obj.subcellIsoCoords      = cM.xCoordsIso;
                obj.cellContainingSubcell = cM.cellContainingSubcell;                
                
            else
                cM = CutMeshProvisionalOthers(cParams);
                cM.compute();
                obj.mesh                  = cM.computeMesh();
                obj.subcellIsoCoords      = cM.xCoordsIso;
                obj.cellContainingSubcell = cM.cellContainingSubcell;
                
            end
            
            obj.computeCutMeshOfSubCellGlobal();
            obj.computeCutMeshOfSubCellLocal();

        end
        
    end
    
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.levelSet       = cParams.levelSet;
            obj.backgroundMesh = cParams.backgroundMesh;
            obj.type           = cParams.type;
            obj.cutElems       = cParams.cutCells;
        end
        
        function computeBackgroundCutMesh(obj)
            cells = obj.cutElems;
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
        
        function computeCutMeshOfSubCellLocal(obj)
            coord = obj.subcellIsoCoords;
            nElem = size(coord,3);
            nNode = size(coord,2);
            nDim  = size(coord,1);
            s.coord = reshape(coord,nDim,[])';
            s.connec = reshape(1:nElem*nNode,nNode,nElem)';
            if isequal(obj.type,'INTERIOR')
                kFace = 0;
            else
                kFace = -1;
            end
            s.kFace = kFace;
            if kFace + obj.backgroundMesh.ndim > 0
                m = Mesh(s);
                obj.cutMeshOfSubCellLocal = m;
            end
        end
        
    end
    
    
end