classdef CutMeshProvisionalOthers < CutMesh
    
    properties (Access = private)
        localMesh
        subMesh
        cutSubMesh
        fullSubCells
        cutSubCells
        coord
        connec
    end
    
    methods (Access = public)
        
        function obj = CutMeshProvisionalOthers(cParams)
            obj.init(cParams);
        end
        
        function compute(obj)
            obj.computeLocalSubMeshHexahedra();
            obj.computeSubcells();
            obj.classifyCells();
            obj.computeCutSubMesh();
            obj.computeCoords();
            obj.computeConnec();
            obj.computeCellContainingSubCell();
            obj.computeMesh();
            obj.computeBoundaryMesh();
            obj.computeBoundaryXCoordsIso();
            obj.computeBoundaryCellContainingSubCell();
            obj.computeInnerCutMesh();
            obj.computeBoundaryCutMesh();
        end
        
    end
  
    methods (Access = private)

        function computeLocalSubMeshHexahedra(obj)
            Xiso     =  [-1 ,-1, -1;...
                        1, -1, -1;...
                        1, 1, -1;...
                        -1, 1, -1;...
                        -1, -1, 1;...
                        1, -1, 1;...
                        1, 1, 1;...
                        -1, 1, 1;];
            connecIso     = delaunay(Xiso);
            s.coord       = Xiso;
            s.connec      = connecIso;
            obj.localMesh = Mesh.create(s);
        end

        function obj = computeSubcells(obj)
            bCutMesh     = obj.backgroundMesh;
            nelem        = bCutMesh.nelem;
            bCutConnec   = bCutMesh.connec;
            connecIso    = obj.localMesh.connec;
            nElemIso     = size(connecIso,1);
            nnodeSubMesh = size(connecIso,2);
            subConnec    = bCutConnec(:,connecIso');
            subConnec    = reshape(subConnec',[nnodeSubMesh,nelem*nElemIso])';
            s.coord      = bCutMesh.coord;
            s.connec     = subConnec;
            obj.subMesh  = Mesh.create(s);
        end

        function classifyCells(obj)
            lsInElem = obj.computeLevelSetInElem();
            isFull  = all(lsInElem<0,2);
            isEmpty = all(lsInElem>0,2);
            isCut = ~isFull & ~isEmpty;
            obj.fullSubCells  = find(isFull);
            obj.cutSubCells   = find(isCut);
        end

        function lsElem = computeLevelSetInElem(obj)
            ls = obj.levelSet;
            nodes = obj.subMesh.connec;
            nnode = size(nodes,2);
            nElem = size(nodes,1);
            lsElem = zeros(nElem,nnode);
            for inode = 1:nnode
                node = nodes(:,inode);
                lsElem(:,inode) = ls(node);
            end
        end

        function computeCutSubMesh(obj)
            s.backgroundMesh = obj.subMesh;
            s.cutCells       = obj.cutSubCells;
            s.levelSet       = obj.levelSet;
            cMesh = CutMesh.create(s);
            cMesh.compute();
            obj.cutSubMesh = cMesh;
        end

        function globalToLocal = computeGlobalToLocal(obj)
            bConnec   = obj.backgroundMesh.connec;
            locConnec = obj.localMesh.connec;
            nSubCells = size(locConnec,1);
            nElem     = size(bConnec,1);
            cell      = repmat((1:nSubCells)',1,nElem);
            globalToLocal = cell(:);
        end

        function computeCoords(obj)
            s.fullCells     = obj.fullSubCells;
            s.cutCells      = obj.cutSubMesh.innerCutMesh.cellContainingSubcell;
            s.globalToLocal = obj.computeGlobalToLocal();
            s.localMesh     = obj.localMesh;
            s.xIsoCutCoord  = obj.cutSubMesh.innerCutMesh.xCoordsIso;
            xC = XcoordIsoComputer(s);
            obj.xCoordsIso = xC.compute();
            obj.coord      = obj.cutSubMesh.innerCutMesh.mesh.coord;
        end
        
        function computeConnec(obj)
            connecCutInterior = obj.cutSubMesh.innerCutMesh.mesh.connec;
            connecFull        = obj.subMesh.connec(obj.fullSubCells,:);
            obj.connec = [connecFull;connecCutInterior];
        end

        function  computeCellContainingSubCell(obj)
            cellSubMesh = obj.cutSubMesh.innerCutMesh.cellContainingSubcell;
            fCells      = obj.fullSubCells;
            cellSubMesh = [fCells;cellSubMesh];
            cell        = obj.computeSubTetrasOfSubCell();
            obj.cellContainingSubcell = cell(cellSubMesh);
        end

        function cell = computeSubTetrasOfSubCell(obj)
            locConnec = obj.localMesh.connec;
            nSubCells = size(locConnec,1);
            cElems    = transpose(obj.cutCells);
            cell      = repmat(cElems,nSubCells,1);
            cell      = cell(:);
        end

        function computeMesh(obj)
            sM.connec = obj.connec;
            sM.coord  = obj.coord;
            sM.kFace  = obj.backgroundMesh.kFace;
            obj.mesh = Mesh.create(sM);
        end

        function computeBoundaryMesh(obj)
            m = obj.cutSubMesh.boundaryCutMesh.mesh;
            obj.boundaryMesh = m;
        end

        function computeBoundaryXCoordsIso(obj)
            xCutIso = obj.cutSubMesh.xCoordsIsoBoundary;
            s.fullCells     = obj.fullSubCells;
            s.cutCells      = obj.cutSubMesh.boundaryCutMesh.cellContainingSubcell;
            s.globalToLocal = obj.computeGlobalToLocal();
            s.localMesh     = obj.localMesh;
            s.xIsoCutCoord  = xCutIso;
            xC = XcoordIsoComputer(s);
            xCutG = xC.computeXSubCut();
            obj.xCoordsIsoBoundary = xCutG;
        end

        function cellCont = computeBoundaryCellContainingSubCell(obj)
            cutC = obj.cutSubMesh.boundaryCutMesh.cellContainingSubcell;
            cell = obj.computeSubTetrasOfSubCell();
            cellCont = cell(cutC);
            obj.cellContainingSubCellBoundary = cellCont;
        end
        
    end
    
    
end