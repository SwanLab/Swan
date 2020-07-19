classdef CutMesh < handle
    
    properties (GetAccess = public, SetAccess = private)        
        cutMeshOfSubCellLocal
        cutMeshOfSubCellGlobal
        
        cellContainingSubcell
        backgroundMesh
        mesh
    end
    
    properties (Access = private)
        subcellIsoCoords
              
        
        
        backgroundEmptyCells
        backgroundCutCells
        backgroundFullCells
        
        
        levelSet_unfitted
        
        index1
        index2
        
        subcellsMesher
        cutPointsCalculator

        memoryManager
        
        coord_iso
        coord_global_raw
        connec_local
        cellContainingNodes
        
        ndimUnf
        
        nCutCells
        ndimBackground
        
        levelSet
        backgroundGeomInterpolation
        
        
        coord 
        connec
        type
       
        
    end
    
    methods (Access = public)
        
        
        
        function obj = CutMesh(cParams)

            obj.initCut2(cParams)
            

            isTriangle = isequal(cParams.backgroundMesh.type,'TRIANGLE');
            isQuad     = isequal(cParams.backgroundMesh.type,'QUAD');
            isLine     = isequal(cParams.backgroundMesh.type,'LINE');
            
            isInterior = isequal(obj.type,'INTERIOR');  
            isBoundary = isequal(obj.type,'BOUNDARY');  
            cutElems   = cParams.cutCells; 
            
            if isTriangle && isInterior  
                 
                ls = obj.levelSet;
                connecCut = obj.backgroundMesh.connec(cutElems,:);            
    
                coord = obj.backgroundMesh.coord;     
                
                s.coord  = coord;
                s.connec = connecCut;
                s.kFace = obj.backgroundMesh.kFace;                
                backgroundCutMesh = Mesh(s);
                
                s.backgroundMesh = backgroundCutMesh;
                s.cutElems = cutElems;
                s.levelSet = ls;
                
                cM = CutMeshComputerProvisional(s);
                
                
                
                cM.compute();
               
                obj.connec = cM.connec;
                obj.coord = zeros(size(cM.coord,1),size(obj.backgroundMesh.coord,2));
                obj.coord  = cM.coord;
            %    obj.coord(:,3) = cParams.backgroundMesh.coord(:,3);

                xCoordIso = cM.xCoordsIso;
                obj.subcellIsoCoords = xCoordIso;                 
                
                obj.cellContainingSubcell = cM.cellContainingSubcell;
                
            elseif isQuad && isInterior 
                
                ls = obj.levelSet;
                connecCut = obj.backgroundMesh.connec(cutElems,:);            
    
                coord = obj.backgroundMesh.coord(:,1:2);     
                
                s.coord  = coord;
                s.connec = connecCut;
                backgroundCutMesh = Mesh(s);
                
                s.backgroundMesh = backgroundCutMesh;
                s.cutElems = cutElems;
                s.levelSet = ls;
                s.lastNode = max(obj.backgroundMesh.connec(:));

                cM = CutMeshProvisionalQuadrilater(s);                
                cM.compute();
               
                obj.connec = cM.connec;
                obj.coord = zeros(size(cM.coord,1),size(obj.backgroundMesh.coord,2));
                obj.coord(:,1:2)  = cM.coord;
                
                xCoordIso = cM.xCoordsIso;
                
                obj.subcellIsoCoords = xCoordIso;                                           
                
                obj.cellContainingSubcell = cM.cellContainingSubcell;                
                
            elseif isTriangle && isBoundary 
                               
                ls = obj.levelSet;
                connecCut = obj.backgroundMesh.connec(cutElems,:);                
                coord = obj.backgroundMesh.coord;     
                
                s.coord  = coord;
                s.connec = connecCut;
                s.kFace = obj.backgroundMesh.kFace;
                backgroundCutMesh = Mesh(s);                
                
                s.backgroundMesh = backgroundCutMesh;
                s.cutElems = cutElems;
                s.levelSet = ls;
                cutMesh = CutMeshComputerProvisional(s);
                cutMesh.compute();
                
                bMesh = cutMesh.computeBoundaryMesh();
                
                obj.connec = bMesh.connec;
                obj.coord  = bMesh.coord;
                obj.subcellIsoCoords = cutMesh.obtainXcutIso();
                obj.cellContainingSubcell = cutElems;

            elseif isQuad && isBoundary            
                

                ls = obj.levelSet;
                connecCut = obj.backgroundMesh.connec(cutElems,:);            
    
                coord = obj.backgroundMesh.coord(:,1:2);     
                
                s.coord  = coord;
                s.connec = connecCut;
                backgroundCutMesh = Mesh(s);
                
                s.backgroundMesh = backgroundCutMesh;
                s.cutElems = cutElems;
                s.levelSet = ls;
                s.lastNode = max(obj.backgroundMesh.connec(:));

                cM = CutMeshProvisionalQuadrilater(s);                
                cM.compute();
               
 
                
                mt = cM.computeBoundaryMesh();
                conn2 = mt.connec;
                coor2 = mt.coord;
                subCel2 = cM.obtainXcutIsoBoundary();
                cellC2 = cM.obtainBoundaryCellContainingSubCell();
                
                
                obj.coord  = coor2;
                obj.connec = conn2;
                obj.subcellIsoCoords = subCel2;
                obj.cellContainingSubcell = cellC2;
                
            elseif isLine && isInterior
                coord = obj.backgroundMesh.coord(:,1:2);     
                connecCut = obj.backgroundMesh.connec(cutElems,:);            

                s.coord  = coord;
                s.connec = connecCut;
                s.kFace = obj.backgroundMesh.kFace;
                backgroundCutMesh = Mesh(s);
                
                s.backgroundMesh = backgroundCutMesh;
                s.cutElems = cutElems;
                s.levelSet = obj.levelSet;
                
                cM = CutMeshProvisionalLine(s);
                
                cM.compute();
                
                sM.connec = cM.connec;
                sM.coord  = cM.coord;
                sM.kFace = obj.backgroundMesh.kFace;                
                m = Mesh(sM);
                
                
                obj.connec = cM.connec;
                obj.coord  = cM.coord;
                obj.subcellIsoCoords = cM.xCoordsIso;
                obj.cellContainingSubcell = cM.cellContainingSubcell;                
                
            else
                        
                obj.initCut(cParams);
                obj.build();
                
                if obj.isLevelSetCrossingZero()
                    obj.computeCutMesh();
                else
                    obj.returnNullMesh();
                end
                
                obj.subcellIsoCoords = permute(obj.subcellIsoCoords,[3 2 1]);                 
                
            end
            
            
            
            obj.computeCutMeshOfSubCellGlobal();
            obj.computeCutMeshOfSubCellLocal();
            
            sM.connec = obj.connec;
            sM.coord  = obj.coord;   
            if isequal(obj.type,'INTERIOR')
               sM.kFace = obj.backgroundMesh.kFace;
            else
               sM.kFace = obj.backgroundMesh.kFace -1;
            end
            obj.mesh = Mesh(sM);
        end
        
    end
    
    
    methods (Access = private)
        
        function createNdimUnf(obj)
            obj.ndimUnf = obj.backgroundMesh.geometryType;
        end
        
        function initCut(obj,cParams)
            obj.backgroundFullCells   = cParams.fullCells;
            obj.backgroundEmptyCells  = cParams.emptyCells;
            obj.backgroundCutCells    = cParams.cutCells;
            obj.nCutCells             = length(obj.backgroundCutCells);            
            obj.backgroundMesh              = cParams.backgroundMesh;
            obj.backgroundGeomInterpolation = cParams.interpolationBackground;
        end
        
        function initCut2(obj,cParams)
            obj.levelSet = cParams.levelSet;                    
            obj.backgroundMesh = cParams.backgroundMesh;
            obj.type           = cParams.type;
        end        
        
        function createSubCellsMesher(obj)
            sS.ndimIso = obj.ndimUnf;
            sS.type = obj.type;            
            sS.posNodes           = obj.backgroundGeomInterpolation.pos_nodes;
            sS.levelSetBackground = obj.levelSet;
            sS.coordsBackground   = obj.backgroundMesh.coord;
            obj.subcellsMesher    = SubcellsMesher.create(sS);
        end
        
        function m = computeCutMeshOfSubCellGlobal(obj)
            coord  = obj.backgroundMesh.coord; 
            connec = obj.backgroundMesh.connec;
            cells  = obj.cellContainingSubcell;            
            s.coord  = coord;
            s.connec = connec(cells,:);              
            s.kFace = obj.backgroundMesh.kFace;            
            m = Mesh(s); 
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
        
        function createMemoryManager(obj)
            s.ndimIso       = obj.ndimUnf;
            s.unfittedType  = obj.type;
            s.nCutCells     = obj.nCutCells;
            s.ndim          = obj.backgroundMesh.ndim;
            obj.memoryManager  = MemoryManager_MeshUnfitted(s);
        end
        
        function build(obj)
            obj.createNdimUnf();
            obj.createSubCellsMesher();
            obj.createMemoryManager();
            obj.cutPointsCalculator  = CutPointsCalculator;
        end
        
        function itIs = isLevelSetCrossingZero(obj)
            itIs = ~isempty(obj.backgroundCutCells);
        end
        
        function computeCutMesh(obj)
            obj.computeSubcells();
            obj.computeGlobalUnfittedMesh();
        end
        
        function computeGlobalUnfittedMesh(obj)
            obj.computeGlobalCoordinates();
            obj.computeGlobalConnectivities();
        end
        
        function obj = computeSubcells(obj)
            obj.memoryManager.allocateMemory();
            
            obj.computeCutPoints();
            
            for icut = 1:obj.nCutCells %Vectorize
                icell = obj.backgroundCutCells(icut);
                
                newSubcells = obj.computeThisCellSubcells(icut,icell);
                
                newCellContainingNodes   = repmat(icell,[newSubcells.nNodes 1]);
                newCellContainingSubcell = repmat(icell,[newSubcells.nSubcells 1]);
                
                nnode = (newSubcells.nNodes);
                obj.memoryManager.saveNewSubcells(newSubcells,newCellContainingNodes,newCellContainingSubcell);
            end
            obj.memoryManager.freeSpareMemory();
            
            
            obj.coord_iso             = obj.memoryManager.coord_iso;
            obj.coord_global_raw      = obj.memoryManager.coord_global_raw;
            obj.subcellIsoCoords      = obj.memoryManager.subcellIsoCoords;
            obj.connec_local          = obj.memoryManager.connec_local;
            obj.connec                = obj.memoryManager.connec;
            obj.levelSet_unfitted     = obj.memoryManager.levelSet_unfitted;
            obj.cellContainingNodes   = obj.memoryManager.cellContainingNodes;
            obj.cellContainingSubcell = obj.memoryManager.cellContainingSubcell;
            
        end
        
        function computeCutPoints(obj)
            s.meshBackground              = obj.backgroundMesh;
            s.levelSet_background         = obj.levelSet;
            s.backgroundCutCells          = obj.backgroundCutCells;
            s.backgroundGeomInterpolation = obj.backgroundGeomInterpolation;
            obj.cutPointsCalculator.init(s);
            obj.cutPointsCalculator.computeCutPoints();
        end
        
        function subcells = computeThisCellSubcells(obj,icut,icell)
            cutPoints_thisCell = obj.cutPointsCalculator.getThisCellCutPoints(icut);
            connec_thisCell = obj.backgroundMesh.connec(icell,:);
            
            
            sS.cellConnec = connec_thisCell;
            sS.cutPoints = cutPoints_thisCell;
            
            obj.subcellsMesher.computeSubcells(sS);
            
            subcells = obj.subcellsMesher.subcells;
        end
        
        function computeGlobalConnectivities(obj)
            nSubcells = size(obj.connec_local,1);
            cellOfSubCell = obj.cellContainingSubcell;
            coordGlobal   = obj.coord_global_raw;
            connecLocal   = obj.connec_local;
            cellContNodes = obj.cellContainingNodes;
            nnode = size(connecLocal,2);
            connecV = zeros(nSubcells,nnode);
            
            
            for isub = 1:nSubcells %Vectorize !!!
                cell = cellContNodes == cellOfSubCell(isub);
                % size(find(cell))
                coordsSubCell = coordGlobal(cell,:);
                indexes = obj.findIndexesComparingCoords(coordsSubCell,obj.coord);
                
                connecV(isub,:) = indexes(connecLocal(isub,:));
                
                
            end
            obj.connec = connecV;
        end
        
        function computeGlobalCoordinates(obj)
            [coord,ind1,ind2] = unique(obj.coord_global_raw,'rows','stable');
            obj.coord = coord;
            obj.index1 = ind1;
            obj.index2 = ind2;
        end
        
        function returnNullMesh(obj)
            obj.coord = zeros(0,obj.ndim);
            obj.connec = [];
        end
        
    end
    
    methods (Access = private, Static)
        
        function I = findIndexesComparingCoords(A,B)
            I = zeros(1,size(A,1));
            for inode = 1:size(A,1)
                match = true(size(B,1),1);
                for idime = 1:size(A,2)
                    match = match & B(:,idime) == A(inode,idime);
                end
                I(inode) = find(match,1);
            end
        end
        
    end
    
end