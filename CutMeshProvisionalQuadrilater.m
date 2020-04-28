classdef CutMeshProvisionalQuadrilater < handle
    
    properties (Access = public)
        connec
        coord
        cellContainingSubcell
        xCoordsIso
    end
    
    properties (Access = private)
       cutMesh 
       subMesh
       subCutSubMesh
       
       subMesher
       
       fullCells
       cutCells
       
       levelSetSubMesh
    end
    
    properties (Access = private)
        backgroundMesh        
        lastNode
        levelSet
        cutElems
    end
    
    methods (Access = public)
        
        function obj = CutMeshProvisionalQuadrilater(cParams)
            obj.init(cParams);
        end
        
        function compute(obj)
            obj.createSubMesher();            
            obj.createSubMesh();
            obj.computeLevelSetInSubMesh();            
            obj.classifyCells();
            obj.computeCutSubMesh();
            obj.computeSubCutSubMesh();
            obj.computeCoord();
            obj.computeConnec();
            obj.computeCellContainingSubCell();
            obj.computeXcoord();
        end
        
    end   
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.backgroundMesh = cParams.backgroundMesh;
            obj.lastNode = cParams.lastNode;            
            obj.levelSet = cParams.levelSet;
            obj.cutElems = cParams.cutElems;
        end
        
        function createSubMesher(obj)
            s.mesh        = obj.backgroundMesh;
            s.lastNode    = obj.lastNode;
            obj.subMesher = SubMesher(s);            
        end
        
        function createSubMesh(obj)
            obj.subMesh = obj.subMesher.computeSubMesh();
        end
        
        function computeLevelSetInSubMesh(obj)
            ls = obj.levelSet;
            
            s.mesh   = obj.backgroundMesh;
            s.fNodes = ls;
            f = FeFunction(s);
            lsSubMesh = f.computeValueInCenterElement();              
            
            obj.levelSetSubMesh = [ls;lsSubMesh];
        end
        
       function classifyCells(obj)
            lsInElem = obj.computeLevelSetInElem();
            isFull  = all(lsInElem<0,2);
            isEmpty = all(lsInElem>0,2);
            isCut = ~isFull & ~isEmpty;           
            obj.fullCells  = find(isFull);
            obj.cutCells   = find(isCut);
       end
       
        function lsElem = computeLevelSetInElem(obj)
            ls = obj.levelSetSubMesh;
            nodes = obj.subMesh.connec;
            nnode = size(nodes,2);
            nElem = size(nodes,1);
            lsElem = zeros(nElem,nnode);
            for inode = 1:nnode 
                node = nodes(:,inode);
                lsElem(:,inode) = ls(node);
            end
        end        
        
        function computeSubCutSubMesh(obj)
            s.backgroundMesh = obj.computeCutSubMesh();
            s.cutElems = obj.cutCells;
            s.levelSet = obj.levelSetSubMesh;
            cMesh = CutMeshComputerProvisional(s);
            cMesh.compute();
            obj.subCutSubMesh = cMesh;
        end
        
        function m = computeCutSubMesh(obj)
            connecCut = obj.subMesh.connec(obj.cutCells,:);            
            s.coord   = obj.subMesh.coord;
            s.connec = connecCut;
            m = Mesh().create(s);
        end
               
        function  computeCellContainingSubCell(obj)
            cellSubMesh = obj.subCutSubMesh.cellContainingSubcell;
            fCells = obj.fullCells;           
            cellSubMesh = [fCells;cellSubMesh];
            
            cell = obj.computeSubTriangleOfSubCell();
            obj.cellContainingSubcell = cell(cellSubMesh);            
        end
        
        function cell = computeSubTriangleOfSubCell(obj)
            nnode  = size(obj.backgroundMesh.connec,2);  
            cElems = transpose(obj.cutElems);
            cell = repmat(cElems,nnode,1);
        end
            
        function computeXcoord(obj)
            cCells = obj.subCutSubMesh.cellContainingSubcell;   
            fCells = obj.fullCells;
            xIso   = obj.subCutSubMesh.xCoordsIso;
            xIso   = permute(xIso,[2 3 1]);                       
            
            
            nFull = size(fCells,1);
            nCut  = size(cCells,1);
            nElem = nFull + nCut;
            iFull = 1:nFull;
            iCut  = nFull + (1:nCut);
            
            allSubCells = zeros(nElem,1);
            
            allSubCells(iFull,1) = fCells;
            allSubCells(iCut,1)  = cCells;            
            
            bConnec = obj.backgroundMesh.connec;
            nnode   = size(bConnec,2);
            nElem   = size(bConnec,1);
            cell = repmat((1:nnode)',1,nElem);
            globalToLocal = cell(:);            
            
            localSubCells = globalToLocal(allSubCells);
                      
            m = obj.subMesher.computeLocalSubMesh(localSubCells);
            
            xNodalAllIso = m.coordElem;
            
            xIsoCutCells = xIso;
            nDim  = size(xIsoCutCells,1);
            nNode = size(xIsoCutCells,2);
            xIsoAll = zeros(nDim,nNode,nElem);
            xIsoFull = xNodalAllIso(:,:,iFull);
            xIsoAll(:,:,iFull) = xIsoFull;
            xIsoAll(:,:,iCut)  = xIsoCutCells;
            
             
            xE = m.computeXgauss(xIsoAll);  
            xCoords = xE;
            
            obj.xCoordsIso = permute(xCoords,[3 1 2]);
        end
        
        function computeConnec(obj)
            connecCutInterior = obj.subCutSubMesh.connec;
            connecFull = obj.subMesh.connec(obj.fullCells,:);
            obj.connec = [connecFull;connecCutInterior];            
        end
        
        function computeCoord(obj)
            obj.coord  = obj.subCutSubMesh.coord;            
        end
     
    end
    
end