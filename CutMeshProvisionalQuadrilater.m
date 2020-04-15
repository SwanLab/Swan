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
            lsSubMesh = obj.subMesher.projectToSubMesh(ls);
            obj.levelSetSubMesh = [ls;lsSubMesh];
        end
        
       function classifyCells(obj)
            lsInElem = obj.computeLevelSetInElem();
            isFull  = all(lsInElem<0,2);
            isEmpty = all(lsInElem>0,2);
            isCut = ~isFull & ~isEmpty;           
            obj.fullCells  = isFull;
            obj.cutCells   = isCut;
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
            s.cutElems = find(obj.cutCells);
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
            fCells = find(obj.fullCells);           
            cellSubMesh = [fCells;cellSubMesh];
            
            cell        = obj.computeSubTriangleOfSubCell();
            obj.cellContainingSubcell = cell(cellSubMesh);            
        end
        
        function cell = computeSubTriangleOfSubCell(obj)
            nnode  = size(obj.backgroundMesh.connec,2);  
            cElems = transpose(obj.cutElems);
            cell = repmat(cElems,nnode,1);
        end
            
        function computeXcoord(obj)
            cellSubMesh = [obj.subCutSubMesh.cellContainingSubcell];   
            fCells = find(obj.fullCells);
            xIso = obj.subCutSubMesh.xCoordsIso;
            xIso = permute(xIso,[2 3 1]);                       
            xCoords = obj.subMesher.projectToIsoSubMesh(cellSubMesh,xIso,fCells); 
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