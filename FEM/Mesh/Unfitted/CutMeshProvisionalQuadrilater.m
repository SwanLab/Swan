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
            obj.computeXcoord();            
            obj.computeCoord();
            obj.computeConnec();
            obj.computeCellContainingSubCell();            

        end
        
        function [xCutG,m] = obtainXcutIsoBoundary(obj)
            xCutIso = obj.subCutSubMesh.obtainXcutIso();
           
            s.fullCells     = obj.fullCells;
            s.cutCells      = obj.cutCells;
            s.globalToLocal = obj.computeGlobalToLocal();
            s.localMesh     = obj.subMesher.localMesh;
            s.xIsoCutCoord  = xCutIso;
            xC = XcoordIsoComputer(s); 
            xCutG = xC.computeXSubCut();
            %xCutG = xC.compute();
            
            nElem = size(xCutG,3);
            nNode = size(xCutG,2);
            nDim  = size(xCutG,1);
            sM.coord = reshape(xCutG,nDim,[])';
            sM.connec = reshape(1:nElem*nNode,nNode,nElem)';            
            m = Mesh(sM);            
        end
        
        function m = computeBoundaryMesh(obj)
            m = obj.subCutSubMesh.computeBoundaryMesh();
        end  
        
        function cellCont = obtainBoundaryCellContainingSubCell(obj)
            cutC = obj.cutCells;
            cell = obj.computeSubTriangleOfSubCell();
            cellCont = cell(cutC);
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
            obj.subMesh = obj.subMesher.subMesh;
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
            m = Mesh(s);
        end
               
        function  computeCellContainingSubCell(obj)
            cellSubMesh = obj.subCutSubMesh.cellContainingSubcell;
            fCells      = obj.fullCells;           
            cellSubMesh = [fCells;cellSubMesh];
            
            cell = obj.computeSubTriangleOfSubCell();
            obj.cellContainingSubcell = cell(cellSubMesh);            
        end
        
        function cell = computeSubTriangleOfSubCell(obj)
            nnode  = size(obj.backgroundMesh.connec,2);  
            cElems = transpose(obj.cutElems);
            cell = repmat(cElems,nnode,1);
            cell = cell(:);
        end        
        
        function globalToLocal = computeGlobalToLocal(obj)
            bConnec = obj.backgroundMesh.connec;
            nnode   = size(bConnec,2);
            nElem   = size(bConnec,1);
            cell = repmat((1:nnode)',1,nElem);
            globalToLocal = cell(:);               
        end
        
        function computeXcoord(obj)
            s.fullCells     = obj.fullCells;
            s.cutCells      = obj.subCutSubMesh.cellContainingSubcell;
            s.globalToLocal = obj.computeGlobalToLocal();
            s.localMesh     = obj.subMesher.localMesh;
            s.xIsoCutCoord  = obj.subCutSubMesh.xCoordsIso;
            xC = XcoordIsoComputer(s);
            obj.xCoordsIso = xC.compute();
        end
     
        function computeConnec(obj)
            connecCutInterior = obj.subCutSubMesh.connec;
            connecFull        = obj.subMesh.connec(obj.fullCells,:);
            obj.connec = [connecFull;connecCutInterior];            
        end
        
        function computeCoord(obj)
            obj.coord  = obj.subCutSubMesh.coord;            
        end
     
    end
    
end