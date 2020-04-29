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
            fCells      = obj.fullCells;           
            cellSubMesh = [fCells;cellSubMesh];
            
            cell = obj.computeSubTriangleOfSubCell();
            obj.cellContainingSubcell = cell(cellSubMesh);            
        end
        
        function cell = computeSubTriangleOfSubCell(obj)
            nnode  = size(obj.backgroundMesh.connec,2);  
            cElems = transpose(obj.cutElems);
            cell = repmat(cElems,nnode,1);
        end
        
        
        function xIsoAll = computeXisoAll(obj,iFull,iCut,xNodalAllIso,nElem)
            xIso   = obj.subCutSubMesh.xCoordsIso;
            xIso   = permute(xIso,[2 3 1]);                
            xIsoCutCells = xIso;
            nDim  = size(xIsoCutCells,1);
            nNode = size(xIsoCutCells,2);
            xIsoAll = zeros(nDim,nNode,nElem);
            xIsoFull = xNodalAllIso(:,:,iFull);
            xIsoAll(:,:,iFull) = xIsoFull;
            xIsoAll(:,:,iCut)  = xIsoCutCells;            
        end
        

         function [iFull,iCut,nElemA] = computeIfullIcut(obj,fCells,cCells)
            nFull = size(fCells,1);
            nCut  = size(cCells,1);
            nElemA = nFull + nCut;            
            iFull = false(nElemA,1);
            iCut  = false(nElemA,1);
            iFull(1:nFull,1) = true;
            iCut(nFull + (1:nCut),1) = true;             
         end
        
        function globalToLocal = computeGlobalToLocal(obj)
            bConnec = obj.backgroundMesh.connec;
            nnode   = size(bConnec,2);
            nElem   = size(bConnec,1);
            cell = repmat((1:nnode)',1,nElem);
            globalToLocal = cell(:);               
        end
        
        function allSubCells = computeAllSubCells(obj,nElemA,iFull,iCut,fCells,cCells)            
            allSubCells = zeros(nElemA,1);            
            allSubCells(iFull,1) = fCells;
            allSubCells(iCut,1)  = cCells;              
        end
        
        function localSubCells = computeLocalSubCells(obj,nElemA,iFull,iCut,fCells,cCells)
            allSubCells   = obj.computeAllSubCells(nElemA,iFull,iCut,fCells,cCells);            
            globalToLocal = obj.computeGlobalToLocal();            
            localSubCells = globalToLocal(allSubCells);            
        end
                
        function computeXcoord(obj)
            cCells = obj.subCutSubMesh.cellContainingSubcell;   
            fCells = obj.fullCells;
            
            [iFull,iCut,nElemA] = obj.computeIfullIcut(fCells,cCells);
     
            localSubCells = obj.computeLocalSubCells(nElemA,iFull,iCut,fCells,cCells);

                      
            m = obj.subMesher.computeLocalSubMesh(localSubCells);
            
            xNodalAllIso = m.coordElem;
            
            xIsoAll = obj.computeXisoAll(iFull,iCut,xNodalAllIso,nElemA);
             
            xE = m.computeXgauss(xIsoAll);  
           
            xCoords = xE;            
            obj.xCoordsIso = permute(xCoords,[3 1 2]);
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