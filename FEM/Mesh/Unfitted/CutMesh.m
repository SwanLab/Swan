classdef CutMesh < handle
    
    properties (Access = protected)
        backgroundCutMesh     
        backgroundMesh        
        levelSet
        cutCells
    end
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            f = CutMeshFactory();
            obj = f.create(cParams);            
        end
        
    end
    
    methods (Abstract, Access = public)
        compute(obj)
        
    end
    
    methods (Abstract, Access = protected)
        
        obtainMesh(obj)
        obtainXcoordIso(obj)
        obtainCellContainingSubCells(obj)
        
        obtainBoundaryMesh(obj)                
        obtainBoundaryXcutIso(obj)
        obtainBoundaryCellContainingSubCell(obj)
    end
        
    
    methods (Access = public)
        
        function c = computeInteriorMesh(obj)
            xCoord                   = obj.obtainXcoordIso();
            cellContSubCell          = obj.obtainCellContainingSubCells();            
            c.mesh                   = obj.obtainMesh();
            c.cellContainingSubcell  = cellContSubCell;
            c.cutMeshOfSubCellGlobal = obj.computeSubMesh(cellContSubCell);
            c.cutMeshOfSubCellLocal  = obj.computeCutMeshOfSubCellLocal(xCoord,'INTERIOR');               
        end
        
        function c = computeBoundaryMesh(obj)
            xCoord                   = obj.obtainBoundaryXcutIso(); 
            cellContSubCell          = obj.obtainBoundaryCellContainingSubCell();
            c.mesh                   = obj.obtainBoundaryMesh();
            c.cellContainingSubcell  = cellContSubCell;
            c.cutMeshOfSubCellGlobal = obj.computeSubMesh(cellContSubCell);
            c.cutMeshOfSubCellLocal  = obj.computeCutMeshOfSubCellLocal(xCoord,'BOUNDARY');               
        end
        
%         function obj = CutMesh(cParams)
%             obj.init(cParams)
%             obj.computeBackgroundCutMesh();
%             
%             bMtype = obj.backgroundMesh.type;
%             isTriangle = isequal(bMtype,'TRIANGLE');
%             isQuad     = isequal(bMtype,'QUAD');
%             isLine     = isequal(bMtype,'LINE');
%             
%             if isTriangle
%                 s.backgroundMesh = obj.backgroundCutMesh;
%                 s.cutCells       = obj.cutCells;
%                 s.levelSet       = obj.levelSet;
%                 obj.cM = CutMeshComputerProvisional(s);
% 
%             elseif isQuad
%                 s.backgroundMesh = obj.backgroundCutMesh;
%                 s.cutCells       = obj.cutCells;
%                 s.levelSet       = obj.levelSet;
%                 s.lastNode       = max(obj.backgroundMesh.connec(:));
%                 obj.cM = CutMeshProvisionalQuadrilater(s);
%                 
%             elseif isLine
%                 s.backgroundMesh = obj.backgroundCutMesh;
%                 s.cutCells       = obj.cutCells;
%                 s.levelSet       = obj.levelSet;
%                 obj.cM = CutMeshProvisionalLine(s);
%             else
%                 obj.cM = CutMeshProvisionalOthers(cParams);
%             end
%                obj.cM.compute();
%             
%         end
        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.levelSet       = cParams.levelSet;
            obj.backgroundMesh = cParams.backgroundMesh;
            obj.backgroundCutMesh = cParams.backgroundCutMesh;            
            obj.cutCells       = cParams.cutCells;
        end
        
    end
    
    methods (Access = private)
        
%         function computeBackgroundCutMesh(obj)
%             cells = obj.cutCells;
%             m = obj.computeSubMesh(cells);
%             obj.backgroundCutMesh = m;
%         end
        
        function m = computeSubMesh(obj,cells)
            s.coord  = obj.backgroundMesh.coord;
            s.connec = obj.backgroundMesh.connec(cells,:);
            s.kFace  = obj.backgroundMesh.kFace;
            m = Mesh(s);            
        end
                
        function m = computeCutMeshOfSubCellGlobal(obj,cells)
            m = obj.computeSubMesh(cells);
        end
        
    end
    
    methods (Access = private, Static)
        
        function m = computeCutMeshOfSubCellLocal(xCoordIso,type)
            coord = xCoordIso;
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
        end
        
    end
    
    
end