classdef UnfittedMesh < handle
    
    properties (GetAccess = public, SetAccess = private)
        unfittedType
        meshBackground
        levelSet_background
        backgroundCutCells
        geometryType
        coord
        connec
        subcellIsoCoords
        cellContainingSubcell
        backgroundFullCells
        
        nActiveMeshes
        nActiveBoxFaces
        activeBoxFaceMeshesList
        boxFaceMeshes
        innerMeshOLD
    end
    
    properties (GetAccess = public, SetAccess = private)
        innerMesh
        innerCutMesh

        globalConnec
        innerMeshConnec
        innerMeshCoord
    end
    
    properties (Access = private)
        oldUnfittedMesh
    end
    
    methods (Access = public)
        
        function obj = UnfittedMesh(cParams)
            obj.oldUnfittedMesh = Mesh_Unfitted.create2(cParams);
        end
        
        function compute(obj,lvlSet)
            obj.oldUnfittedMesh.computeMesh(lvlSet);
            obj.computeInnerMesh();
            obj.computeInnerCutMesh();
        end
        
        function plot(obj)
            obj.oldUnfittedMesh.plot();
        end
        
    end
    
    methods (Access = private)
        
        function computeInnerMesh(obj)
            obj.computeGlobalConnec();
            obj.computeInnerMeshCoords();
            obj.computeInnerMeshConnec();
            obj.innerMesh = Mesh().create(obj.innerMeshCoord,obj.innerMeshCoord);
        end
        
        function computeGlobalConnec(obj)
            fullCells = obj.oldUnfittedMesh.backgroundFullCells;
            obj.globalConnec = obj.meshBackground.connec(fullCells,:);
        end
        
         function computeInnerMeshCoords(obj)
            coordElem = [];
             for inode = 1:obj.meshBackground.nnode
                nodes = obj.globalConnec(:,inode);
                coordElem = [coordElem; obj.meshBackground.coord(nodes,:)];
            end
            subCoords = unique(coordElem,'rows','stable');
            obj.innerMeshCoord = subCoords;
         end  
         
         function computeInnerMeshConnec(obj)
            connec = obj.globalConnec;
            coords = obj.meshBackground.coord;
            subCoords = obj.innerMeshCoord;
            nnode = size(connec,2);
             for inode = 1:nnode 
                    coord = coords(connec(:,inode),:);
                    I = obj.findIndexesComparingCoords(coord,subCoords);                    
                    subConnec(:,inode) = I;
             end
            obj.innerMeshConnec = subConnec;
         end
        
        function computeInnerCutMesh(obj)
            cParams.coord = obj.oldUnfittedMesh.coord;
            cParams.connec = obj.oldUnfittedMesh.connec;
            cParams.backgroundMesh = obj.oldUnfittedMesh.meshBackground;
            cParams.subcellIsoCoords = obj.subcellIsoCoords;
            cParams.cellContainingSubcell = obj.cellContainingSubcell;
            obj.innerCutMesh = CutMesh(cParams);
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
    
    
    methods (Access = public)
        
        function m = computeMass(obj)
            m = obj.oldUnfittedMesh.computeMass();
        end
        
        function aMeshes = getActiveMeshes(obj)
            aMeshes = obj.oldUnfittedMesh.getActiveMeshes();
        end
        
        function add2plot(obj,ax,removedDim,removedCoord)
            obj.oldUnfittedMesh.add2plot(ax,removedDim,removedCoord);
        end
        
    end
    
    
    
    methods
        
        function type = get.unfittedType(obj)
            type = obj.oldUnfittedMesh.unfittedType;
        end
        
        function mB = get.meshBackground(obj)
            mB = obj.oldUnfittedMesh.meshBackground;
        end
        
        function lvlSet = get.levelSet_background(obj)
            lvlSet = obj.oldUnfittedMesh.levelSet_background;
        end
        
        function backCutCells = get.backgroundCutCells(obj)
            backCutCells = obj.oldUnfittedMesh.backgroundCutCells;
        end
        
        function gType = get.geometryType(obj)
            gType = obj.oldUnfittedMesh.geometryType;
        end
        
        function coord = get.coord(obj)
            coord = obj.oldUnfittedMesh.coord;
        end
        
        function connec = get.connec(obj)
            connec = obj.oldUnfittedMesh.connec;
        end
        
        function subcellIsoCoords = get.subcellIsoCoords(obj)
            subcellIsoCoords = obj.oldUnfittedMesh.subcellIsoCoords;
        end
        
        function cellContainingSubcell = get.cellContainingSubcell(obj)
            cellContainingSubcell = obj.oldUnfittedMesh.cellContainingSubcell;
        end
        
        function backgroundFullCells = get.backgroundFullCells(obj)
            backgroundFullCells = obj.oldUnfittedMesh.backgroundFullCells;
        end
        
        function nActiveMeshes = get.nActiveMeshes(obj)
            nActiveMeshes = obj.oldUnfittedMesh.nActiveMeshes;
        end
        
        function nActiveBoxFaces = get.nActiveBoxFaces(obj)
            nActiveBoxFaces = obj.oldUnfittedMesh.nActiveBoxFaces;
        end
        
        function activeBoxFaceMeshesList = get.activeBoxFaceMeshesList(obj)
            activeBoxFaceMeshesList = obj.oldUnfittedMesh.activeBoxFaceMeshesList;
        end
        
        function boxFaceMeshes = get.boxFaceMeshes(obj)
            boxFaceMeshes = obj.oldUnfittedMesh.boxFaceMeshes;
        end
        
        function innerMeshOLD = get.innerMeshOLD(obj)
            innerMeshOLD = obj.oldUnfittedMesh.innerMeshOLD;
        end
                
    end
    
end