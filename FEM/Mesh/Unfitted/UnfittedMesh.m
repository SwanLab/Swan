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
        nodesInBoxFaces
    end
    
    properties (GetAccess = public, SetAccess = private)
        innerMesh
        innerCutMesh
        
        globalConnec
    end
    
    properties (Access = private)
        oldUnfittedMesh
        type
        
    end
    
    methods (Access = public)
        
        function obj = UnfittedMesh(cParams)
            obj.meshBackground = cParams.meshBackground;
            obj.oldUnfittedMesh = Mesh_Unfitted.create2(cParams);
            if isobject(cParams)
                if (isempty(cParams.type))
                obj.type = 'INTERIOR';
                else 
                    obj.type = cParams.type;
                end
            else
            if isfield(cParams,'type') 
                obj.type = cParams.type;
            else
                obj.type = 'INTERIOR';
            end
            end
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
                obj.computeInnerGlobalConnec();
                s.backgroundCoord = obj.meshBackground.coord;
                s.globalConnec = obj.globalConnec;
                s.type = obj.type;
                obj.innerMesh = InnerMesh(s);
            end
            
            function computeInnerGlobalConnec(obj)
                fullCells = obj.oldUnfittedMesh.backgroundFullCells;
                obj.globalConnec = obj.meshBackground.connec(fullCells,:);
            end
            
            function computeInnerCutMesh(obj)
                cParams.coord  = obj.oldUnfittedMesh.coord;
                cParams.connec = obj.oldUnfittedMesh.connec;
                cParams.type   = obj.oldUnfittedMesh.typeMesh;
                cParams.backgroundMesh = obj.meshBackground;
                cParams.subcellIsoCoords = obj.subcellIsoCoords;
                cParams.cellContainingSubcell = obj.cellContainingSubcell;
                obj.innerCutMesh = CutMesh(cParams);
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
                %obj.oldUnfittedMesh.add2plot(ax);
            end
            
        end
        
        methods
            
            function type = get.unfittedType(obj)
                type = obj.oldUnfittedMesh.unfittedType;
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
            
            function nodesInBoxFaces = get.nodesInBoxFaces(obj)
                nodesInBoxFaces = obj.oldUnfittedMesh.nodesInBoxFaces;
            end
            
        end
        
    end