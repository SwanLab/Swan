classdef UnfittedMesh < handle
    
    
    properties (GetAccess = public, SetAccess = private)
        innerMesh
        innerCutMesh
        boundaryCutMesh
        
        
        %TopOpt
        backgroundCutCells
        geometryType
        
                
        nodesInBoxFaces
        nActiveBoxFaces
        coord
        connec
        subcellIsoCoords
        cellContainingSubcell
        activeBoxFaceMeshesList
        boxFaceMeshes
        
        globalConnec
        backgroundFullCells
        unfittedType
        meshBackground
    end
    
    properties (Access = private)
        oldUnfittedMesh
        type
    end
    
    methods (Access = public)
        
        function obj = UnfittedMesh(cParams)
            obj.meshBackground = cParams.meshBackground;
            obj.oldUnfittedMesh = Mesh_Unfitted.create2(cParams);
            obj.unfittedType = obj.oldUnfittedMesh.unfittedType;

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
            
            obj.backgroundFullCells = obj.oldUnfittedMesh.backgroundFullCells;

            if isprop(obj.oldUnfittedMesh,'geometryType')
                obj.geometryType = obj.oldUnfittedMesh.geometryType;                
            end
            
            if isprop(obj.oldUnfittedMesh,'nActiveBoxFaces')
                obj.nActiveBoxFaces = obj.oldUnfittedMesh.nActiveBoxFaces;                
            end    
            
            if isprop(obj.oldUnfittedMesh,'boxFaceMeshes')
                obj.boxFaceMeshes = obj.oldUnfittedMesh.boxFaceMeshes;                
            end  
            
            if isprop(obj.oldUnfittedMesh,'activeBoxFaceMeshesList')
                obj.activeBoxFaceMeshesList = obj.oldUnfittedMesh.activeBoxFaceMeshesList;                
            end              
            
            if isprop(obj.oldUnfittedMesh,'backgroundCutCells')
                obj.backgroundCutCells = obj.oldUnfittedMesh.backgroundCutCells;                
            end             
            
            if isprop(obj.oldUnfittedMesh,'nodesInBoxFaces')
                obj.nodesInBoxFaces = obj.oldUnfittedMesh.nodesInBoxFaces;                
            end                     
   

            obj.coord  = obj.oldUnfittedMesh.coord;
            obj.connec = obj.oldUnfittedMesh.connec;
            obj.subcellIsoCoords = obj.oldUnfittedMesh.subcellIsoCoords;
            obj.cellContainingSubcell = obj.oldUnfittedMesh.cellContainingSubcell;
            
            obj.computeInnerMesh();
            obj.computeInnerCutMesh();
            obj.computeBoundaryCutMesh();
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
            cParams.subcellIsoCoords = obj.oldUnfittedMesh.subcellIsoCoords;
            cParams.cellContainingSubcell = obj.oldUnfittedMesh.cellContainingSubcell;
            obj.innerCutMesh = CutMesh(cParams);
        end
        
        function computeBoundaryCutMesh(obj)
            cParams.coord  = obj.oldUnfittedMesh.coord;
            cParams.connec = obj.oldUnfittedMesh.connec;
            cParams.type   = 'BOUNDARY';
            cParams.backgroundMesh = obj.meshBackground;
            cParams.subcellIsoCoords = obj.oldUnfittedMesh.subcellIsoCoords;
            cParams.cellContainingSubcell = obj.oldUnfittedMesh.cellContainingSubcell;
            obj.boundaryCutMesh = CutMesh(cParams);
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
    
    
end