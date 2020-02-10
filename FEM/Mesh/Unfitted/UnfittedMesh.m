classdef UnfittedMesh < handle
    
    
    properties (GetAccess = public, SetAccess = private)
        innerMesh
        innerCutMesh
        boundaryCutMesh
        
        %Both
        backgroundFullCells
        globalConnec
        unfittedType
        meshBackground
        
        
        %TopOpt
        backgroundCutCells
        subcellIsoCoords
        cellContainingSubcell
        geometryType
        coord
        connec
        
        
        %unfitted
        nActiveBoxFaces
        boxFaceMeshes
        activeBoxFaceMeshesList
        nodesInBoxFaces
        
        oldUnfittedMeshBoundary  
        oldUnfittedMeshInterior
        
        
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
            
            cParams.unfittedType = 'INTERIOR';
            obj.oldUnfittedMeshInterior = Mesh_Unfitted.create2(cParams);
                        
            cParams.unfittedType = 'BOUNDARY';
            obj.oldUnfittedMeshBoundary = Mesh_Unfitted.create2(cParams);            
            
        end
        
        function compute(obj,lvlSet)
            
            obj.oldUnfittedMesh.computeMesh(lvlSet);            
            obj.oldUnfittedMeshInterior.computeMesh(lvlSet)
            obj.oldUnfittedMeshBoundary.computeMesh(lvlSet);
            
            
            %Both
            obj.backgroundFullCells = obj.oldUnfittedMesh.backgroundFullCells;
            
            obj.computeInnerMesh();
            obj.computeInnerCutMesh();
            obj.computeBoundaryCutMesh();            
            
            % caseUnf = 'TopOpt';
            %caseUnf = 'Unfi';
            %caseUnf = '';
            caseUnf = 'Both';
            
            switch  caseUnf
                case 'TopOpt'
                    obj.updateParamsforTopOpt();                    
                case 'Unfi'
                    obj.updateParamsForUnfittedTest();
                case 'Both'
                    obj.updateParamsforTopOpt();
                    obj.updateParamsForUnfittedTest();
            end
            

        end
        
        function plot(obj)
            obj.oldUnfittedMesh.plot();
        end
        
    end
    
    methods (Access = private)
        
        function updateParamsForUnfittedTest(obj)
            
            mesh = obj.oldUnfittedMesh;            
            
            if isprop(mesh,'nActiveBoxFaces')
                obj.nActiveBoxFaces = mesh.nActiveBoxFaces;
            end
            
            if isprop(mesh,'boxFaceMeshes')
                obj.boxFaceMeshes = mesh.boxFaceMeshes;
            end
            
            if isprop(mesh,'activeBoxFaceMeshesList')
                obj.activeBoxFaceMeshesList = mesh.activeBoxFaceMeshesList;
            end
            
            if isprop(mesh,'nodesInBoxFaces')
                obj.nodesInBoxFaces = mesh.nodesInBoxFaces;
            end
        end
        
        function updateParamsforTopOpt(obj)
            mesh = obj.oldUnfittedMesh;
            
            if isprop(mesh,'geometryType')
                obj.geometryType = mesh.geometryType;
            end
            
            if isprop(mesh,'backgroundCutCells')
                obj.backgroundCutCells = mesh.backgroundCutCells;
            end
            
            obj.subcellIsoCoords      = mesh.subcellIsoCoords;
            obj.cellContainingSubcell = mesh.cellContainingSubcell;
            obj.coord  = mesh.coord;
            obj.connec = mesh.connec;
            
        end
        
        function computeInnerMesh(obj)
            obj.computeInnerGlobalConnec();
            s.backgroundCoord = obj.meshBackground.coord;
            s.globalConnec    = obj.globalConnec;
            s.type            = 'INTERIOR';
            obj.innerMesh = InnerMesh(s);
        end
        
        function computeInnerGlobalConnec(obj)
            fullCells = obj.oldUnfittedMeshInterior.backgroundFullCells;
            obj.globalConnec = obj.meshBackground.connec(fullCells,:);
        end
        
        function computeInnerCutMesh(obj)
            s.coord                 = obj.oldUnfittedMeshInterior.coord;
            s.connec                = obj.oldUnfittedMeshInterior.connec;
            % s.type                  = obj.oldUnfittedMesh.typeMesh;            
            
%            cells = obj.oldUnfittedMesh.cellContainingSubcell;
%            connec = obj.meshBackground.connec(cells,:);
            
%            s.coord                 = obj.meshBackground.coord;
%            s.connec                = connec;
            
            s.type = 'INTERIOR';
            
            
            s.backgroundMesh        = obj.meshBackground;
            s.subcellIsoCoords      = obj.oldUnfittedMeshInterior.subcellIsoCoords;
            s.cellContainingSubcell = obj.oldUnfittedMeshInterior.cellContainingSubcell;
            
                      
            obj.innerCutMesh = CutMesh(s);
        end
        
        function computeBoundaryCutMesh(obj)
            cParams.coord                 = obj.oldUnfittedMeshBoundary.coord;
            cParams.connec                = obj.oldUnfittedMeshBoundary.connec;
            cParams.type                  = 'BOUNDARY';
            cParams.backgroundMesh        = obj.meshBackground;
            cParams.subcellIsoCoords      = obj.oldUnfittedMeshBoundary.subcellIsoCoords;
            cParams.cellContainingSubcell = obj.oldUnfittedMeshBoundary.cellContainingSubcell;
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