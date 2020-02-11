classdef UnfittedMesh < handle
    
    
    properties (GetAccess = public, SetAccess = private)
        innerMesh
        innerCutMesh
        boundaryCutMesh
        
        backgroundEmptyCells
        
        %TopOpt + Unfitted
        backgroundFullCells
        globalConnec
        unfittedType
        meshBackground
        
        
        %TopOpt
        backgroundCutCells

        cellContainingSubcell
        geometryType
        coord
        connec
        
        
        oldUnfittedMeshBoundary
        
        subcellIsoCoords        

        
    end
    
    properties (Access = private)
        %oldUnfittedMesh
        oldUnfittedMeshInterior
        levelSet
        
        type
    end
    
    methods (Access = public)
        
        function obj = UnfittedMesh(cParams)
            obj.meshBackground = cParams.meshBackground;
            
            obj.unfittedType = cParams.unfittedType;
            
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
            
            obj.levelSet = lvlSet;
            obj.oldUnfittedMeshInterior.computeMesh(lvlSet)
            obj.oldUnfittedMeshBoundary.computeMesh(lvlSet);

            cellsClassifier = CellsClassifier;
            [F,E,C] = cellsClassifier.classifyCells(lvlSet,obj.meshBackground.connec);
            
            
            %1.subcellIsoCoords
            %2.cellContainingSubcell
            %3.backgroundFullCells
            %4.backgroundCutCells
            
            %Both
            obj.backgroundFullCells  = F;
            obj.backgroundEmptyCells = E;            
            obj.backgroundCutCells   = C;
         
            
            obj.computeInnerMesh();
            obj.computeInnerCutMesh();
            obj.computeBoundaryCutMesh();
            obj.computeUnfittedBoxMesh();
            
            obj.updateParamsforTopOpt();

            
        end
        
        function plot(obj)
            switch obj.unfittedType
                case 'BOUNDARY'
                    obj.oldUnfittedMeshBoundary.plot();
                case 'INTERIOR'
                    obj.oldUnfittedMeshInterior.plot();
            end
        end
        
    end
    
    methods (Access = private)
        
        function updateParamsforTopOpt(obj)
            
            switch obj.unfittedType
                case 'BOUNDARY'
                    mesh = obj.oldUnfittedMeshBoundary;
                case 'INTERIOR'
                    mesh = obj.oldUnfittedMeshInterior;
            end
                        
            if isprop(mesh,'geometryType')
                obj.geometryType = mesh.geometryType;
            end
            
            obj.subcellIsoCoords      = mesh.subcellIsoCoords;
            obj.cellContainingSubcell = mesh.cellContainingSubcell;
            obj.coord  = mesh.coord; %Why?? Not necessary
            obj.connec = mesh.connec; %Why?? Not necessary
            
        end
        
        function computeInnerMesh(obj)
            obj.computeInnerGlobalConnec();
            s.backgroundCoord = obj.meshBackground.coord;
            s.globalConnec    = obj.globalConnec;
            s.type            = 'INTERIOR';
            obj.innerMesh = InnerMesh(s);
        end
        
        function computeInnerGlobalConnec(obj)
            fullCells = obj.backgroundFullCells;
            obj.globalConnec = obj.meshBackground.connec(fullCells,:);
        end
        
        function computeInnerCutMesh(obj)
            s.type                  = 'INTERIOR';
            s.coord                 = obj.oldUnfittedMeshInterior.coord;
            s.connec                = obj.oldUnfittedMeshInterior.connec;
            s.backgroundMesh        = obj.meshBackground;
            s.subcellIsoCoords      = obj.oldUnfittedMeshInterior.subcellIsoCoords;
            s.cellContainingSubcell = obj.oldUnfittedMeshInterior.cellContainingSubcell;
            obj.innerCutMesh = CutMesh(s);
        end
        
        function computeBoundaryCutMesh(obj)
            s.type                  = 'BOUNDARY';
            s.coord                 = obj.oldUnfittedMeshBoundary.coord;
            s.connec                = obj.oldUnfittedMeshBoundary.connec;
            s.backgroundMesh        = obj.meshBackground;
            s.subcellIsoCoords      = obj.oldUnfittedMeshBoundary.subcellIsoCoords;
            s.cellContainingSubcell = obj.oldUnfittedMeshBoundary.cellContainingSubcell;
            obj.boundaryCutMesh = CutMesh(s);
        end
        
        function computeUnfittedBoxMesh(obj)
            obj.meshBackground;
        end
        
    end
    
    methods (Access = public)
        
        function m = computeMass(obj)
            switch obj.unfittedType
                case 'BOUNDARY'
                   m = obj.oldUnfittedMeshBoundary.computeMass();
                  %   cParams.mesh = obj;
                  %   cParams.type = obj.unfittedType;
                  %   integrator = Integrator.create(cParams);
                  %   nnodesBackground = size(obj.levelSet);
                  %   fInt = integrator.integrate(ones(nnodesBackground));
                  %   m = sum(fInt);                       
                case 'INTERIOR'
                    %m = obj.oldUnfittedMeshInterior.computeMass();
                    cParams.mesh = obj;
                    cParams.type = obj.unfittedType;
                    integrator = Integrator.create(cParams);
                    nnodesBackground = size(obj.levelSet);
                    fInt = integrator.integrate(ones(nnodesBackground));
                    m = sum(fInt);                    
                    
            end
        end
        
        function aMeshes = getActiveMeshes(obj)
            aMeshes = obj.oldUnfittedMeshInterior.getActiveMeshes();
        end
        
        function add2plot(obj,ax,removedDim,removedCoord)
            switch obj.unfittedType
                case 'BOUNDARY'
                    obj.oldUnfittedMeshBoundary.add2plot(ax,removedDim,removedCoord);
                case 'INTERIOR'
                    obj.oldUnfittedMeshInterior.add2plot(ax,removedDim,removedCoord);
            end
            
        end
        
    end
    
    
end