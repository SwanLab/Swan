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
        
        
        
        subcellIsoCoords        

        
    end
    
    properties (Access = private)
        oldUnfittedMeshInterior
        oldUnfittedMeshBoundary        
        
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
                    %obj.oldUnfittedMeshBoundary.plot();
                    figure
                    hold on                    
                    obj.boundaryCutMesh.plot;
                    light
                    axis equal off
                    hold off                      
                case 'INTERIOR'
                    figure;
                    hold on
                    obj.innerMesh.plot;
                    obj.innerCutMesh.plot;
                    light
                    axis equal off
                    hold off                    
            end
        end
        
    end
    
    methods (Access = private)
        
        function updateParamsforTopOpt(obj)
            
            switch obj.unfittedType
                case 'BOUNDARY'
                    mesh = obj.boundaryCutMesh;
                case 'INTERIOR'
                    mesh = obj.innerCutMesh;                    
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
            s.unfittedType            = 'INTERIOR';
            s.meshBackground          = obj.meshBackground;
            s.interpolationBackground = Interpolation.create(obj.meshBackground,'LINEAR');
            s.backgroundFullCells     = obj.backgroundFullCells;
            s.backgroundEmptyCells    = obj.backgroundEmptyCells;
            s.backgroundCutCells      = obj.backgroundCutCells;     
            s.levelSet = obj.levelSet;
            obj.innerCutMesh = CutMesh(s);


        end
        
        function computeBoundaryCutMesh(obj)
            s.unfittedType            = 'BOUNDARY';
            s.meshBackground          = obj.meshBackground;
            s.interpolationBackground = Interpolation.create(obj.meshBackground,'LINEAR');
            s.backgroundFullCells     = obj.backgroundFullCells;
            s.backgroundEmptyCells    = obj.backgroundEmptyCells;
            s.backgroundCutCells      = obj.backgroundCutCells;     
            s.levelSet = obj.levelSet;
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
                     cParams.mesh = obj.boundaryCutMesh;
                     cParams.type = obj.unfittedType;
                     integrator = Integrator.create(cParams);
                     nnodesBackground = size(obj.levelSet);
                     fInt = integrator.integrate(ones(nnodesBackground));
                     m = sum(fInt);                       
                case 'INTERIOR'
                    cParams.mesh = obj;
                    cParams.type = obj.unfittedType;
                    integrator = Integrator.create(cParams);
                    nnodesBackground = size(obj.levelSet);
                    fInt = integrator.integrate(ones(nnodesBackground));
                    m = sum(fInt);                                        
            end
        end

        function add2plot(obj,ax,removedDim,removedCoord)
            switch obj.unfittedType
                case 'BOUNDARY'
                    obj.oldUnfittedMeshBoundary.add2plot(ax,removedDim,removedCoord);
                case 'INTERIOR'
                    %obj.oldUnfittedMeshInterior.add2plot(ax,removedDim,removedCoord);
                    
                    %obj.boundaryCutMesh.add2plot(ax);
                    obj.innerMesh.add2plot(ax);
                    obj.innerCutMesh.add2plot(ax,removedDim,removedCoord);
                    
            end
            
        end
        
    end
    
    
end