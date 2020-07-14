classdef UnfittedMesh < handle
    
    properties (GetAccess = public, SetAccess = private)
        innerMesh
        innerCutMesh        
        boundaryCutMesh       
        unfittedBoundaryMesh    
        backgroundMesh
    end
    
    properties (Access = private)
        unfittedType  
        fullCells                        
        cutCells        
        emptyCells    
        plotter
    end
    
    properties (Access = private)
        boundaryMesh        
        isInBoundary
        levelSet        
    end
    
    methods (Access = public)
        
        function obj = UnfittedMesh(cParams)
            obj.init(cParams);
        end
        
        function compute(obj,lSet)
            obj.levelSet = lSet;
            obj.classifyCells();
            obj.computeInnerMesh();
            obj.computeInnerCutMesh();
            obj.computeBoundaryCutMesh();
            obj.computeUnfittedBoxMesh();
            obj.createPlotter()            
        end
        
        function createPlotter(obj)
            s.uMesh = obj;
            obj.plotter = UnfittedMeshPlotter(s);
        end
        
        function plotBoundary(obj)
            obj.plotter.plotBoundary();
        end
        
        function plot(obj)
            obj.plotter.plotDomain();
        end
        
        function plotAll(obj)
            obj.plotter.plotAll();           
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.backgroundMesh = cParams.backgroundMesh;
            obj.boundaryMesh   = cParams.boundaryMesh;
            obj.unfittedType   = cParams.unfittedType;
            obj.isInBoundary   = cParams.isInBoundary;
        end
        
        function classifyCells(obj)            
            nodes         = obj.backgroundMesh.connec;
            allCells(:,1) = 1:obj.backgroundMesh.nelem;            
            lsNodes  = obj.levelSet(nodes);
            isLsNeg  = lsNodes < 0;            
            full  = all(isLsNeg,2);
            empty = all(~isLsNeg,2);
            cut   = ~or(full,empty);            
            obj.fullCells  = allCells(full);
            obj.emptyCells = allCells(empty);
            obj.cutCells   = allCells(cut);
        end       
        
        function computeInnerMesh(obj)
            s.backgroundMesh = obj.backgroundMesh;
            s.isInBoundary   = obj.isInBoundary;
            s.fullCells      = obj.fullCells;
            obj.innerMesh = InnerMesh(s);
        end
                  
        function computeInnerCutMesh(obj)
            s.type                   = 'INTERIOR';
            s.backgroundMesh          = obj.backgroundMesh;
            s.interpolationBackground = Interpolation.create(obj.backgroundMesh,'LINEAR');
            s.fullCells     = obj.fullCells;
            s.emptyCells    = obj.emptyCells;
            s.cutCells      = obj.cutCells;
            s.isInBoundary  = obj.isInBoundary;
            s.levelSet = obj.levelSet;
            obj.innerCutMesh = CutMesh(s);
        end
        
        function computeBoundaryCutMesh(obj)
            if ~obj.isInBoundary
                s.type                    = 'BOUNDARY';
                s.backgroundMesh          = obj.backgroundMesh;
                s.interpolationBackground = Interpolation.create(obj.backgroundMesh,'LINEAR');
                s.fullCells     = obj.fullCells;
                s.emptyCells    = obj.emptyCells;
                s.cutCells      = obj.cutCells;
                s.isInBoundary  = obj.isInBoundary;
                s.levelSet = obj.levelSet;
                obj.boundaryCutMesh = CutMesh(s);
            end
        end
        
        function computeUnfittedBoxMesh(obj)
            s.boundaryMesh = obj.boundaryMesh;            
            obj.unfittedBoundaryMesh = UnfittedBoundaryMesh(s);
            if ~obj.backgroundMesh.isInBoundary
               ls   = obj.levelSet;
               obj.unfittedBoundaryMesh.compute(ls);
            end
        end
        
    end
    
    methods (Access = public)
        
        function mass = computeMass(obj)
            npnod = obj.backgroundMesh.npnod;
            f = ones(npnod,1);
            s.mesh = obj;
            s.type = 'Unfitted';
            integrator = Integrator.create(s);            
            fInt = integrator.integrateInDomain(f);
            %%Now to check IntegrateNodal, later by obj.mesh.computeMass
            mass = sum(fInt);
        end
        
        function mass = computePerimeter(obj)
            npnod = obj.backgroundMesh.npnod;
            f = ones(npnod,1);
            s.mesh = obj;
            s.type = 'Unfitted';
            integrator = Integrator.create(s);   
            fInt = integrator.integrateInBoundary(f);
            %%Now to check IntegrateNodal, later by obj.mesh.computeMass
            mass = sum(fInt);
        end

    end

end