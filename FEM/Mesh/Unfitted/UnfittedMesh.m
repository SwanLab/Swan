classdef UnfittedMesh < handle
    
    properties (GetAccess = public, SetAccess = private)
        innerMesh
        innerCutMesh
        boundaryCutMesh
        unfittedBoundaryMesh
        backgroundMesh
    end
    
    properties (Access = private)
        cutMesh
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
            if ~isempty(obj.cutCells)                 
                obj.computeCutMesh();
                obj.computeInnerCutMesh();
                obj.computeBoundaryCutMesh();
            else 
                obj.cutMesh         = [];
                obj.innerCutMesh    = [];
                obj.boundaryCutMesh = [];
            end
            obj.computeUnfittedBoxMesh();            
            obj.createPlotter();            
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
        
        function plotNormals(obj)
            if ~isempty(obj.boundaryCutMesh)
                obj.boundaryCutMesh.mesh.plotNormals();
            end            
        end
        
        function dv = computeDvolume(obj,quad)
            if ~isempty(obj.innerMesh)            
                dvI = obj.computeInnerDvolume(quad);
            else 
                dvI = 0;
            end
            dvC = obj.computeInnerCutDvolume(quad);    
            dv = dvC + dvI;
        end
        
        
    end
    
    methods (Access = private)
        
        function dvolume = computeInnerDvolume(obj,quad)
            nelem = obj.backgroundMesh.nelem;
            ngaus = quad.ngaus;
            dv    = obj.innerMesh.mesh.computeDvolume(quad);
            iCells = obj.innerMesh.fullCells;
            dvolume = zeros(nelem,ngaus);
            dvolume(iCells,:) = dv;
        end
        
        function dvolume = computeInnerCutDvolume(obj,quad)
            nelem = obj.backgroundMesh.nelem;
            ngaus = quad.ngaus;  
            dvolume = zeros(nelem,ngaus);   
            if ~isempty(obj.innerCutMesh)
            cCells = obj.innerCutMesh.cellContainingSubcell;                
            dv = obj.innerCutMesh.mesh.computeDvolume(quad);                        
            for igaus = 1:ngaus
               dvolume(:,igaus) = accumarray(cCells,dv(igaus,:)',[nelem,1],@sum,0);            
            end
            end
        end
        
        function init(obj,cParams)
            obj.backgroundMesh = cParams.backgroundMesh;
            obj.boundaryMesh   = cParams.boundaryMesh;
        end
        
        function classifyCells(obj)
            allCells(:,1) = 1:obj.backgroundMesh.nelem;
            lsNodes  = obj.computelsNodes();
            isLsNeg  = lsNodes < 0;
            full  = all(isLsNeg,2);
            empty = all(~isLsNeg,2);
            cut   = ~or(full,empty);
            obj.fullCells  = allCells(full);
            obj.emptyCells = allCells(empty);
            obj.cutCells   = allCells(cut);
        end
        
        function ls = computelsNodes(obj)
            nodes = obj.backgroundMesh.connec;
            ls = zeros(size(nodes));
            for iNode = 1:size(nodes,2)
                ls(:,iNode) = obj.levelSet(nodes(:,iNode));            
            end
        end
        
        function computeInnerMesh(obj)
            if ~isempty(obj.fullCells)
                s.backgroundMesh = obj.backgroundMesh;
                s.fullCells      = obj.fullCells;
                obj.innerMesh = InnerMesh(s);
            else
                obj.innerMesh = [];
            end
        end
        
        function computeCutMesh(obj)
            s.backgroundMesh = obj.backgroundMesh;
            s.cutCells       = obj.cutCells;
            s.levelSet       = obj.levelSet;
            obj.cutMesh = CutMesh.create(s);
            obj.cutMesh.compute();
        end
        
        function computeInnerCutMesh(obj)
            obj.innerCutMesh = obj.cutMesh.innerCutMesh;
        end
        
        function computeBoundaryCutMesh(obj)
            if ~isequal(obj.backgroundMesh.geometryType,'Line')
                obj.boundaryCutMesh = obj.cutMesh.boundaryCutMesh;
            end
            
        end
        
        function computeUnfittedBoxMesh(obj)
            s.boundaryMesh = obj.boundaryMesh;
            obj.unfittedBoundaryMesh = UnfittedBoundaryMesh(s);
            if ~isempty(obj.unfittedBoundaryMesh.meshes)
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
            %disp('Interior')
            %sum(fInt<0)/size(fInt,1)
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
            %disp('Boundary')
            %sum(fInt<0)/size(fInt,1)           
            mass = sum(fInt);
        end
        
    end
    
end