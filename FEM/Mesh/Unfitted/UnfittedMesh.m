classdef UnfittedMesh < handle
    
    properties (GetAccess = public, SetAccess = private)
        innerMesh
        innerCutMesh        
        boundaryCutMesh       
        unfittedBoxMeshes
    
        %%%%%%ehhh
        backgroundMesh
        
    end
    
    properties (Access = private)
        unfittedType  
        fullCells                        
        cutCells        
        emptyCells
    end
    
    properties (Access = private)
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
        end
        
        function plotBoundary(obj)
            figure
            hold on
            obj.plotMesh(obj.backgroundMesh);
            obj.plotMesh(obj.boundaryCutMesh);
            for imesh = 1:length(obj.unfittedBoxMeshes.isBoxFaceMeshActive)
                if obj.unfittedBoxMeshes.isBoxFaceMeshActive(imesh)
                    uBoxMesh = obj.unfittedBoxMeshes.boxFaceMeshes{imesh};
                    uBoxMesh.plotAll();
                end
            end
        end
        
        function plot(obj)
            figure
            hold on
            obj.plotAll()
        end
        
        function plotAll(obj)
            obj.plotMesh(obj.backgroundMesh);
            obj.plotMesh(obj.innerMesh);
            obj.plotMesh(obj.innerCutMesh);
            obj.plotMesh(obj.boundaryCutMesh);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.backgroundMesh = cParams.meshBackground;
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
            if isequal(class(obj.backgroundMesh),'Mesh_Total')
                m = obj.backgroundMesh;
                fMeshes = m.boxFaceMeshes;
                fNodes  = m.nodesInBoxFaces;
                fGlobalConnec = m.globalConnectivities;
                sides = 2;
                nboxFaces = sides*m.ndim;
                isBoxFaceMeshActive = false([1 nboxFaces]);
                iFace = 0;
                for idime = 1:m.ndim
                    for iside = 1:sides
                        iFace = iFace + 1;
                        mesh = fMeshes{iFace};
                        nodesInBoxFace = fNodes{iFace};
                        interp = Interpolation.create(mesh,'LINEAR');
                        s.type = 'INTERIOR';
                        s.meshBackground = mesh;
                        s.interpolationBackground = interp;
                        cParams = SettingsMeshUnfitted(s);
                        cParams.isInBoundary = true;
                        boxFaceMesh = UnfittedMesh(cParams);
                        
                        
                        lsBoxFace = obj.levelSet(nodesInBoxFace);
                        if any(sign(lsBoxFace)<0)
                            boxFaceMesh.compute(lsBoxFace);
                            isBoxFaceMeshActive(iFace) = true;
                        end
                        
                        boxFaceMeshes{iFace}        = boxFaceMesh;
                        nodesInBoxFaces{iFace}      = nodesInBoxFace;
                        globalConnectivities{iFace} = fGlobalConnec{iFace};
                        
                    end
                end
                obj.unfittedBoxMeshes.boxFaceMeshes        = boxFaceMeshes;
                obj.unfittedBoxMeshes.isBoxFaceMeshActive  = isBoxFaceMeshActive;
                obj.unfittedBoxMeshes.nodesInBoxFaces      = nodesInBoxFaces;
                obj.unfittedBoxMeshes.globalConnectivities = globalConnectivities;
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
    
    methods (Access = private)
        
      
    end
    
    methods (Access = private, Static)
        
        function plotMesh(mesh)
            s.mesh = mesh;
            mP = MeshPlotter(s);
            mP.plot();
        end
        
    end
    
    
end