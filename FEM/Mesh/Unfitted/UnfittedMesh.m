classdef UnfittedMesh < handle
    
    properties (Access = private)
        innerMesh
        innerCutMesh
        boundaryCutMesh        
        unfittedBoxMeshes   
        
        globalConnec

        backgroundFullCells
        backgroundEmptyCells
        backgroundCutCells

        backgroundMesh           
    end
    
    properties (Access = private)
        unfittedType
        isInBoundary     
        levelSet
    end
    
    methods (Access = public)
        
        function obj = UnfittedMesh(cParams)
            obj.init(cParams);
        end
        
        function compute(obj,lvlSet)
            
            obj.levelSet = lvlSet;
            
            cellsClassifier = CellsClassifier;
            [F,E,C] = cellsClassifier.classifyCells(lvlSet,obj.backgroundMesh.connec);
            
            obj.backgroundFullCells  = F;
            obj.backgroundEmptyCells = E;
            obj.backgroundCutCells   = C;
            
            
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
        
        function computeInnerMesh(obj)
            obj.computeInnerGlobalConnec();
            s.backgroundMesh = obj.backgroundMesh;
            s.globalConnec    = obj.globalConnec;
            s.isInBoundary    = obj.isInBoundary;
            obj.innerMesh = InnerMesh(s);
        end
        
        function computeInnerGlobalConnec(obj)
            fullCells = obj.backgroundFullCells;
            obj.globalConnec = obj.backgroundMesh.connec(fullCells,:);
        end
        
        function computeInnerCutMesh(obj)
            s.type                   = 'INTERIOR';
            s.backgroundMesh          = obj.backgroundMesh;
            s.interpolationBackground = Interpolation.create(obj.backgroundMesh,'LINEAR');
            s.backgroundFullCells     = obj.backgroundFullCells;
            s.backgroundEmptyCells    = obj.backgroundEmptyCells;
            s.backgroundCutCells      = obj.backgroundCutCells;
            s.isInBoundary            = obj.isInBoundary;
            s.levelSet = obj.levelSet;
            obj.innerCutMesh = CutMesh(s);
        end
        
        function computeBoundaryCutMesh(obj)
            if ~obj.isInBoundary
                s.type                    = 'BOUNDARY';
                s.backgroundMesh          = obj.backgroundMesh;
                s.interpolationBackground = Interpolation.create(obj.backgroundMesh,'LINEAR');
                s.backgroundFullCells     = obj.backgroundFullCells;
                s.backgroundEmptyCells    = obj.backgroundEmptyCells;
                s.backgroundCutCells      = obj.backgroundCutCells;
                s.isInBoundary            = obj.isInBoundary;
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
        
        function mass = computeMass2(obj)
            npnod = obj.backgroundMesh.npnod;
            f = ones(npnod,1);
            fInt = obj.integrateNodalFunction(f);
            mass = sum(fInt);
            mass2 = obj.computeMass();
        end        
        
        function m = computeMass(obj)
            s.mesh = obj;
            s.type = 'COMPOSITE';
            s = obj.createInteriorParams(s,s.mesh);
            integrator = Integrator.create(s);
            nodalF = ones(size(obj.levelSet));
            fInt = integrator.integrateAndSum(nodalF);
            m = sum(fInt);
        end
        
        function m = computePerimeter(obj)
            cParams.mesh = obj.boundaryCutMesh;
            cParams.type = obj.unfittedType;
            integrator = Integrator.create(cParams);
            nnodesBackground = size(obj.levelSet);
            fInt = integrator.integrate(ones(nnodesBackground));
            m = sum(fInt);
        end
        
        function int = integrateNodalFunction(obj,f)
            uMesh = obj;
            cParams.mesh = uMesh;
            cParams.type = 'COMPOSITE';
            switch obj.unfittedType%    meshType
                case 'INTERIOR'
                    cParams = obj.createInteriorParams(cParams,uMesh);
                case 'BOUNDARY'
                    cParams.compositeParams = cell(0);
                    if ~isempty(uMesh.unfittedBoxMeshes)
                        meshes = uMesh.unfittedBoxMeshes;
                        cParams.boxFaceToGlobal = meshes.nodesInBoxFaces;
                        iActive = 1;
                        for iMesh = 1:length(meshes.isBoxFaceMeshActive)
                            isActive = meshes.isBoxFaceMeshActive(iMesh);
                            if isActive
                                boxFaceMesh = meshes.boxFaceMeshes{iMesh};
                                s = obj.createCompositeParams(boxFaceMesh,uMesh);
                                s.boxFaceToGlobal = meshes.nodesInBoxFaces{iMesh};
                                s.meshBackground = uMesh.backgroundMesh;
                                cParams.compositeParams{iActive} = s;
                                iActive = iActive + 1;
                            end
                        end
                    end
                    cParamsInnerBoundaryCut = obj.createBoundaryCutParams(uMesh);
                    cParams.compositeParams{end+1} = cParamsInnerBoundaryCut;
            end
            integratorC = Integrator.create(cParams);
            int = integratorC.integrateAndSum(f);
        end                
        
    end
    
    methods (Access = private)
        
        function cParams = createInteriorParams(obj,cParams,mesh)
            if mesh.innerCutMesh.nelem ~= 0
                cParamsInnerCut = obj.createInnerCutParams(mesh);
                cParams.compositeParams{1} = cParamsInnerCut;
                if mesh.innerMesh.nelem ~= 0
                    cParamsInner = obj.createInnerParams(mesh);
                    cParams.compositeParams{end+1} = cParamsInner;
                end
            else
                if mesh.innerMesh.nelem ~= 0
                    cParamsInner = obj.createInnerParams(mesh);
                    cParams.compositeParams{1} = cParamsInner;
                else
                    cParams.compositeParams = cell(0);
                end
            end
        end
        
        function cParams = createInnerParams(obj,mesh)
            cParams.mesh = mesh.innerMesh;
            cParams.type = 'SIMPLE';
            cParams.globalConnec = mesh.globalConnec;
            cParams.npnod = mesh.innerMesh.npnod;
            cParams.backgroundMesh = mesh.backgroundMesh;
            cParams.innerToBackground = mesh.backgroundFullCells;
        end
        
        function cParams = createInnerCutParams(obj,mesh)
            cParams.mesh = mesh.innerCutMesh;
            cParams.type = 'CutMesh';
            cParams.meshBackground = mesh.backgroundMesh;
        end
        
        function cParams = createBoundaryCutParams(obj,mesh)
            cParams.mesh = mesh.boundaryCutMesh;
            cParams.type = 'CutMesh';
            cParams.meshBackground = mesh.backgroundMesh;
        end
        
        function cParams = createCompositeParams(obj,mesh,thisMesh)
            cParams.mesh = mesh;
            cParams.type = 'COMPOSITE';
            cParams.meshBackground = obj.backgroundMesh; %%
            cParams.globalConnec = obj.globalConnec;
            cParams.innerToBackground = [];
            cParams.npnod = thisMesh.backgroundMesh.npnod;
            cParams = obj.createInteriorParams(cParams,mesh);
        end
    end
    
    methods (Access = private, Static)
        
        function plotMesh(mesh)
            s.mesh = mesh;
            mP = MeshPlotter(s);
            mP.plot();
        end
        
    end
    
    
end