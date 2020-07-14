classdef Integrator_Unfitted < Integrator
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        integrators
    end
    
    methods (Access = public)
        
        function obj = Integrator_Unfitted(cParams)
            obj.init(cParams);
        end
        
        function int = integrateInDomain(obj,F)
            obj.computeInteriorIntegrators();
            int = obj.integrators.integrateAndSum(F);
        end
        
        function int = integrateInBoundary(obj,F)
            obj.computeBoundaryIntegrators();
            int = obj.integrators.integrateAndSum(F);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
        end
        
        function computeInteriorIntegrators(obj)
            s = obj.createInteriorParams(obj.mesh);
            obj.integrators = Integrator.create(s);
        end
        
        function computeBoundaryIntegrators(obj)
            uMesh  = obj.mesh;
            s = obj.createUnfittedBoxMeshesParams(uMesh);
            s.compositeParams{end+1} = obj.createCutParams(uMesh.boundaryCutMesh,uMesh.backgroundMesh);
            obj.integrators = Integrator.create(s);
        end
        
        function s = createUnfittedBoxMeshesParams(obj,uMesh)            
            s.npnod = uMesh.backgroundMesh.npnod;
            s.type = 'COMPOSITE';
            s.compositeParams = cell(0);
            uBoundary = uMesh.unfittedBoundaryMesh;
            iActive = 1;
            for iMesh = 1:numel(uBoundary.meshes)
                m = uBoundary.meshes{iMesh};
                isActive = m.isBoxFaceMeshActive(iMesh);
                if isActive
                    boxFaceMesh = m.boxFaceMeshes;
                    sB = obj.createInteriorParams(boxFaceMesh);
                    sB.boxFaceToGlobal = m.nodesInBoxFaces;
                    s.compositeParams{iActive} = sB;
                    iActive = iActive + 1;
                end
            end
        end
        
        function s = createInnerParams(obj,mesh,backgroundMesh)
            s.mesh = mesh;
            s.type = 'SIMPLE';
            s.globalConnec      = backgroundMesh.connec(mesh.fullCells,:);
            s.npnod             = backgroundMesh.npnod;
            s.geometryType      = mesh.geometryType;
        end
        
        function s = createCutParams(obj,mesh,backgroundMesh)
            s.mesh = mesh;
            s.type = 'CutMesh';
            s.globalConnec      = backgroundMesh.connec(:,:);
            s.npnod             = backgroundMesh.npnod;
            s.geometryType      = mesh.geometryType;
        end
        
        function s = createInteriorParams(obj,mesh)
            s.type = 'COMPOSITE';
            s.npnod = mesh.backgroundMesh.npnod;
            s.compositeParams = cell(0);
            if mesh.innerMesh.nelem ~= 0
                s.compositeParams{1} = obj.createInnerParams(mesh.innerMesh,mesh.backgroundMesh);
            end
            if mesh.innerCutMesh.nelem ~= 0
                s.compositeParams{end+1} = obj.createCutParams(mesh.innerCutMesh,mesh.backgroundMesh);
            end
        end
        
    end
    
end