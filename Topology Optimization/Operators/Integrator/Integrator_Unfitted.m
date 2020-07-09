classdef Integrator_Unfitted < Integrator
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        integratorC
    end
    
    properties (Access = private)
       unfittedMesh 
    end
    
    methods (Access = public)
        
        function obj = Integrator_Unfitted(cParams)
            obj.init(cParams);
        end
        
        function int = integrateInDomain(obj,F)
            obj.computeInteriorIntegrators();                    
            int = obj.integratorC.integrateAndSum(F);                        
        end
        
        function int = integrateInBoundary(obj,F)
            obj.computeBoundaryIntegrators();                    
            int = obj.integratorC.integrateAndSum(F);                 
        end
        
    end
    
    methods (Access = private)
               
        function computeInteriorIntegrators(obj)
            uMesh = obj.mesh;
            cParams.mesh = uMesh;
            cParams.type = 'COMPOSITE';
            cParams = obj.createInteriorParams(cParams,uMesh);
            obj.integratorC = Integrator.create(cParams);                                    
        end
        
        function computeBoundaryIntegrators(obj)
            uMesh = obj.mesh;
            cParams.mesh = uMesh;
            cParams.type = 'COMPOSITE';
            cParams.compositeParams = cell(0);
            if ~isempty(uMesh.unfittedBoxMeshes)
                meshes = uMesh.unfittedBoxMeshes;
                cParams.boxFaceToGlobal = meshes.nodesInBoxFaces;
                iActive = 1;
                for iMesh = 1:length(meshes.isBoxFaceMeshActive)
                    isActive = meshes.isBoxFaceMeshActive(iMesh);
                    if isActive
                        boxFaceMesh = meshes.boxFaceMeshes{iMesh};
                        s = obj.createCompositeParams(boxFaceMesh);
                        s.boxFaceToGlobal = meshes.nodesInBoxFaces{iMesh};
                        cParams.compositeParams{iActive} = s;
                        iActive = iActive + 1;
                    end
                end
            end
            cParamsInnerBoundaryCut = obj.createBoundaryCutParams(uMesh);
            cParams.compositeParams{end+1} = cParamsInnerBoundaryCut;
            obj.integratorC = Integrator.create(cParams);                                    
        end
        
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
            cParams.globalConnec      = mesh.connecFullCells;
            cParams.npnod             = mesh.innerMesh.npnod;
            cParams.backgroundMesh    = mesh.backgroundMesh;
            cParams.innerToBackground = mesh.fullCells;
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
        
        function cParams = createCompositeParams(obj,mesh)
            cParams.mesh = mesh;
            cParams.type = 'COMPOSITE';
            cParams.meshBackground    = obj.mesh.backgroundMesh; %%
            cParams.globalConnec      = obj.mesh.connecFullCells;
            cParams.innerToBackground = [];
            cParams.npnod = obj.mesh.backgroundMesh.npnod;
            cParams = obj.createInteriorParams(cParams,mesh);
            cParams.meshBackground = obj.mesh.backgroundMesh;            
        end        
        
    end
    
end