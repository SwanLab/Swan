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
            s.type = 'COMPOSITE';      
            s.npnod = uMesh.backgroundMesh.npnod; 
            s.compositeParams = obj.computeCompositeParams();
            obj.integrators = Integrator.create(s);
        end
        
        function s = computeCompositeParams(obj)
            uMesh  = obj.mesh.unfittedBoundaryMesh;        
            [sUnfitted,nMeshes] = obj.createUnfittedBoundaryMeshParams(uMesh);     
            s = cell(nMeshes+1,1);
            for iMesh = 1:nMeshes
                s{iMesh} = sUnfitted{iMesh};
            end
            sCut = obj.createCutParams(obj.mesh.boundaryCutMesh,obj.mesh.backgroundMesh);            
            s{nMeshes+1} = sCut;            
        end
        
        function [s,nMeshes] = createUnfittedBoundaryMeshParams(obj,uMesh)            
            uMeshes = uMesh.getActiveMesh();
            gConnec = uMesh.getGlobalConnec();
            nMeshes = numel(uMeshes);
            s = cell(nMeshes,1);
            for iMesh = 1:nMeshes
              s{iMesh} = obj.createInteriorParams(uMeshes{iMesh});
              s{iMesh}.boxFaceToGlobal = gConnec{iMesh};
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
            if ~isempty(mesh.innerCutMesh)
                s.compositeParams{end+1} = obj.createCutParams(mesh.innerCutMesh,mesh.backgroundMesh);
            end
        end
        
    end
    
end