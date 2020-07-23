classdef Integrator_Unfitted < Integrator
    
    properties (Access = private)
        integrators
    end
    
    methods (Access = public)
        
        function obj = Integrator_Unfitted(cParams)
            obj.mesh = cParams.mesh;
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
        
        function computeInteriorIntegrators(obj)
            s = obj.createInteriorParams(obj.mesh,obj.mesh.backgroundMesh.connec);
            obj.integrators = Integrator.create(s);
        end
        
        function computeBoundaryIntegrators(obj)
            uMesh  = obj.mesh;
            s.type = 'COMPOSITE';      
            s.npnod = uMesh.backgroundMesh.npnod; 
            s.compositeParams = obj.computeBoundaryParams();
            obj.integrators = Integrator.create(s);
        end
        
        function s = computeBoundaryParams(obj)
            s{1} = obj.computeBoundaryCutParams();                    
            [sU,nMeshes] = obj.computeUnfittedBoundaryMeshParams();
            if nMeshes > 0
               s(1+(1:nMeshes)) = sU; 
            end
        end
        
        function [s,nMeshes] = computeUnfittedBoundaryMeshParams(obj)            
            uMesh   = obj.mesh.unfittedBoundaryMesh;        
            uMeshes = uMesh.getActiveMesh();
            gConnec = uMesh.getGlobalConnec();
            nMeshes = numel(uMeshes);
            s = cell(nMeshes,1);
            for iMesh = 1:nMeshes
              s{iMesh} = obj.createInteriorParams(uMeshes{iMesh},gConnec{iMesh});
            end
        end
               
       function s = computeBoundaryCutParams(obj)
            boundaryCutMesh = obj.mesh.boundaryCutMesh;

            gConnec   = obj.mesh.backgroundMesh.connec;            
            connec = gConnec(boundaryCutMesh.cellContainingSubcell,:);
            
            
            sM.coord  = rand(size(obj.mesh.backgroundMesh.coord));
            sM.connec = connec;
            sM.kFace  = obj.mesh.backgroundMesh.kFace;
            backgroundCutMesh = Mesh(sM); 
            
            coord = boundaryCutMesh.xCoordsIso;
            nElem = size(coord,3);
            nNode = size(coord,2);
            nDim  = size(coord,1);
            s.coord = reshape(coord,nDim,[])';
            s.connec = reshape(1:nElem*nNode,nNode,nElem)';
            s.kFace = -1;          
            localCutMesh = Mesh(s);            
            
            s.type                  = 'CutMesh';            
            s.mesh                  = boundaryCutMesh.mesh;
            s.cutMeshOfSubCellLocal = localCutMesh;
            s.cellContainingSubcell = boundaryCutMesh.cellContainingSubcell;            
            s.globalConnec          = obj.mesh.backgroundMesh.connec;
            s.npnod                 = backgroundCutMesh.npnod;
            s.connec                = connec;            
            s.meshType                  = backgroundCutMesh.type;
        end        
        
        function s = createInnerCutParams(obj,gConnec,mesh)
            innerCutMesh = mesh.innerCutMesh;
            
            connec = gConnec(innerCutMesh.cellContainingSubcell,:);
            
            sM.connec = connec;
            sM.coord  = rand(size(obj.mesh.backgroundMesh.coord));
            sM.kFace  = mesh.backgroundMesh.kFace;
            backgroundCutMesh = Mesh(sM);
            
            coord = innerCutMesh.xCoordsIso;
            nElem = size(coord,3);
            nNode = size(coord,2);
            nDim  = size(coord,1);
            s.coord = reshape(coord,nDim,[])';
            s.connec = reshape(1:nElem*nNode,nNode,nElem)';
            s.kFace = 0;
            m = Mesh(s);
            
            
            s.type                  = 'CutMesh';
            s.mesh                  = innerCutMesh.mesh;
            s.cutMeshOfSubCellLocal = m;
            s.cellContainingSubcell = innerCutMesh.cellContainingSubcell;
            s.globalConnec          = gConnec;
            s.npnod                 = backgroundCutMesh.npnod;
            s.connec                = connec;
            s.meshType                  = backgroundCutMesh.type;
        end
        
        function s = createInnerParams(obj,gConnec,mesh)
            s.mesh = mesh;
            s.type = 'SIMPLE';
            s.globalConnec      = gConnec;
            s.npnod             = obj.mesh.backgroundMesh.npnod;
        end        
        
        function s = createInteriorParams(obj,mesh,connec)
            s.type = 'COMPOSITE';
            s.npnod = mesh.backgroundMesh.npnod;
            s.compositeParams = cell(0);
            if ~isempty(mesh.innerMesh)
                fullCells = mesh.innerMesh.fullCells;
                gConnec   = connec(fullCells,:);
                s.compositeParams{1} = obj.createInnerParams(gConnec,mesh.innerMesh.mesh);
            end
            if ~isempty(mesh.innerCutMesh)
                gConnec = connec;
                innerCutParams = obj.createInnerCutParams(gConnec,mesh);
                s.compositeParams{end+1} = innerCutParams;                
            end
        end
        
    end
    
end