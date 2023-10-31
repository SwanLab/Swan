classdef RHSintegrator_ShapeFunctionUnfitted < handle

    properties (Access = private)
        globalConnec
        mesh
        integrators
        quadType
    end

    methods (Access = public)

        function obj = RHSintegrator_ShapeFunctionUnfitted(cParams)
            obj.init(cParams);
        end

        function int = integrateInDomain(obj, unfFun, test)
            noIMesh  = isempty(unfFun.unfittedMesh.innerMesh);
            noICMesh = isempty(unfFun.unfittedMesh.innerCutMesh);
            if noIMesh && noICMesh
                int = zeros(size(ls));
            else
                obj.createInteriorIntegrators(test);
                int = obj.integrators.integrateAndSum(unfFun);
            end
        end

        function int = integrateInBoundary(obj,F,test)
            obj.createBoundaryIntegrators(test);
            int = obj.integrators.integrateAndSum(F);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.quadType = cParams.quadType;
        end

        function createInteriorIntegrators(obj,test)
            s = obj.createInteriorParams(obj.mesh,obj.mesh.backgroundMesh.connec);
            s.quadType = obj.quadType;
            s.test = test;
            obj.integrators = RHSintegrator.create(s);
        end
        
        function s = createInteriorParams(obj,mesh,connec)
            s.type = 'Composite';
            s.npnod = obj.mesh.backgroundMesh.nnodes;
            s.compositeParams = cell(0);
            s.unfittedMesh = mesh;
            if ~isempty(mesh.innerMesh)
                s.compositeParams{1} = obj.createInnerParams(mesh.innerMesh);
            end
            if ~isempty(mesh.innerCutMesh)
                gConnec = connec;
                innerCutParams = obj.createInnerCutParams(gConnec,mesh);
                s.compositeParams{end+1} = innerCutParams;
            end
        end

        function s = createInnerParams(obj,innerMesh)
            s.type = 'ShapeFunction';
            s.mesh = innerMesh.mesh;
        end

        function s = createInnerCutParams(obj,gConnec,mesh)
            innerCutMesh = mesh.innerCutMesh;
            s.type                  = 'CutMesh';
            s.mesh                  = innerCutMesh.mesh;
            s.xCoordsIso            = innerCutMesh.xCoordsIso;
            s.cellContainingSubcell = innerCutMesh.cellContainingSubcell;
            s.globalConnec          = gConnec;
            s.npnod                 = obj.mesh.backgroundMesh.nnodes;
            s.backgroundMeshType    = mesh.backgroundMesh.type;
        end

        function createBoundaryIntegrators(obj,test)
            uMesh  = obj.mesh;
            s.type = 'Composite';
            s.npnod = uMesh.backgroundMesh.nnodes;
            s.unfittedMesh = uMesh;
            s.compositeParams = obj.createBoundaryParams();
            s.quadType = 'LINEAR';
            s.test = test;
            obj.integrators = RHSintegrator.create(s);
        end

        function s = createBoundaryParams(obj)
            gConnec = obj.mesh.backgroundMesh.connec;
            s{1} = obj.createBoundaryCutParams(gConnec);
            [sU,nMeshes] = obj.createUnfittedBoundaryMeshParams();
            if nMeshes > 0
                s(1+(1:nMeshes)) = sU;
            end
        end

        function [s,nMeshes] = createUnfittedBoundaryMeshParams(obj)
            uMesh   = obj.mesh.unfittedBoundaryMesh;
            uMeshes = uMesh.getActiveMesh();
            gConnec = uMesh.getGlobalConnec();
            nMeshes = numel(uMeshes);
            s = cell(nMeshes,1);
            for iMesh = 1:nMeshes
                s{iMesh} = obj.createInteriorParams(uMeshes{iMesh},gConnec{iMesh});
            end
        end

        function s = createBoundaryCutParams(obj,gConnec)
            boundaryCutMesh = obj.mesh.boundaryCutMesh;
            s.type                  = 'CutMesh';
            s.mesh                  = boundaryCutMesh.mesh;
            s.xCoordsIso            = boundaryCutMesh.xCoordsIso;
            s.cellContainingSubcell = boundaryCutMesh.cellContainingSubcell;
            s.globalConnec          = gConnec;
            s.npnod                 = obj.mesh.backgroundMesh.nnodes;
            s.backgroundMeshType    = obj.mesh.backgroundMesh.type;
        end

    end

end