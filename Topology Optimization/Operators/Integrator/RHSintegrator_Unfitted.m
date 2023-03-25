classdef RHSintegrator_Unfitted < handle

    properties (Access = private)
        globalConnec
        mesh
        integrators
    end

    methods (Access = public)

        function obj = RHSintegrator_Unfitted(cParams)
            obj.init(cParams);
        end

        function int = integrateInDomain(obj,F)
            s.mesh = obj.mesh.backgroundMesh;
            s.fValues = F;
            p1f  = P1Function(s);

            connecGlobalInner = obj.mesh.innerMesh.globalConnec;
            innerDofs = unique(connecGlobalInner);

            s.mesh = obj.mesh.innerMesh.mesh;
            s.fValues = F(innerDofs);
            p1finner  = P1Function(s);


            a.mesh = obj.mesh.innerMesh.mesh;
            a.type = 'ShapeFunctionFun';
            rhss = RHSintegrator.create(a);
            p1innerInt = rhss.compute(p1finner);


%             b.mesh = obj.mesh.innerCutMesh.mesh;
%             s.cellContainingSubcell = obj.mesh.innerCutMesh.cellContainingSubcell;
%             b.type = 'CutMeshFun';
%             rhss = RHSintegrator.create(b);
%             p1cutInt = rhss.compute(p1finner);

            obj.computeInteriorIntegrators();
            int = obj.integrators.integrateAndSum(F);

            s.fValues = int;
            s.mesh = obj.mesh.backgroundMesh;
            p1int = P1Function(s);


        end

        function int = integrateInBoundary(obj,F)
            obj.computeBoundaryIntegrators();
            int = obj.integrators.integrateAndSum(F);
            
            s.fValues = int;
            s.mesh = obj.mesh.backgroundMesh;
            p1int = P1Function(s);

            s.fValues = F;
            p1f  = P1Function(s);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh = cParams.mesh;
        end

        function computeInteriorIntegrators(obj)
            s = obj.createInteriorParams(obj.mesh,obj.mesh.backgroundMesh.connec);
            obj.integrators = RHSintegrator.create(s);
        end

        function s = createInteriorParams(obj,mesh,connec)
            s.type = 'COMPOSITE';
            s.npnod = mesh.backgroundMesh.nnodes;
            s.compositeParams = cell(0);
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
            s.type         = 'ShapeFunction';
            s.mesh         = innerMesh.mesh;
            s.globalConnec = innerMesh.globalConnec;
            s.npnod        = obj.mesh.backgroundMesh.nnodes;
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

        function computeBoundaryIntegrators(obj)
            uMesh  = obj.mesh;
            s.type = 'COMPOSITE';
            s.npnod = uMesh.backgroundMesh.nnodes;
            s.compositeParams = obj.computeBoundaryParams();
            obj.integrators = RHSintegrator.create(s);
        end

        function s = computeBoundaryParams(obj)
            gConnec = obj.mesh.backgroundMesh.connec;
            s{1} = obj.computeBoundaryCutParams(gConnec);
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

        function s = computeBoundaryCutParams(obj,gConnec)
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