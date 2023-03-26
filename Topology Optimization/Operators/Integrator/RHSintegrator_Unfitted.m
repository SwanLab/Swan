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

            sAF.fHandle = @(x) x(1,:,:);
            sAF.ndimf   = 1;
            sAF.mesh    = obj.mesh.backgroundMesh;
            xFun = AnalyticalFunction(sAF);
            SAMPLEP1 = xFun.project('P1');

            
            globalConnecInner = obj.mesh.innerMesh.globalConnec;
            localConnecInner  = obj.mesh.innerMesh.mesh.connec;
            innerDofsGlobal = unique(globalConnecInner(:)); % nein
            innerDofs = unique(localConnecInner); % nein
            innerlocal2innerglobal(localConnecInner(:)) = globalConnecInner(:);


%             fV_global = zeros(861,1);
%             fV_global(innerDofsGlobal) = SAMPLEP1.fValues(innerDofsGlobal);
% 
%             fV_local = zeros(733,1);
%             fV_local(innerDofs) = fV_global(innerlocal2innerglobal(innerDofs));


            fV_global = zeros(length(F),1);
            fV_global(innerDofsGlobal) = F(innerDofsGlobal);

            fV_local = zeros(length(innerDofs),1);
            fV_local(innerDofs) = fV_global(innerlocal2innerglobal(innerDofs));

            s.mesh = obj.mesh.innerMesh.mesh;
            s.fValues = fV_local;
            p1finner  = P1Function(s);


            a.mesh = obj.mesh.innerMesh.mesh;
            a.type = 'ShapeFunctionFun';
            rhss = RHSintegrator.create(a);
            p1innerInt = rhss.compute(p1finner);


            b.mesh = obj.mesh.innerCutMesh.mesh;
            b.cellContainingSubcell = obj.mesh.innerCutMesh.cellContainingSubcell;
            b.type = 'CutMeshFun';
            b.xCoordsIso            = obj.mesh.innerCutMesh.xCoordsIso;
            b.globalConnec          = obj.mesh.backgroundMesh.connec;
            b.npnod                 = obj.mesh.backgroundMesh.nnodes;
            b.backgroundMeshType    = obj.mesh.backgroundMesh.type;
            b.backgroundMesh        = obj.mesh.backgroundMesh;
            rhss = RHSintegrator.create(b);
            cutFVals = rhss.compute(p1f.fValues);

            cc.mesh = obj.mesh.backgroundMesh;
            cc.fValues = cutFVals;
            p1cutInt = P1Function(cc);

            fVintInnerinGlobal = zeros(length(F),1);
            fVintInnerinGlobal(innerlocal2innerglobal(innerDofs)) = p1innerInt.fValues(innerDofs);
            fValuesInt = p1cutInt.fValues + fVintInnerinGlobal;

            zzz.mesh = obj.mesh.backgroundMesh;
            zzz.fValues = fValuesInt;
            myIntp1 = P1Function(zzz);

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