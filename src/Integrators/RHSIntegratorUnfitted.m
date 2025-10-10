classdef RHSIntegratorUnfitted < handle

    properties (Access = private)
        quadType
        unfittedMesh
    end

    methods (Access = public)
        function obj = RHSIntegratorUnfitted(cParams)
            obj.init(cParams);
        end

        function int = computeUnfitted(obj,f,test)
            intInner    = obj.integrateInnerMeshFunction(@(v) f(v).innerMeshFunction,test);
            iCMesh      = obj.unfittedMesh.innerCutMesh;
            intInnerCut = obj.integrateCutMeshFunction(@(v) f(v).innerCutMeshFunction,test,iCMesh);
            int         = intInner+intInnerCut;
        end

        function int = computeUnfittedBoundary(obj,f,test)
            bCMesh      = obj.unfittedMesh.boundaryCutMesh;
            intBoundCut = obj.integrateCutMeshFunction(@(v) f(v).boundaryCutMeshFunction,test,bCMesh);
            intUnfBound = obj.integrateUnfittedBoundaryMeshFunction(@(v) f(v).unfittedBoundaryMeshFunction,test);
            int         = intBoundCut+intUnfBound;
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.quadType     = cParams.quadType;
            obj.unfittedMesh = cParams.mesh;
        end

        function int = integrateInnerMeshFunction(obj,fI,test)
            dofs  = size(test.fValues,1);
            int   = zeros(dofs,1);
            iMesh = obj.unfittedMesh.innerMesh;
            if ~isempty(iMesh)
                fullCells = iMesh.fullCells;
                testLoc   = LagrangianFunction.create(iMesh.mesh,test.ndimf,test.order);
                intLoc    = IntegrateRHS(fI,testLoc,iMesh.mesh,'Domain',obj.quadType);
                dofG      = test.getDofConnec();
                dofL      = testLoc.getDofConnec();
                l2g(dofL) = dofG(fullCells,:);
                int(l2g)  = intLoc;
            end
        end

        function int = integrateCutMeshFunction(obj,f,test,cutMesh)
            if ~isempty(cutMesh)
                int = IntegrateRHSCutMesh(f,test,cutMesh,obj.quadType);
            else
                dofs = size(test.fValues,1);
                int  = zeros(dofs,1);
            end
        end

        function int = integrateUnfittedBoundaryMeshFunction(obj,f,test)
            dofs = size(test.fValues,1);
            int  = zeros(dofs,1);
            if ~isempty(obj.unfittedMesh.unfittedBoundaryMesh)
                uMeshBound = obj.unfittedMesh.unfittedBoundaryMesh.getActiveMesh();
                nFun       = length(uMeshBound);
                for i = 1:nFun
                    uFi          = @(v) subsref(f(v), struct('type','{}','subs',{{i}}));
                    intLoci      = obj.integrateExternalBoundaryLocal(uFi,test,uMeshBound{i});
                    conG         = obj.unfittedMesh.unfittedBoundaryMesh.getGlobalConnec{i};
                    conL         = uMeshBound{i}.backgroundMesh.connec;
                    l2g(conL(:)) = conG(:);
                    int(l2g)     = int(l2g)+ intLoci;
                    l2g          = [];
                end 
            end
        end

        function intLoc = integrateExternalBoundaryLocal(obj,uFi,test,uMeshi)
            mi      = uMeshi.backgroundMesh;
            testLoc = LagrangianFunction.create(mi,test.ndimf,test.order);
            intLoc  = IntegrateRHS(uFi,testLoc,mi,'Domain',obj.quadType);
        end
    end
end