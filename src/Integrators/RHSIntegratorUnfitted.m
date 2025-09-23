classdef RHSIntegratorUnfitted < handle

    properties (Access = private)
        quadType
        unfittedMesh
    end

    properties (Access = private)
        innerIntegrator
        innerCutQuad
        boundCutQuad
   end

    methods (Access = public)
        function obj = RHSIntegratorUnfitted(cParams)
            obj.init(cParams);
            obj.createQuadratureInnerCut();
            obj.createQuadratureBoundaryCut();
        end

        function int = compute(obj,uFun,test)
            switch class(uFun)
                case 'UnfittedFunction'
                    intInner    = obj.integrateInnerMeshFunction(uFun,test);
                    fInnerCut   = uFun.innerCutMeshFunction;
                    iCMesh      = obj.unfittedMesh.innerCutMesh;
                    qICMesh     = obj.innerCutQuad;
                    intInnerCut = obj.integrateCutMeshFunction(fInnerCut,test,iCMesh,qICMesh);
                    int         = intInner+intInnerCut;
                    %s.mesh=obj.unfittedMesh.backgroundMesh;s.order='P1';s.fValues=intInner;fP1=LagrangianFunction(s);fP1.plot;s.fValues=intInnerCut;fP1=LagrangianFunction(s);fP1.plot;s.fValues=int;fP1=LagrangianFunction(s);fP1.plot
                case 'UnfittedBoundaryFunction'
                    fBoundCut   = uFun.boundaryCutMeshFunction;
                    bCMesh      = obj.unfittedMesh.boundaryCutMesh;
                    qBCMesh     = obj.boundCutQuad;
                    intBoundCut = obj.integrateCutMeshFunction(fBoundCut,test,bCMesh,qBCMesh);
                    intUnfBound = obj.integrateUnfittedBoundaryMeshFunction(uFun,test);
                    int         = intBoundCut+intUnfBound;
            end
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.quadType     = cParams.quadType;
            obj.unfittedMesh = cParams.mesh;
        end

        function createQuadratureInnerCut(obj)
            if ~isempty(obj.unfittedMesh.innerCutMesh)
                m = obj.unfittedMesh.innerCutMesh.mesh;
                q = Quadrature.create(m,obj.quadType);
                obj.innerCutQuad = q;
            end
        end

        function createQuadratureBoundaryCut(obj)
            if ~isempty(obj.unfittedMesh.boundaryCutMesh)
                m = obj.unfittedMesh.boundaryCutMesh.mesh;
                q = Quadrature.create(m,obj.quadType);
                obj.boundCutQuad = q;
            end
        end

        function integrator = createIntegratorUnfittedBoundary(obj,uMesh)
            s.mesh     = uMesh;
            s.quadType = obj.quadType;
            integrator = RHSIntegratorUnfitted(s);
        end

        function int = integrateInnerMeshFunction(obj,uFun,test)
            dofs  = size(test.fValues,1);
            int   = zeros(dofs,1);
            iMesh = obj.unfittedMesh.innerMesh;
            if ~isempty(iMesh)
                fullCells = iMesh.fullCells;
                fInner    = uFun.innerMeshFunction;
                testLoc   = LagrangianFunction.create(iMesh.mesh,test.ndimf,test.order);
                intLoc    = IntegrateRHS(@(v) DP(v,fInner),testLoc,iMesh.mesh,obj.quadType);
                dofG      = test.getDofConnec();
                dofL      = testLoc.getDofConnec();
                l2g(dofL) = dofG(fullCells,:);
                int(l2g)  = intLoc;
            end
        end

        function int = integrateCutMeshFunction(obj,f,test,cutMesh,quad)
            if ~isempty(cutMesh)
                xVLoc    = quad.posgp;
                fG       = f.evaluate(xVLoc);
                dV       = cutMesh.mesh.computeDvolume(quad);
                isoMesh  = obj.obtainIsoparametricMesh(cutMesh);
                xV       = isoMesh.evaluate(xVLoc);
                globCell = cutMesh.cellContainingSubcell;
                N        = test.computeShapeFunctions(xV);
                nDofElem = size(N,1);
                nElem    = test.mesh.nelem;
                nGaus    = quad.ngaus;
                intElem  = zeros(nElem,nDofElem);
                for iDof = 1:nDofElem
                    for iGaus = 1:nGaus
                        dVg(:,1) = dV(iGaus, :);
                        fV   = squeeze(fG(1,iGaus,:));
                        Ni   = squeeze(N(iDof,iGaus,:));
                        fNdV = Ni.*fV.*dVg;
                        intElem(:,iDof) = intElem(:,iDof) + accumarray(globCell,fNdV,[nElem,1],@sum,0);
                    end
                end
                int = obj.assembleIntegrand(test,intElem);


               % testLoc = UnfittedFunction.create(uMesh,test);
               % f2 = @(v) DomainFunction.create(@(xV) DP(v.evaluate(isoMesh.evaluate(xV)),f),cutMesh.mesh,1);
               % f3 = @(v) DomainFunction.create(@(xV) accumarray(globCell,f2(v).evaluate(xV),[1,nGaus,nElem],@sum,0),cutMesh.mesh,1);
               % intElem2 = IntegrateRHS(f3,testLoc.innerCutMeshFunction,cutMesh.mesh,obj.quadType); % Or create IntegrateRHSCutMesh
            else
                dofs = size(test.fValues,1);
                int  = zeros(dofs,1);
            end
        end

        function int = integrateUnfittedBoundaryMeshFunction(obj,uFun,test)
            dofs = size(test.fValues,1);
            int  = zeros(dofs,1);
            if ~isempty(obj.unfittedMesh.unfittedBoundaryMesh)
                uMeshBound = obj.unfittedMesh.unfittedBoundaryMesh.getActiveMesh();
                nFun       = length(uMeshBound);
                for i = 1:nFun
                    uFi          = obj.computeUnfittedFunctionAtExternalBoundary(uFun,i);
                    intLoci      = obj.integrateExternalBoundaryLocal(uFi,test,uMeshBound{i});
                    conG         = obj.unfittedMesh.unfittedBoundaryMesh.getGlobalConnec{i};
                    conL         = uMeshBound{i}.backgroundMesh.connec;
                    l2g(conL(:)) = conG(:);
                    int(l2g)     = int(l2g)+ intLoci;
                    l2g          = [];
                end 
            end
        end

        function uF = computeUnfittedFunctionAtExternalBoundary(obj,uFun,i)
            uFunBound  = uFun.unfittedBoundaryMeshFunction.activeFuns;
            uMeshBound = obj.unfittedMesh.unfittedBoundaryMesh.getActiveMesh();
            funi       = uFunBound{i}.backgroundFunction;
            s.fun      = funi;
            s.uMesh    = uMeshBound{i};
            uF         = UnfittedFunction(s);
        end

        function intLoc = integrateExternalBoundaryLocal(obj,uFi,test,uMeshi)
            mi         = uMeshi.backgroundMesh;
            integrator = obj.createIntegratorUnfittedBoundary(uMeshi);
            testLoc    = LagrangianFunction.create(mi,test.ndimf,test.order);
            intLoc     = integrator.compute(uFi,testLoc);
        end

    end

    methods (Static, Access = private)

        function m = obtainIsoparametricMesh(cutMesh)
            coord      = cutMesh.xCoordsIso;
            nDim       = size(coord,1);
            nNode      = size(coord,2);
            nElem      = size(coord,3);
            msh.connec = reshape(1:nElem*nNode,nNode,nElem)';
            msh.type   = cutMesh.mesh.type;
            s.fValues  = reshape(coord,nDim,[])';
            s.mesh     = msh;
            s.order    = 'P1';
            m          = LagrangianFunction(s);
        end

        function f = assembleIntegrand(test,rhsElem)
            integrand = rhsElem;
            connec    = test.getDofConnec();
            ndofs     = max(max(connec));
            nDofElem  = size(connec,2);
            f         = zeros(ndofs,1);
            for iDof = 1:nDofElem
                int = integrand(:,iDof);
                con = connec(:,iDof);
                f   = f + accumarray(con,int,[ndofs,1],@sum,0);
            end
        end
    end
end