classdef RHSIntegratorUnfitted < handle

    properties (Access = private)
        quadType
        unfittedMesh
    end

    properties (Access = private)
        innerIntegrator
        innerCutQuad
   end

    methods (Access = public)
        function obj = RHSIntegratorUnfitted(cParams)
            obj.init(cParams);
            obj.createIntegratorInner();
            obj.createQuadratureInnerCut();
        end

        function int = compute(obj,uMeshFun,test)
            fInner      = uMeshFun.innerMeshFunction;
            fInnerCut   = uMeshFun.innerCutMeshFunction;
            intInner    = obj.integrateInnerMeshFunction(fInner,test);
            intInnerCut = obj.integrateInnerCutMeshFunction(fInnerCut,test);
            int         = intInner+intInnerCut;
            %s.mesh=obj.unfittedMesh.backgroundMesh;s.order='P1';s.fValues=intInner;fP1=LagrangianFunction(s);fP1.plot;s.fValues=intInnerCut;fP1=LagrangianFunction(s);fP1.plot;s.fValues=int;fP1=LagrangianFunction(s);fP1.plot
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.quadType     = cParams.quadType;
            obj.unfittedMesh = cParams.mesh;
        end

        function createIntegratorInner(obj)
            s.mesh     = obj.unfittedMesh.innerMesh.mesh;
            s.type     = 'ShapeFunction';
            s.quadType = obj.quadType;
            int        = RHSintegrator.create(s);
            obj.innerIntegrator = int;
        end

        function createQuadratureInnerCut(obj)
            m = obj.unfittedMesh.innerCutMesh.mesh;
            q = Quadrature.set(m.type);
            q.computeQuadrature(obj.quadType);
            obj.innerCutQuad = q;
        end

        function int = integrateInnerMeshFunction(obj,f,test)
            intLoc        = obj.integrateInnerMeshLocal(f,test);
            iMesh         = obj.unfittedMesh.innerMesh;
            dofs          = size(test.fValues,1);
            int           = zeros(dofs,1);
            conG          = iMesh.globalConnec;
            conBL         = iMesh.mesh.connec;
            l2g(conBL(:)) = conG(:);
            int(l2g)      = intLoc;
        end

        function intLoc = integrateInnerMeshLocal(obj,f,test)
            m       = obj.unfittedMesh.innerMesh.mesh;
            testLoc = LagrangianFunction.create(m,test.ndimf,test.order);
            intLoc  = obj.innerIntegrator.compute(f,testLoc);
        end

        function int = integrateInnerCutMeshFunction(obj,f,test)
            iCMesh   = obj.unfittedMesh.innerCutMesh;
            q        = obj.innerCutQuad;
            xVLoc    = q.posgp;
            fG       = f.evaluate(xVLoc);
            dV       = iCMesh.mesh.computeDvolume(q);
            isoMesh  = obj.obtainIsoparametricMesh();
            xV       = isoMesh.computeXgauss(xVLoc);
            globCell = iCMesh.cellContainingSubcell;
            N        = test.computeShapeFunctions(xV);
            nDofElem = size(N,1);
            nElem    = test.mesh.nelem;
            nGaus    = obj.innerCutQuad.ngaus;
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
        end

        function m = obtainIsoparametricMesh(obj)
            iCMesh   = obj.unfittedMesh.innerCutMesh;
            coord    = iCMesh.xCoordsIso;
            nDim     = size(coord,1);
            nNode    = size(coord,2);
            nElem    = size(coord,3);
            s.coord  = reshape(coord,nDim,[])';
            s.connec = reshape(1:nElem*nNode,nNode,nElem)';
            m        = Mesh.create(s); 
        end
    end

    methods (Static, Access = private)
        function f = assembleIntegrand(test,rhsElem)
            integrand = rhsElem;
            connec    = test.computeDofConnectivity()';
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