classdef LHSintegratorMassDifferentMeshes < handle

    properties (Access = private)
        uMesh
        testOrder, trialOrder
        testI, testIC, trial
        quadrature
        quadratureOrder
    end

    methods (Access = public)

        function obj = LHSintegratorMassDifferentMeshes(cParams)
            obj.init(cParams);
            obj.createTestTrial();
            obj.createQuadrature();
        end

        function LHS = compute(obj)
            lhsI  = obj.computeElementalLHSInner();
            lhsIC = obj.computeElementalLHSInnerCut();
            LHS   = lhsI+lhsIC;
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.uMesh           = cParams.uMesh;
            obj.testOrder       = cParams.testOrder;
            obj.trialOrder      = cParams.trialOrder;
            obj.quadratureOrder = cParams.quadratureOrder;
        end

        function createQuadrature(obj)
            m = obj.uMesh.innerCutMesh.mesh;
            q = Quadrature.create(m,obj.quadratureOrder);
            obj.quadrature = q;
        end

        function createTestTrial(obj)
            obj.trial  = LagrangianFunction.create(obj.uMesh.backgroundMesh, 1, obj.testOrder);
            obj.testI  = LagrangianFunction.create(obj.uMesh.innerMesh.mesh, 1, obj.trialOrder);
            obj.testIC = LagrangianFunction.create(obj.uMesh.innerCutMesh.mesh, 1, obj.trialOrder);
        end

        function lhs = assembleMatrix(obj, M)
            nElem         = obj.uMesh.innerMesh.mesh.nelem;
            iCMesh        = obj.uMesh.innerCutMesh;
            subDomain     = obj.uMesh.createInnerMesh();
            globalConnec  = subDomain.connec(nElem+1:end,:);
            nDofTest      = subDomain.nnodes;
            nDofTrial     = obj.trial.nDofs;
            connecTest    = obj.testIC.computeDofConnectivity()';
            connecTrial   = obj.trial.computeDofConnectivity()';
            nDofElemTest  = size(connecTest,2);
            nDofElemTrial = size(connecTrial,2);
            l2g(iCMesh.mesh.connec(:)) = globalConnec(:);
            lhs = sparse(nDofTest,nDofTrial);
            for iDof = 1:nDofElemTest
                for jDof = 1:nDofElemTrial
                    conTestReal = zeros(size(connecTrial,1),1);
                    int      = M(iDof,jDof,:);
                    conTest  = connecTest(:,iDof);
                    conTestReal(l2g) = conTest;
                    conTrial = connecTrial(:,jDof);
                    lhs   = lhs + sparse(conTestReal,conTrial,squeeze(int),nDofTest,nDofTrial);
                end
            end
        end

        function lhs = computeElementalLHSInner(obj)
            subDomain = obj.uMesh.createInnerMesh();
            nDofTest  = subDomain.nnodes;
            nDofTrial = obj.trial.nDofs;
            lhs       = sparse(nDofTest,nDofTrial);
            iMesh     = obj.uMesh.innerMesh.mesh;
            setNodes  = 1:iMesh.nnodes;
            l2g(iMesh.connec(:)) = obj.uMesh.innerMesh.globalConnec(:);
            s.mesh    = iMesh;
            s.test    = obj.testI;
            s.trial   = LagrangianFunction.create(obj.uMesh.innerMesh.mesh, 1, obj.trialOrder);
            s.quadratureOrder = 2;
            s.type     = 'MassMatrix';
            lhsLoc     = LHSintegrator.create(s);
            MLoc       = lhsLoc.compute();
            lhs(setNodes,l2g) = MLoc;
        end

        function lhs = computeElementalLHSInnerCut(obj)
            iCMesh    = obj.uMesh.innerCutMesh;
            globCell  = iCMesh.cellContainingSubcell;
            quad      = obj.quadrature;
            xVLoc     = quad.posgp;
            dVolu     = iCMesh.mesh.computeDvolume(quad);
            isoMesh   = obj.obtainIsoparametricMesh(iCMesh);
            xV        = isoMesh.evaluate(xVLoc);
            Ni        = obj.testIC.computeShapeFunctions(xVLoc);
            Nj        = obj.trial.computeShapeFunctions(xV);
            nDofTestLoc = size(Ni,1);
            nDofTrialLoc = size(Nj,1);
            nGaus     = quad.ngaus;
            nElem     = obj.uMesh.backgroundMesh.nelem;
            M         = zeros(nDofTestLoc,nDofTrialLoc,nElem);
            for igauss = 1 :nGaus
                for idof= 1:nDofTestLoc
                    for jdof= 1:nDofTrialLoc
                        dvol = dVolu(igauss,:);
                        Nig  = Ni(idof,igauss,:);
                        Njg  = Nj(jdof,igauss,:);
                        v    = squeeze(Nig.*Njg);
                        MLoc = v(:).*dvol';
                        int  = accumarray(globCell,MLoc,[nElem,1],@sum,0);
                        M(idof,jdof,:) = M(idof,jdof,:) + reshape(int,[1,1,nElem]);
                    end
                end
            end
            lhs = obj.assembleMatrix(M);
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
    end

end