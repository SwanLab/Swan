classdef LHSIntegratorLaplacian < handle

    properties (Access = private)
        mesh
        test, trial
        quadrature
    end

    methods (Access = public)

        function obj = LHSIntegratorLaplacian(cParams)
            obj.initLaplacian(cParams);
            obj.createQuadrature();
        end

        function LHS = compute(obj)
            lhs = obj.computeElementalLHS();
            LHS = obj.assembleMatrix(lhs);
        end

    end

    methods (Access = protected)

        function LHSe = computeElementalLHS(obj)
            xV     = obj.quadrature.posgp;
            dNdxTs = obj.test.evaluateCartesianDerivatives(xV);
            dNdxTr = obj.trial.evaluateCartesianDerivatives(xV);
            dNTs   = obj.reshapeGradient(dNdxTs);
            dNTr   = obj.reshapeGradient(dNdxTr);
            dV     = obj.mesh.computeDvolume(obj.quadrature);

            nNodETs = size(dNdxTs,2);
            nDofETs = nNodETs*obj.test.ndimf;
            nNodETr = size(dNdxTr,2);
            nDofETr = nNodETr*obj.trial.ndimf;
            nElem   = size(dV,2);

            lhs = zeros(nDofETs, nDofETr, nElem);
            for i = 1:nDofETs
                for j = 1:nDofETr
                    dNi   = squeezeParticular(dNTs(:,:,i,:,:),3);
                    dNj   = squeezeParticular(dNTr(:,:,j,:,:),3);
                    dNiNj = pagetensorprod(dNi,dNj,[1 2],[1 2],2,2);
                    dK    = dNiNj.*dV;
                    lhs(i,j,:) = squeeze(lhs(i,j,:))' + sum(dK,1);
                end
            end
            LHSe = lhs;
        end

    end

    methods (Access = private)

        function initLaplacian(obj, cParams)
            obj.mesh  = cParams.mesh;
            obj.test  = cParams.test;
            obj.trial = cParams.trial;
        end

        function createQuadrature(obj)
            orderTr = obj.trial.getOrderNum();
            orderTe = obj.test.getOrderNum();
            order = orderTr + orderTe;
            quad = Quadrature.create(obj.mesh,order);
            obj.quadrature = quad;
        end

        function dN = reshapeGradient(obj,dNdx)
            ndim   = obj.mesh.ndim;
            ndofE  = obj.test.nDofsElem;
            nGauss = obj.quadrature.ngaus;
            nElem  = obj.mesh.nelem;
            dN  = zeros(ndim,ndim,ndofE,nGauss,nElem);
            for i=1:ndim
                dN(i,:,i:ndim:end,:,:) = dNdx;
            end
        end

        function LHS = assembleMatrix(obj, lhs)
            s.fun    = []; % !!!
            assembler = AssemblerFun(s);
            LHS = assembler.assemble(lhs, obj.test, obj.trial);
        end

    end

end