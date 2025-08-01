classdef LHSIntegratorAnisotropicStiffness < LHSIntegrator

    properties (Access = private)
        A
    end

    methods (Access = public)

        function obj = LHSIntegratorAnisotropicStiffness(cParams)
            obj@LHSIntegrator(cParams)
            obj.A = cParams.A;
        end

        function LHS = compute(obj)
            lhs = obj.computeElementalLHS();
            LHS = obj.assembleMatrix(lhs);
        end

    end

    methods (Access = protected)

        function lhs = computeElementalLHS(obj)
            xV     = obj.quadrature.posgp;
            dN   = ShapeDer(obj.test).evaluate(xV);
            C      = obj.A.evaluate(xV);
            nnodeE = obj.mesh.nnodeElem;
            ndim   = obj.trial.ndimf;
            ndofE  = nnodeE*ndim;
            nElem  = obj.mesh.nelem;
            lhs    = zeros(ndofE,ndofE,nElem);
            dV     = obj.mesh.computeDvolume(obj.quadrature);
            for i = 1:ndofE
                for j = 1:ndofE
                    dTrial     = squeezeParticular(dN(:,:,i,:,:),3);
                    dTest      = squeezeParticular(dN(:,:,j,:,:),3);
                    sigN       = pagetensorprod(C,dTrial,[2],[2],2,2);
                    de         = pagetensorprod(dTest,sigN,[2],[1],2,2);
                    dK         = squeeze(de).*dV;
                    lhs(i,j,:) = squeeze(lhs(i,j,:))' + sum(dK,1);
                end
            end
        end

    end
end