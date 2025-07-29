classdef LHSIntegratorStiffnessElastic < LHSIntegrator

    properties (Access = private)
        material
    end

    methods (Access = public)
        function obj = LHSIntegratorStiffnessElastic(cParams)
            obj@LHSIntegrator(cParams)
            obj.material = cParams.material;
        end

        function LHS = compute(obj)
            lhs = obj.computeElementalLHS();
            LHS = obj.assembleMatrix(lhs);
        end
    end

    methods (Access = protected)
        function lhs = computeElementalLHS(obj)
            xV     = obj.quadrature.posgp;
            dSymN  = ShapeDerSym(obj.test);
            symN   = dSymN.evaluate(xV);
            % C = obj.material.obtainTensor
            C      = obj.material.evaluate(xV);
            nnodeE = obj.mesh.nnodeElem;
            ndim   = obj.mesh.ndim;
            ndofE  = nnodeE*ndim;
            nElem  = obj.mesh.nelem;
            lhs    = zeros(ndofE,ndofE,nElem);
            dV     = obj.mesh.computeDvolume(obj.quadrature);
            for i = 1:ndofE
                for j = 1:ndofE
                    symTrial   = squeezeParticular(symN(:,:,i,:,:),3);
                    symTest    = squeezeParticular(symN(:,:,j,:,:),3);
                    sigN       = pagetensorprod(C,symTrial,[3 4],[1 2],4,2);
                    de         = pagetensorprod(symTest,sigN,[1 2],[1 2],2,2);
                    dK         = de.*dV;
                    lhs(i,j,:) = squeeze(lhs(i,j,:))' + sum(dK,1);
                end
            end
        end
    end
end