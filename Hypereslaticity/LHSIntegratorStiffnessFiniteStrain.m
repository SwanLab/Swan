classdef LHSIntegratorStiffnessFiniteStrain < LHSIntegrator


    methods (Access = public)
        function obj = LHSIntegratorStiffnessFiniteStrain(cParams)
            obj@LHSIntegrator(cParams)
        end

        function LHS = compute(obj,mat)
            lhs = obj.computeElementalLHS(mat);
            LHS = obj.assembleMatrix(lhs);
        end
    end

    methods (Access = protected)
        function lhs = computeElementalLHS(obj,mat)
            xV     = obj.quadrature.posgp;
            dN     = ShapeDerTensor(obj.test);
            dN     = dN.evaluate(xV);
            C      = mat.evaluate(xV);
            nnodeE = obj.mesh.nnodeElem;
            ndim   = obj.mesh.ndim;
            ndofE  = nnodeE*ndim;
            nElem  = obj.mesh.nelem;
            lhs    = zeros(ndofE,ndofE,nElem);
            dV     = obj.mesh.computeDvolume(obj.quadrature);
            for i = 1:ndofE
                for j = 1:ndofE
                    dTrial     = squeezeParticular(dN(:,:,i,:,:),3);
                    dTest      = squeezeParticular(dN(:,:,j,:,:),3);
                    sigN       = pagetensorprod(C,dTrial,[3 4],[1 2],4,2);
                    de         = pagetensorprod(dTest,sigN,[1 2],[1 2],2,2);
                    dK         = de.*dV;
                    lhs(i,j,:) = squeeze(lhs(i,j,:))' + sum(dK,1);
                end
            end
        end
    end
end