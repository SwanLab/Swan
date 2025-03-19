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
            xV   = obj.quadrature.posgp;
            dNdx = obj.test.evaluateCartesianDerivatives(xV);
            symN = obj.computeSymGradShapeFunctions(dNdx);
            C    = obj.material.evaluate(xV);

            nnodeE = obj.mesh.nnodeElem;
            ndim   = obj.mesh.ndim;
            ndofE  = nnodeE*ndim;
            nElem  = obj.mesh.nelem;
           
            lhs = zeros(ndofE,ndofE,nElem);
            dV  = obj.mesh.computeDvolume(obj.quadrature);
            for i=1:ndofE
                for j=1:ndofE
                    symTrial = squeezeParticular(symN(:,:,i,:,:),3);
                    symTest  = squeezeParticular(symN(:,:,j,:,:),3);
                    sigN = pagetensorprod(C,symTrial,[3 4],[1 2],4,2);
                    Kelem = pagetensorprod(symTest,sigN,[1 2],[1 2],2,2);
                    Kelem = Kelem.*dV;
                    lhs(i,j,:) = squeeze(lhs(i,j,:))' + sum(Kelem,1);
                end
            end
        end

    end

    methods (Access = private)

        function symN = computeSymGradShapeFunctions(obj, dNdx)
            nnodeE = obj.mesh.nnodeElem;
            ndim   = obj.mesh.ndim;
            ndofE  = nnodeE*ndim;
            nGauss = size(obj.quadrature.posgp,2);
            nElem  = obj.mesh.nelem;

            gradN = zeros(ndim,ndim,ndofE,nGauss,nElem);
            for i=1:ndim
                gradN(i,:,i:ndim:end,:,:) = dNdx;
            end
            symN = 0.5*(gradN + pagetranspose(gradN));
        end

    end

end