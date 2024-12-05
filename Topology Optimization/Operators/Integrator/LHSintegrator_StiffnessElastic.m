classdef LHSintegrator_StiffnessElastic < LHSintegrator

    properties (Access = private)
        material
    end

    methods (Access = public)

        function obj = LHSintegrator_StiffnessElastic(cParams)
            obj@LHSintegrator(cParams)
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
            B    = obj.computeB(dNdx);            
            dV   = obj.mesh.computeDvolume(obj.quadrature);
            dV   = permute(dV,[3, 4, 2, 1]);
            Cmat = obj.material.evaluate(xV);
            Cmat = permute(Cmat,[1 2 4 3]);
            Bt   = permute(B,[2 1 3 4]);
            BtC  = pagemtimes(Bt, Cmat);
            BtCB = pagemtimes(BtC, B);
            lhs  = sum(BtCB .* dV, 4);
        end

    end

    methods (Access = private)

        function B = computeB(obj, dNdx)
            nGaus = obj.quadrature.ngaus;
            nNode = size(dNdx, 2);
            nElem = obj.mesh.nelem;
            nDim  = obj.mesh.ndim;
            nDimF = obj.test.ndimf;
            nDofE = nNode*nDimF;            
            nVoigt = nDim * (nDim + 1) / 2;
            j = nDimF * reshape((1:nNode) - 1, 1, nNode, 1);
            B = zeros(nVoigt, nDofE, nElem, nGaus);
            d = permute(dNdx, [1, 2, 4, 3]);
            if nDim == 2
                B(1, j + 1, :, :) = d(1, :, :, :);
                B(2, j + 2, :, :) = d(2, :, :, :);
                B(3, j + 1, :, :) = d(2, :, :, :);
                B(3, j + 2, :, :) = d(1, :, :, :);
            elseif nDim == 3
                B(1, j + 1, :, :) = d(1, :, :, :);
                B(2, j + 2, :, :) = d(2, :, :, :);
                B(3, j + 3, :, :) = d(3, :, :, :);
                B(4, j + 1, :, :) = d(2, :, :, :);
                B(4, j + 2, :, :) = d(1, :, :, :);
                B(5, j + 1, :, :) = d(3, :, :, :);
                B(5, j + 3, :, :) = d(1, :, :, :);
                B(6, j + 2, :, :) = d(3, :, :, :);
                B(6, j + 3, :, :) = d(2, :, :, :);
            end
        end

    end

end