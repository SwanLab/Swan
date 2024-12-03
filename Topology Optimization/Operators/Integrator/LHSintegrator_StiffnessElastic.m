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
            xV = obj.quadrature.posgp;
            dNdx = obj.test.evaluateCartesianDerivatives(xV);
            dVolu = obj.mesh.computeDvolume(obj.quadrature);
            nGaus = obj.quadrature.ngaus;
            nElem = obj.mesh.nelem;
            nNodE = size(dNdx, 2);
            nDofE = nNodE * obj.test.ndimf;
            lhs = zeros(nDofE, nDofE, nElem);
        %    Bcomp = obj.createBComputer(dNdx);
            Cmat = obj.material.evaluate(xV);
            deriv = dNdx;
            nNode = size(deriv, 2);
            nElem = size(deriv, 4);
            nVoigt = 6;
            for iGaus = 1:nGaus
                j = reshape((0:nNode-1) * obj.test.ndimf, 1, nNode, 1);
                d1 = permute(deriv(1, :, iGaus, :), [1, 2, 4, 3]);
                d2 = permute(deriv(2, :, iGaus, :), [1, 2, 4, 3]);
                d3 = permute(deriv(3, :, iGaus, :), [1, 2, 4, 3]);
                B = zeros(nVoigt, obj.test.ndimf * nNode, nElem);
                B(1, j + 1, :) = d1;
                B(2, j + 2, :) = d2;
                B(3, j + 3, :) = d3;
                B(4, j + 1, :) = d2;
                B(4, j + 2, :) = d1;
                B(5, j + 1, :) = d3;
                B(5, j + 3, :) = d1;
                B(6, j + 2, :) = d3;
                B(6, j + 3, :) = d2;


                C = squeeze(Cmat(:, :, iGaus, :));
                dV = permute(dVolu(iGaus, :), [3, 1, 2]);
                BtCB = pagemtimes(pagemtimes(permute(B, [2, 1, 3]), C), B);
                lhs = lhs + BtCB .* dV;
            end

        end

    end

    methods (Access = private)

        function Bcomp = createBComputer(obj, dNdx)
            s.fun  = obj.test;
            s.dNdx = dNdx;
            Bcomp = BMatrixComputer(s);
        end

    end

end