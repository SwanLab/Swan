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

            % nDim = obj.mesh.ndim;
            % nVoigt = nDim * (nDim + 1) / 2;
            % j = reshape((0:nNode-1) * obj.test.ndimf, 1, nNode, 1);
            % B = zeros(nVoigt, obj.test.ndimf * nNode, nElem);
            % 

            nDim = obj.mesh.ndim;
            nVoigt = nDim * (nDim + 1) / 2;
            j = reshape((0:nNode-1) * obj.test.ndimf, 1, nNode, 1);
            B = zeros(nVoigt, obj.test.ndimf * nNode, nElem, nGaus);

            d = permute(deriv(1:nDim, :, :, :), [1, 2, 4, 3]);

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

            for iGaus = 1:nGaus           
                BG = B(:,:,:,iGaus);
                C = squeeze(Cmat(:, :, iGaus, :));
                dV = permute(dVolu(iGaus, :), [3, 1, 2]);
                BtCB = pagemtimes(pagemtimes(permute(BG, [2, 1, 3]), C), BG);
                t = BtCB .* dV;
                lhs = lhs + t;
            end



            %lhs2 = lhs;



% C = permute(Cmat, [1, 2, 4, 3]);
% dV = reshape(dVolu, [1, 1, nElem, nGaus]);  % Ensure alignment with BtCB
% 
% BtCB = pagemtimes(pagemtimes(permute(B, [2, 1, 3, 4]), C), B);
% 
% lhs = sum(BtCB .* dV,4);



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