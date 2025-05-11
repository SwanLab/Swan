classdef LHSintegrator_NonLinearNS < LHSintegrator

    methods (Access = public)

        function obj = LHSintegrator_NonLinearNS(cParams)
            obj@LHSintegrator(cParams)
        end

        function LHS = compute(obj, velocityField)
            lhs = obj.computeElementalLHS(velocityField);
            LHS = obj.assembleMatrix(lhs);
        end

    end

    methods (Access = private)

        function lhs = computeElementalLHS(obj, velocityField)
            xV          = obj.quadrature.posgp;             
            dNdxTest    = obj.test.evaluateCartesianDerivatives(xV);
            dNdxTrial   = obj.trial.evaluateCartesianDerivatives(xV);
            dVolu       = obj.mesh.computeDvolume(obj.quadrature);

            NTest  = obj.test.computeShapeFunctions(xV);
            NTrial = obj.trial.computeShapeFunctions(xV);

            nGaus       = obj.quadrature.ngaus;
            nElem       = size(dVolu, 2);
            nNodeTest   = size(dNdxTest, 2);
            nNodeTrial  = size(dNdxTrial, 2);
            nDim        = size(dNdxTrial,1);
            nDofTest    = nNodeTest * obj.test.ndimf;
            nDofTrial   = nNodeTrial * obj.trial.ndimf;

            C = zeros(nDofTest, nDofTrial, nElem);
            uValues = velocityField.evaluate(xV);

            
            for igauss = 1:nGaus
                u_gp = squeeze(uValues(:, igauss, :));  % [nDim x nElem]
                dV   = dVolu(igauss, :);                % [1 x nElem]
            
                for idof = 1:nNodeTest
                    phi_i = NTest(idof, igauss);        % scalar
            
                    for jdof = 1:nNodeTrial
                        grad_phi_j = squeeze(dNdxTrial(:, jdof, igauss, :)); % [nDim x nElem]
            
                        for alpha = 1:nDim
                            idof2 = nDim * (idof - 1) + alpha;
            
                            for beta = 1:nDim
                                jdof2 = nDim * (jdof - 1) + beta;
                                conv_term = u_gp(beta, :) .* grad_phi_j(beta, :);  % [1 x nElem]
            
                                % update matrix entry
                                C(idof2, jdof2, :) = C(idof2, jdof2, :) + ...
                                    reshape(phi_i * conv_term .* dV, [1 1 nElem]);
                            end
                        end
                    end
                end
            end


            lhs = C;
        end

    end
end
