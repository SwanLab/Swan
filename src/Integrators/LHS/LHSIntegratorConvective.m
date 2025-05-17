classdef LHSIntegratorConvective < LHSIntegrator

    methods (Access = public)

        function obj = LHSIntegratorConvective(cParams)
            obj@LHSIntegrator(cParams)
        end

        function LHS = compute(obj, velocityField)
            lhs = obj.computeElementalLHS(velocityField);
            LHS = obj.assembleMatrix(lhs);
        end

    end

    methods (Access = private)

        function lhs = computeElementalLHS(obj,velocityField)
            xV          = obj.quadrature.posgp;  
            NTest       = obj.test.computeShapeFunctions(xV);
            dNdxTest    = obj.test.evaluateCartesianDerivatives(xV);
            dNdxTrial   = obj.trial.evaluateCartesianDerivatives(xV);
            dVolu       = obj.mesh.computeDvolume(obj.quadrature);            

            nGaus       = obj.quadrature.ngaus;
            nElem       = size(dVolu, 2);

            nNodeTest   = size(dNdxTest, 2);
            nNodeTrial  = size(dNdxTrial, 2);
            nDofTest    = nNodeTest * obj.test.ndimf;
            nDofTrial   = nNodeTrial * obj.trial.ndimf;

            lhs = zeros(nDofTest, nDofTrial, nElem);
            uValues = velocityField.evaluate(xV);

            
            for iGaus = 1:nGaus
                uGp = squeeze(uValues(:, iGaus, :));  % [nDim x nElem]
                dV   = dVolu(iGaus, :);                % [1 x nElem]
            
                for iNode = 1:nNodeTest
                    phiI = NTest(iNode, iGaus);        % scalar
            
                    for jNode = 1:nNodeTrial
                        gradPhiJ = squeeze(dNdxTrial(:, jNode, iGaus, :)); % [nDim x nElem]
            
                        for iDim = 1:obj.test.ndimf
                            idof = obj.test.ndimf * (iNode - 1) + iDim;
            
                            for jDim = 1:obj.trial.ndimf
                                jdof = obj.trial.ndimf * (jNode - 1) + jDim;
                                convTerm = uGp(jDim, :) .* gradPhiJ(jDim, :);  % [1 x nElem]
            
                                % update matrix entry
                                lhs(idof, jdof, :) = lhs(idof, jdof, :) + ...
                                    reshape(phiI * convTerm .* dV, [1 1 nElem]);
                            end
                        end
                    end
                end
            end
        end

    end
end