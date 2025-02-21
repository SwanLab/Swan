classdef LHSintegrator_NonLinearNS < LHSintegrator

    methods (Access = public)

        function obj = LHSintegrator_NonLinearNS(cParams)
            obj@LHSintegrator(cParams)
        end

        function LHS = compute(obj)
            lhs = obj.computeElementalLHS();
            LHS = obj.assembleMatrix(lhs);
        end

    end

    methods (Access = private)

        function lhs = computeElementalLHS(obj)
            xV        = obj.quadrature.posgp;
            dNdxTest  = obj.test.evaluateCartesianDerivatives(xV);
            dNdxTrial = obj.trial.evaluateCartesianDerivatives(xV);
            dVolu     = obj.mesh.computeDvolume(obj.quadrature);
        
            nGaus      = obj.quadrature.ngaus;
            nElem      = size(dVolu,2);
            nNodeTest  = size(dNdxTest,2);
            nDofTest   = nNodeTest*obj.test.ndimf;
            nNodeTrial = size(dNdxTrial,2);
            nDofTrial  = nNodeTrial*obj.trial.ndimf;
        
            C = zeros(nDofTest,nDofTrial,nElem);
            
            uValues = obj.test.evaluate(xV);
        
            for igauss = 1 : nGaus
                for inode = 1 : nNodeTest
                    for jnode = 1 : nNodeTrial
                        for iDimf = 1 : obj.test.ndimf
                            idof = obj.test.ndimf*(inode-1) + iDimf;
                            jdof = obj.trial.ndimf*(jnode-1) + iDimf;
                            dvol = dVolu(igauss,:);
                            dNi  = dNdxTest(:,inode,igauss,:);
                            dNi  = permute(dNi, [2 1 3 4]);
                            dNj  = dNdxTrial(:,jnode,igauss,:);
                            
                            % Compute the convective (nonlinear) term
                            % The convective term is (u * grad(u)) or (u_i * (grad u)_i)
                            uGrad = squeeze(uValues(iDimf, igauss,:)); % Get velocity field at gauss point
                            convTerm = uGrad .* squeeze(pagemtimes(dNi, dNj));  % Compute the convective term
        
                            % Now, add this convective term into the matrix K
                            C(idof, jdof, :) = squeeze(C(idof,jdof,:)) + convTerm(:) .* dvol';
                        end
                    end
                end
            end
            lhs = C;
        end
    end
end
