classdef LHSIntegratorStiffness < LHSIntegrator

    methods (Access = public)

        function obj = LHSIntegratorStiffness(cParams)
            obj@LHSIntegrator(cParams)
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

            K = zeros(nDofTest,nDofTrial,nElem);
            for igauss = 1 :nGaus
                for inode= 1:nNodeTest
                    for jnode= 1:nNodeTrial
                        for iDimf = 1:obj.test.ndimf
                            idof = obj.test.ndimf*(inode-1)+iDimf;
                            jdof = obj.trial.ndimf*(jnode-1)+iDimf;
                            dvol = dVolu(igauss,:);
                            dNi  = dNdxTest(:,inode,igauss,:);
                            dNi  = permute(dNi,[2 1 3 4]);
                            dNj  = dNdxTrial(:,jnode,igauss,:);
                            v    = squeeze(pagemtimes(dNi,dNj));
                            K(idof, jdof, :)= squeeze(K(idof,jdof,:)) ...
                                + v(:).*dvol';
                        end
                    end
                end
            end
            lhs = K;
        end
    end
end