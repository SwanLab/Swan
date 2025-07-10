classdef LHSIntegratorAnisotropicStiffness < LHSIntegrator

    properties (Access = private)
        A
    end

    methods (Access = public)

        function obj = LHSIntegratorAnisotropicStiffness(cParams)
            obj@LHSIntegrator(cParams)
            obj.A = cParams.A;
        end

        function LHS = compute(obj)
            lhs = obj.computeElementalLHS();
            LHS = obj.assembleMatrix(lhs);
        end

    end

    methods (Access = protected)

        function lhs = computeElementalLHS(obj)
            xV = obj.quadrature.posgp;
            dNTs = ShapeDer(obj.test);
            dNTs = dNTs.evaluate(xV);
            dNTr = ShapeDer(obj.trial);
            dNTr = dNTr.evaluate(xV);
            dVolu = obj.mesh.computeDvolume(obj.quadrature);
            nGaus = obj.quadrature.ngaus;
            nElem = size(dVolu,2);

            nNodETs = size(dNTs,2);
            nDofETs = nNodETs*obj.test.ndimf;
            nNodETr = size(dNTr,2);
            nDofETr = nNodETr*obj.trial.ndimf;

            AxV = obj.A.evaluate(xV);
            K   = zeros(nDofETs,nDofETr,nElem);
            for iGaus = 1:nGaus
                for inode = 1:nNodETs
                    for jnode = 1:nNodETr
                        for iDimf = 1:obj.test.ndimf
                            idof = obj.test.ndimf*(inode-1)+iDimf;
                            jdof = obj.trial.ndimf*(jnode-1)+iDimf;
                            dV  = dVolu(iGaus,:)';
                            dNi = dNTs(:,inode,iGaus,:);
                            dNi = permute(dNi,[2 1 3 4]);
                            dNj = dNTr(:,jnode,iGaus,:);
                            Ag  = AxV(:,:,iGaus,:);
                            v   = pagemtimes(dNi,Ag);
                            v    = squeeze(pagemtimes(v,dNj));
                            K(idof, jdof, :)= squeeze(K(idof,jdof,:)) ...
                                + v(:).*dV;
                        end
                    end
                end
            end
            lhs = K;
        end
    end
end