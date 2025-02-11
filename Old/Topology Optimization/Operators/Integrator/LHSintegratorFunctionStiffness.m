classdef LHSintegratorFunctionStiffness < LHSintegrator

    properties (Access = private)     
        fun
    end

    methods (Access = public)

        function obj = LHSintegratorFunctionStiffness(cParams)
            obj@LHSintegrator(cParams)
            obj.fun = cParams.function;
        end

        function LHS = compute(obj)
            lhs = obj.computeElementalLHS();
            LHS = obj.assembleMatrix(lhs);
        end

    end

    methods (Access = protected)

        function lhs = computeElementalLHS(obj)
            xV = obj.quadrature.posgp;
            dNdxTs = obj.test.evaluateCartesianDerivatives(xV);
            dNdxTr = obj.trial.evaluateCartesianDerivatives(xV);
            dVolu = obj.mesh.computeDvolume(obj.quadrature);
            nGaus = obj.quadrature.ngaus;
            nElem = size(dVolu,2);

            nNodETs = size(dNdxTs,2);
            nDofETs = nNodETs*obj.trial.ndimf;
            nNodETr = size(dNdxTr,2);
            nDofETr = nNodETr*obj.trial.ndimf;
            lhs = zeros(nDofETs,nDofETr,nElem);

            fV  = obj.fun.evaluate(xV);
            for iGaus = 1:nGaus
                for inode= 1:nNodETs
                    for jnode= 1:nNodETr
                        for iunkn= 1:obj.test.ndimf
                            for idim = 1:obj.mesh.ndim
                                idof = obj.fun.ndimf*(inode-1)+iunkn;
                                jdof = obj.fun.ndimf*(jnode-1)+iunkn;
                                dvol = dVolu(iGaus,:);
                                dNi = squeeze(dNdxTs(idim,inode,iGaus,:));
                                dNj = squeeze(dNdxTr(idim,jnode,iGaus,:));
                                fVG = squeeze(fV(1,iGaus,:));
                                v = squeeze(fVG.*dNi.*dNj);
                                lhs(idof, jdof, :)= squeeze(lhs(idof,jdof,:)) ...
                                    + v(:).*dvol';
                            end
                        end
                    end
                end
            end
        end

    end


end