<<<<<<<< HEAD:src/Integrators/LHS/LHSintegratorFunctionMass.m
classdef LHSintegratorFunctionMass < LHSIntegrator
========
classdef LHSIntegratorFunctionMass < LHSIntegrator
>>>>>>>> master:src/Integrators/LHS/LHSIntegratorFunctionMass.m

    properties (Access = private)
        fun
    end

    methods (Access = public)

<<<<<<<< HEAD:src/Integrators/LHS/LHSintegratorFunctionMass.m
        function obj = LHSintegratorFunctionMass(cParams)
            obj@LHSIntegrator(cParams)
            obj.fun = cParams.fun;
========
        function obj = LHSIntegratorFunctionMass(cParams)
            obj@LHSIntegrator(cParams);
            obj.fun = cParams.function;            
>>>>>>>> master:src/Integrators/LHS/LHSIntegratorFunctionMass.m
        end

        function LHS = compute(obj)
            lhs = obj.computeElementalLHS();
            LHS = obj.assembleMatrix(lhs);
        end

    end

    methods (Access = protected)

        function lhs = computeElementalLHS(obj)
            quad = obj.quadrature;
            xV   = quad.posgp;
            shapesTest  = obj.test.computeShapeFunctions(xV);
            shapesTrial = obj.trial.computeShapeFunctions(xV);
            dVolu  = obj.mesh.computeDvolume(quad);
            nGaus  = obj.quadrature.ngaus;
            nElem  = size(dVolu,2);

            nNodeTest  = size(shapesTest,1);
            nNodeTrial = size(shapesTrial,1);
            nDofTest   = nNodeTest*obj.test.ndimf;
            nDofTrial  = nNodeTrial*obj.trial.ndimf;

            fG = obj.fun.evaluate(quad.posgp);
            fG = squeezeParticular(fG,1);
            M = zeros(nDofTest, nDofTrial, nElem);
            % lhs = zeros(nDofTest/2, nDofTrial/2, nElem);
            % for igaus = 1:nGaus
            %     dv(1,1,:) = dVolu(igaus,:);
            %     Nv = shapesTest(:,igaus);
            %     Nu = shapesTrial(:,igaus);
            %     NvNu = Nv*Nu';
            %     Aij = bsxfun(@times,NvNu,dv);
            %     lhs = lhs + Aij;
            % end
            for igauss = 1 :nGaus
                for inode= 1:nNodeTest
                    for jnode= 1:nNodeTrial
                        for iunkn= 1:obj.test.ndimf
                       %     for junkn= 1:obj.trial.ndimf
                                fdv = fG(igauss,:).*dVolu(igauss,:);
                                idof = obj.test.ndimf*(inode-1)+iunkn;
                                jdof = obj.trial.ndimf*(jnode-1)+iunkn;
                                Ni = shapesTest(inode,igauss,:);
                                Nj = shapesTrial(jnode,igauss,:);
                                v = squeeze(Ni.*Nj.*fdv);
                                M(idof, jdof, :)= squeeze(M(idof,jdof,:)) ...
                                    + v(:);
                       %     end
                        end
                    end
                end
            end
            lhs = M;
        end

    end
<<<<<<<< HEAD:src/Integrators/LHS/LHSintegratorFunctionMass.m

end
========
 
end
>>>>>>>> master:src/Integrators/LHS/LHSIntegratorFunctionMass.m
