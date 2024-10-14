classdef LHSintegratorAnisotropicMass < LHSintegrator

    properties (Access = private)
        A % change name to not get confused
    end

    methods (Access = public)

        function obj = LHSintegratorAnisotropicMass(cParams)
            obj@LHSintegrator(cParams);
            obj.A = cParams.A;
        end

        function LHS = compute(obj)
            lhs = obj.computeElementalLHS();
            LHS = obj.assembleMatrix(lhs);
        end

    end

    methods (Access = private)
        
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

            M = zeros(nDofTest, nDofTrial, nElem);
            Aij  = obj.A;

            for igauss = 1 :nGaus
                for inode= 1:nNodeTest
                    for jnode= 1:nNodeTrial
                        Ni = shapesTest(inode,igauss,:);
                        Nj = shapesTrial(jnode,igauss,:);
                        Nveci = [Ni;Ni];
                        Nvecj = [Nj;Nj];
                        for iunkn= 1:obj.test.ndimf
                            for junkn= 1:obj.trial.ndimf
                                idof = obj.test.ndimf*(inode-1)+iunkn;
                                jdof = obj.trial.ndimf*(jnode-1)+junkn;
                                dvol = dVolu(igauss,:);
                                v = Nveci(iunkn)*Aij(iunkn,junkn)*Nvecj(junkn);
                                %v   = Nveci'*obj.A*Nvecj;
                                M(idof, jdof, :)= squeeze(M(idof,jdof,:)) ...
                                    + v(:).*dvol';
                            end
                        end
                    end
                end
            end
            lhs = M;

        end

    end

end