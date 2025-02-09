classdef LHSintegrator_Mass_RT < LHSintegrator

    methods (Access = public)

        function obj = LHSintegrator_Mass_RT(cParams)
            obj@LHSintegrator(cParams)
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
            nDim   = obj.mesh.ndim;
            nNodeTest  = size(shapesTest,1);
            nNodeTrial = size(shapesTrial,1);
            nDofTest   = nNodeTest*obj.test.ndimf;
            nDofTrial  = nNodeTrial*obj.trial.ndimf;
            
            shapesTestMapped  = obj.test.mapFunction(shapesTest, xV);
            shapesTrialMapped = obj.trial.mapFunction(shapesTrial, xV);

            M = zeros(nDofTest, nDofTrial, nElem);
            for igauss = 1 :nGaus
                for inode= 1:nNodeTest
                    for jnode= 1:nNodeTrial
                        for iunkn= 1:obj.test.ndimf
                            idof = obj.test.ndimf*(inode-1)+iunkn;
                            jdof = obj.trial.ndimf*(jnode-1)+iunkn;
                            dvol = abs(dVolu(igauss,:))';
                            Ni = reshape(squeeze(shapesTestMapped(inode,igauss,:,:))',nDim,[])';
                            Nj = reshape(squeeze(shapesTrialMapped(jnode,igauss,:,:))',nDim,[])';
                            v  = sum(Ni.*Nj,2);
                            M(idof, jdof, :)= squeeze(M(idof,jdof,:)) + v(:).*dvol;
                        end
                    end
                end
            end
            lhs = M;
        end

    end

end