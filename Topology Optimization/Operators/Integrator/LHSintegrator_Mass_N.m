classdef LHSintegrator_Mass_N < LHSintegrator

    methods (Access = public)

        function obj = LHSintegrator_Mass_N(cParams)
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
            nNodeTest  = size(shapesTest,1);
            nNodeTrial = size(shapesTrial,1);
            nDofTest   = nNodeTest*obj.test.ndimf;
            nDofTrial  = nNodeTrial*obj.trial.ndimf;
            
            sides = obj.test.computeSidesOrientation();
            J = obj.mesh.computeJacobian(0);
            R = [0 -1; 1 0];
            JGlob = pagemtimes(pagemtimes(R,J),R');
            Jdet = obj.mesh.computeJacobianDeterminant(xV);

            M = zeros(nDofTest, nDofTrial, nElem);
            for igauss = 1 :nGaus
                for inode= 1:nNodeTest
                    for jnode= 1:nNodeTrial
                        for iunkn= 1:obj.test.ndimf
                            idof = obj.test.ndimf*(inode-1)+iunkn;
                            jdof = obj.trial.ndimf*(jnode-1)+iunkn;
                            dvol = abs(dVolu(igauss,:))';
                            Jd(1,1,1,:) = 1./Jdet(igauss,:);
                            Ni = pagemtimes(squeeze(shapesTest(inode,igauss,:,:))',JGlob).*Jd;
                            Nj = pagemtimes(squeeze(shapesTrial(jnode,igauss,:,:))',JGlob).*Jd;
                            v = squeeze(pagemtimes(Ni,pagetranspose(Nj)));
                            M(idof, jdof, :)= squeeze(M(idof,jdof,:)) ...
                                + v(:).*dvol.*sides(:,inode).*sides(:,jnode);
                        end
                    end
                end
            end
            lhs = M;
        end

    end

end