classdef LHSIntegratorStiffnessModal < LHSIntegrator


    methods (Access = public)
        function obj = LHSIntegratorStiffnessModal(cParams)
            obj@LHSIntegrator(cParams)
        end

        function LHS = compute(obj,test,trial,mat)
            LHS = obj.computeElementalLHS(test,trial,mat);
%             LHS = obj.assembleMatrix(lhs);
        end
    end

    methods (Access = protected)
        %         function lhs = computeElementalLHS(obj,mat)
        %             xV     = obj.quadrature.posgp;
        %             C      = mat.evaluate(xV);
        %             nnodeE = obj.mesh.nnodeElem;
        %             ndim   = obj.mesh.ndim;
        %             ndofE  = nnodeE*ndim;
        %             lhs    = zeros(obj.trial.nbasis,obj.test.nbasis);
        %             dV     = obj.mesh.computeDvolume(obj.quadrature);
        %             for iBasis = 1: obj.trial.nbasis
        %                 uI     = obj.trial.basisFunctions{iBasis};
        %                 dNi     = ShapeDerTensor(uI);
        %                 dNi     = dNi.evaluate(xV);
        %                 for jBasis = 1: obj.test.nbasis
        %                     uJ     = obj.test.basisFunctions{jBasis};
        %                     dNj     = ShapeDerTensor(uJ);
        %                     dNj     = dNj.evaluate(xV);
        %                     for i = 1:ndofE
        %                         for j = 1:ndofE
        %                             dTrial     = squeezeParticular(dNi(:,:,i,:,:),3);
        %                             dTest      = squeezeParticular(dNj(:,:,j,:,:),3);
        %                             sigN       = pagetensorprod(C,dTrial,[3 4],[1 2],4,2);
        %                             de         = pagetensorprod(dTest,sigN,[1 2],[1 2],2,2);
        %                             dK         = de.*dV;
        %                             lhs(iBasis,jBasis) = lhs(iBasis,jBasis) + sum(dK,"all");
        %                         end
        %                     end
        %                 end
        %             end
        %         end
        function lhs = computeElementalLHS(obj,test,trial,mat)
            xV     = obj.quadrature.posgp;
            C      = mat.evaluate(xV);
            nnodeE = obj.mesh.nnodeElem;
            ndim   = obj.mesh.ndim;
            lhs    = zeros(obj.trial.nbasis,obj.test.nbasis);
            dV     = obj.mesh.computeDvolume(obj.quadrature);
            for iBasis = 1: obj.trial.nbasis
                uI      = trial.basisFunctions{iBasis};
                dNi     = Grad(uI);
                dNi     = dNi.evaluate(xV);
                for jBasis = 1: obj.test.nbasis
                    uJ      = test.basisFunctions{jBasis};
                    dNj     = Grad(uJ);
                    dNj     = dNj.evaluate(xV);

                    %                             dTrial     = squeezeParticular(dNi(:,:,i,:,:),3);
                    %                             dTest      = squeezeParticular(dNj(:,:,j,:,:),3);
                    sigN       = pagetensorprod(C,dNi,[3 4],[1 2],4,2);
                    de         = pagetensorprod(dNj,sigN,[1 2],[1 2],2,2);
                    dK         = de.*dV;
                    lhs(iBasis,jBasis) = lhs(iBasis,jBasis) + sum(dK,"all");
                end
            end
        end
    end
end