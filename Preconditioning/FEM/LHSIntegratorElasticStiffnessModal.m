classdef LHSIntegratorElasticStiffnessModal < LHSIntegrator


    methods (Access = public)
        function obj = LHSIntegratorElasticStiffnessModal(cParams)
            obj@LHSIntegrator(cParams)
        end

        function LHS = compute(obj,test,trial,mat)
            LHS = obj.computeElementalLHS(test,trial,mat);
%             LHS = obj.assembleMatrix(lhs);
        end
    end

    methods (Access = protected)
       
        function lhs = computeElementalLHS(obj,test,trial,mat)
            xV     = obj.quadrature.posgp;
            C      = mat.evaluate(xV);
            lhs    = zeros(obj.trial.nbasis,obj.test.nbasis);
            dV     = obj.mesh.computeDvolume(obj.quadrature);
            for iBasis = 1: obj.trial.nbasis
                uI      = trial.basisFunctions{iBasis};
                dNi     = SymGrad(uI);
                dNi     = dNi.evaluate(xV);
                for jBasis = 1: obj.test.nbasis
                    uJ      = test.basisFunctions{jBasis};
                    dNj     = SymGrad(uJ);
                    dNj     = dNj.evaluate(xV);
                    sigN       = pagetensorprod(C,dNi,[3 4],[1 2],4,2);
                    de         = pagetensorprod(dNj,sigN,[1 2],[1 2],2,2);
                    dK         = de.*dV;
                    lhs(iBasis,jBasis) = lhs(iBasis,jBasis) + sum(dK,"all");
                end
            end
        end

    end
end