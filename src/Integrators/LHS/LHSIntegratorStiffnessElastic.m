classdef LHSIntegratorStiffnessElastic < LHSIntegrator

    properties (Access = private)
        material
    end

    methods (Access = public)
        function obj = LHSIntegratorStiffnessElastic(cParams)
            obj@LHSIntegrator(cParams)
            obj.material = cParams.material;
        end

        function LHS = compute(obj)
            lhs = obj.computeElementalLHS();
            LHS = obj.assembleMatrix(lhs);
        end
    end

    methods (Access = protected)
        function lhs = computeElementalLHS(obj)
            xV    = obj.quadrature.posgp;
            SymN  = @(i) ShapeDerSym(obj.test,i);
            C     = obj.material;
            dE    = @(i,j) DDP(SymN(i),DDP(C,SymN(j)));
            dV     = obj.mesh.computeDvolume(obj.quadrature);
            nnodeE = obj.mesh.nnodeElem;
            ndim   = obj.mesh.ndim;
            ndofE  = nnodeE*ndim;
            nElem  = obj.mesh.nelem;
            lhs    = zeros(ndofE,ndofE,nElem);
            for i = 1:ndofE
                for j = 1:ndofE
                    dEval = dE(i,j).evaluate(xV);
                    dK    = dEval.*dV;
                    lhs(i,j,:) = squeeze(lhs(i,j,:))' + sum(dK,1);
                end
            end
        end
    end
end