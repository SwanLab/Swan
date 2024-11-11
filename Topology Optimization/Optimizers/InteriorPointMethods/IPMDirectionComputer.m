classdef IPMDirectionComputer < handle

    properties (Access = public)
        gradients
        updatedHessian
    end

    properties (Access = private)
        cost
        constraint
        designVariable
        slack
        dualVariable
        baseVariables
        hessian
        nConstr
        nnode
        nSlack
        bounds
    end

    properties (Access = private)
        dLX
        dUX
        lowerSigma
        upperSigma
        LHS
        RHS
        sol
    end

    methods (Access = public)
        function obj = IPMDirectionComputer(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.computeBoundsMargins();
            obj.computeSigma();
            obj.updateHessian();
            obj.computeLHS();
            obj.computeRHS();
            obj.solveLinearSystem();
            obj.saveNewDirections();
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.cost           = cParams.cost;
            obj.constraint     = cParams.constraint;
            obj.designVariable = cParams.designVariable;
            obj.slack          = cParams.slack;
            obj.dualVariable   = cParams.dualVariable;
            obj.baseVariables  = cParams.baseVariables;
            obj.hessian        = cParams.hessian;
            obj.nConstr        = length(cParams.constraint.value);
            obj.nnode          = cParams.designVariable.fun.mesh.nnodes;
            obj.nSlack         = cParams.nSlack;
            obj.bounds         = cParams.bounds;
        end

        function computeBoundsMargins(obj)
            x       = obj.designVariable.fun.fValues';
            s       = obj.slack;
            obj.dLX = [x-obj.bounds.xLB s-obj.bounds.sLB];
            obj.dUX = [obj.bounds.xUB-x obj.bounds.sUB-s];
        end

        function computeSigma(obj)
            obj.lowerSigma = diag(obj.bounds.zLB./obj.dLX);
            obj.upperSigma = diag(obj.bounds.zUB./obj.dUX);
        end

        function updateHessian(obj)
            obj.hessian(obj.nnode+1:obj.nnode+obj.nSlack) = 0;
            obj.updatedHessian = obj.hessian + obj.lowerSigma + obj.upperSigma;
        end

        function computeLHS(obj)
            nSF     = obj.nConstr;
            H       = obj.updatedHessian;
            Dg      = obj.constraint.gradient;
            lhsX    = [H,Dg];
            lhsLam  = [Dg',zeros(nSF)];
            obj.LHS = [lhsX;lhsLam];
        end

        function computeRHS(obj)
            nSF     = obj.nConstr;
            l       = obj.dualVariable.fun.fValues;
            mu      = obj.baseVariables.mu;
            g       = obj.constraint.value;
            DJ      = obj.cost.gradient;
            Dg      = obj.constraint.gradient;
            rhsX    = (DJ + l*Dg'- mu./obj.dLX + mu./obj.dUX)';
            rhsLam  = g(1:nSF)';
            obj.RHS = [rhsX;rhsLam];
        end

        function solveLinearSystem(obj)
            s.type  = 'DIRECT';
            sLS     = Solver.create(s);
            obj.sol = sLS.solve(-obj.LHS,obj.RHS);
        end

        function saveNewDirections(obj)
            mu                 = obj.baseVariables.mu;
            lZ                 = obj.bounds.zLB';
            uZ                 = obj.bounds.zUB';
            lSig               = obj.lowerSigma;
            uSig               = obj.upperSigma;
            obj.gradients.dx   = obj.sol(1:obj.nnode,1);
            obj.gradients.ds   = obj.sol(obj.nnode+1:obj.nnode+obj.nSlack,1);
            obj.gradients.dlam = obj.sol(obj.nnode+obj.nSlack+1:end,1);
            desVarGrad         = [obj.gradients.dx; obj.gradients.ds];
            obj.gradients.dzL  = mu*obj.dLX'-lZ-lSig*desVarGrad;
            obj.gradients.dzU  = mu*obj.dUX'-uZ+uSig*desVarGrad;
        end
    end
end