classdef IPMDirectionComputer < handle

    properties (Access = public)
        gradients
        lowerSigma
        upperSigma
        H
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
        lowerBounds
        upperBounds
    end

    properties (Access = private)
        invLowerDesignVarMargin
        invUpperDesignVarMargin
        LHS
        RHS
        sol
    end

    methods (Access = public)
        function obj = IPMDirectionComputer(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.computeSigma();
            obj.updateHessian();
            obj.computeLHS();
            obj.computeRHS();
            obj.solveLinearSystem();
            obj.searchNewDirections();
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
            obj.nConstr        = cParams.nConstr;
            obj.nnode          = cParams.nX;
            obj.nSlack         = cParams.nSlack;
            obj.lowerBounds    = cParams.lowerBounds;
            obj.upperBounds    = cParams.upperBounds;
        end

        function computeSigma(obj)
            x   = obj.designVariable.value';
            s   = obj.slack;
            dLX = [x-obj.lowerBounds.X s-obj.lowerBounds.s];
            dUX = [obj.upperBounds.X-x obj.upperBounds.s-s];

            obj.invLowerDesignVarMargin = diag(1./dLX);
            obj.invUpperDesignVarMargin = diag(1./dUX);
            obj.lowerSigma = diag(obj.lowerBounds.Z./dLX);
            obj.upperSigma = diag(obj.upperBounds.Z./dUX);
        end

        function updateHessian(obj)
            obj.hessian(obj.nnode+1:obj.nnode+obj.nSlack) = 0;
            obj.H = obj.hessian + obj.lowerSigma + obj.upperSigma;
        end

        function computeLHS(obj)
            s.H          = obj.H;
            s.m          = obj.nConstr;
            s.constraint = obj.constraint;
            s.type       = 'IPMSymmetric';
            cLHS         = LHSintegrator.create(s);
            obj.LHS      = cLHS.compute();
        end

        function computeRHS(obj)
            s.baseVariables = obj.baseVariables;
            s.e             = ones(obj.nnode + obj.nSlack,1);
            s.nX            = obj.nnode;
            s.nSlack        = obj.nSlack;
            s.m             = obj.nConstr;
            s.cost          = obj.cost;
            s.constraint    = obj.constraint;
            s.lambda        = obj.dualVariable.value;
            s.lowerZ        = obj.lowerBounds.Z;
            s.upperZ        = obj.upperBounds.Z;
            s.invDiagdL     = obj.invLowerDesignVarMargin;
            s.invDiagdU     = obj.invUpperDesignVarMargin;
            s.type          = 'IPMSymmetric';
            cRHS            = RHSintegrator.create(s);
            obj.RHS         = cRHS.compute();
        end

        function solveLinearSystem(obj)
            s.type   = 'DIRECT';
            sLS      = Solver.create(s);
            obj.sol      = sLS.solve(-obj.LHS,obj.RHS);
        end

        function searchNewDirections(obj)
            mu                 = obj.baseVariables.mu;
            e                  = ones(obj.nnode + obj.nSlack,1);
            invDiagdLX         = obj.invLowerDesignVarMargin;
            invDiagdUX         = obj.invUpperDesignVarMargin;
            obj.gradients.dx   = obj.sol(1:obj.nnode,1);
            obj.gradients.ds   = obj.sol(obj.nnode+1:obj.nnode+obj.nSlack,1);
            obj.gradients.dlam = obj.sol(obj.nnode + obj.nSlack + 1:obj.nnode + obj.nSlack + obj.nConstr,1);
            obj.gradients.dzL  = mu*invDiagdLX*e-obj.lowerBounds.Z'-obj.lowerSigma*[obj.gradients.dx; obj.gradients.ds];
            obj.gradients.dzU  = mu*invDiagdUX*e-obj.upperBounds.Z'+obj.upperSigma*[obj.gradients.dx; obj.gradients.ds];
        end
    end
end