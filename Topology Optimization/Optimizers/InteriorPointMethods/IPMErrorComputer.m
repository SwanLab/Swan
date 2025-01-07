classdef IPMErrorComputer < handle

    properties (Access = public)
        error
    end

    properties (Access = private)
        sMax
        cost
        constraint
        designVariable
        dualVariable
        slack
        nSlack
        bounds
    end

    properties (Access = private)
        gradRef
        fieldRef
        errorGrad
        errorConstr
        errorDesVarLB
        errorDesVarUB
    end

    methods (Access = public)
        function obj = IPMErrorComputer(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.computeGradientReference();
            obj.computeFieldReference();
            obj.computeErrorDueToGradients();
            obj.computeErrorDueToConstraint();
            obj.computeErrorDueToLowerBoundMargins();
            obj.computeErrorDueToUpperBoundMargins();
            obj.computeLinfinityNorm();
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.sMax           = cParams.sMax;
            obj.cost           = cParams.cost;
            obj.constraint     = cParams.constraint;
            obj.designVariable = cParams.designVariable;
            obj.dualVariable   = cParams.dualVariable;
            obj.slack          = cParams.slack;
            obj.nSlack         = length(obj.slack);
            obj.bounds         = cParams.bounds;
        end

        function computeGradientReference(obj)
            lb          = obj.sMax;
            l           = obj.dualVariable.fun.fValues;
            lZ          = obj.bounds.zLB;
            uZ          = obj.bounds.zUB;
            nConstr     = length(obj.constraint.value);
            nnode       = obj.designVariable.fun.mesh.nnodes;
            den         = nConstr+2*(nnode+obj.nSlack);
            obj.gradRef = max(lb,(sum(abs(l))+sum(abs(lZ))+sum(abs(uZ)))/den);
        end

        function computeFieldReference(obj)
            lb           = obj.sMax;
            lZ           = obj.bounds.zLB;
            uZ           = obj.bounds.zUB;
            nnode        = obj.designVariable.fun.mesh.nnodes;
            den          = 2*(nnode+obj.nSlack);
            obj.fieldRef = max(lb,(sum(abs(uZ))+sum(abs(lZ)))/den);
        end

        function computeErrorDueToGradients(obj)
            DJ            = obj.cost.gradient';
            Dg            = obj.constraint.gradient;
            l             = obj.dualVariable.fun.fValues';
            lZ            = obj.bounds.zLB';
            uZ            = obj.bounds.zUB';
            sD            = obj.gradRef;
            obj.errorGrad = max(abs(DJ + Dg*l - lZ + uZ))/sD;
        end

        function computeErrorDueToConstraint(obj)
            g               = obj.constraint.value;
            obj.errorConstr = max(abs(g));
        end

        function computeErrorDueToLowerBoundMargins(obj)
            x                 = obj.designVariable.fun.fValues';
            lX                = obj.bounds.xLB;
            s                 = obj.slack;
            lS                = obj.bounds.sLB;
            lZ                = obj.bounds.zLB;
            nnode             = obj.designVariable.fun.mesh.nnodes;
            e                 = ones(nnode+obj.nSlack,1);
            sC                = obj.fieldRef;
            obj.errorDesVarLB = max(abs(diag([x-lX s-lS])*diag(lZ)*e))/sC;
        end

        function computeErrorDueToUpperBoundMargins(obj)
            x                 = obj.designVariable.fun.fValues';
            uX                = obj.bounds.xUB;
            s                 = obj.slack;
            uS                = obj.bounds.sUB;
            uZ                = obj.bounds.zUB;
            nnode             = obj.designVariable.fun.mesh.nnodes;
            e                 = ones(nnode+obj.nSlack,1);
            sC                = obj.fieldRef;
            obj.errorDesVarUB = max(abs(diag([uX-x uS-s])*diag(uZ)*e))/sC;
        end

        function computeLinfinityNorm(obj)
            e1        = obj.errorGrad;
            e2        = obj.errorConstr;
            e3        = obj.errorDesVarLB;
            e4        = obj.errorDesVarUB;
            e         = [e1,e2,e3,e4];
            obj.error = max(e);
        end
    end
end