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
        lowerBounds
        upperBounds
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
            obj.lowerBounds    = cParams.lowerBounds;
            obj.upperBounds    = cParams.upperBounds;
        end

        function computeGradientReference(obj)
            lb          = obj.sMax;
            l           = obj.dualVariable.value;
            lZ          = obj.lowerBounds.Z;
            uZ          = obj.upperBounds.Z;
            nConstr     = obj.constraint.nSF;
            nnode       = obj.designVariable.mesh.nnodes;
            den         = nConstr+2*(nnode+obj.nSlack);
            obj.gradRef = max(lb,(sum(abs(l))+sum(abs(lZ))+sum(abs(uZ)))/den);
        end

        function computeFieldReference(obj)
            lb           = obj.sMax;
            lZ           = obj.lowerBounds.Z;
            uZ           = obj.upperBounds.Z;
            nnode        = obj.designVariable.mesh.nnodes;
            den          = 2*(nnode+obj.nSlack);
            obj.fieldRef = max(lb,(sum(abs(uZ))+sum(abs(lZ)))/den);
        end

        function computeErrorDueToGradients(obj)
            DJ            = obj.cost.gradient';
            Dg            = obj.constraint.gradient;
            l             = obj.dualVariable.value';
            lZ            = obj.lowerBounds.Z';
            uZ            = obj.upperBounds.Z';
            sD            = obj.gradRef;
            obj.errorGrad = max(abs(DJ + Dg*l - lZ + uZ))/sD;
        end

        function computeErrorDueToConstraint(obj)
            g               = obj.constraint.value;
            obj.errorConstr = max(abs(g));
        end

        function computeErrorDueToLowerBoundMargins(obj)
            x                 = obj.designVariable.value';
            lX                = obj.lowerBounds.X;
            s                 = obj.slack;
            lS                = obj.lowerBounds.S;
            lZ                = obj.lowerBounds.Z;
            nnode             = obj.designVariable.mesh.nnodes;
            e                 = ones(nnode+obj.nSlack,1);
            sC                = obj.fieldRef;
            obj.errorDesVarLB = max(abs(diag([x-lX s-lS])*diag(lZ)*e))/sC;
        end

        function computeErrorDueToUpperBoundMargins(obj)
            x                 = obj.designVariable.value';
            uX                = obj.upperBounds.X;
            s                 = obj.slack;
            uS                = obj.upperBounds.S;
            uZ                = obj.upperBounds.Z;
            nnode             = obj.designVariable.mesh.nnodes;
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