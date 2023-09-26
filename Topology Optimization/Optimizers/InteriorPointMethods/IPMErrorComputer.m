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
    end

    methods (Access = public)
        function obj = IPMErrorComputer(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.computeGradientReference();
            obj.computeFieldReference();
            obj.computeErrorDueToGradients();
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
    end
end