classdef IPMVariablesPusher < handle

    properties (Access = private)
        tau
        designVariable
        dualVariable
        slack
        bounds
        tol = 1e-2
    end

    methods (Access = public)
        function obj = IPMVariablesPusher(cParams)
            obj.init(cParams)
        end

        function updateBounds(obj,b)
            obj.bounds = b;
        end

        function x = replaceOutOfDesignVarBounds(obj)
            s.x   = obj.designVariable.fun.fValues;
            s.xUB = obj.bounds.xUB;
            s.xLB = obj.bounds.xLB;
            x     = obj.computeOutOfBounds(s);
        end

        function x = replaceOutOfSlackBounds(obj)
            s.x   = obj.slack;
            s.xUB = obj.bounds.sUB;
            s.xLB = obj.bounds.sLB;
            x     = obj.computeOutOfBounds(s);
        end

        function xNew = pushDesignVariable(obj,x)
            s.x   = x;
            s.xUB = obj.bounds.xUB;
            s.xLB = obj.bounds.xLB;
            xNew  = obj.pushVarsFunction(s);
        end

        function sNew = pushSlack(obj)
            s.x   = obj.slack;
            s.xUB = obj.bounds.sUB;
            s.xLB = obj.bounds.sLB;
            sNew  = obj.pushVarsFunction(s);
        end

        function bounds = pushDualVariableBounds(obj)
            bounds                 = obj.bounds;
            t                      = obj.tau;
            zLB                    = bounds.zNewLB;
            zUB                    = bounds.zNewUB;
            isZLNeg                = zLB<0;
            isZUNeg                = zUB<0;
            bounds.zNewLB(isZLNeg) = t*(zLB(isZLNeg));
            bounds.zNewUB(isZUNeg) = t*(zUB(isZUNeg));
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.tau            = cParams.tau;
            obj.designVariable = cParams.designVariable;
            obj.dualVariable   = cParams.dualVariable;
            obj.slack          = cParams.slack;
            obj.bounds         = cParams.bounds;
        end

        function x = computeOutOfBounds(obj,s)
            x          = s.x;
            xLB        = s.xLB;
            xUB        = s.xUB;
            isLower    = x<=xLB;
            isUpper    = x>=xUB;
            x(isLower) = min(xUB(isLower),xLB(isLower)+obj.tol);
            x(isUpper) = max(xLB(isUpper),xUB(isUpper)-obj.tol);
        end

        function x = pushVarsFunction(obj,s)
            t          = obj.tau;
            x          = s.x;
            xLB        = s.xLB;
            xUB        = s.xUB;
            isLower    = x<xLB;
            isUpper    = x>xUB;
            dxLB       = x(isLower) - xLB(isLower);
            dxUB       = xUB(isUpper) - x(isUpper);
            x(isLower) = xLB(isLower) + t*dxLB;
            x(isUpper) = xUB(isUpper) - t*dxUB;
        end

    end
end