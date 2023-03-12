classdef bp_verify_hes < handle
    properties (Access = public)
        okHessian
    end
    properties (Access = private)
        numericalHessian
        exactHessian
        x 
        s 
        lambda 
        bL
        bU
    end

    methods (Access = public)
        function obj = bp_verify_hes(cParams)
            obj.init(cParams)
        end

        function compute(obj)
            obj.computeNumericalDerivative();
            obj.computeExactDerivative();
            obj.checkDeviation();
        end
    end
    methods (Access = private)
        function init(obj,cParams)
            obj.x = cParams.x;
            obj.s = cParams.s;
            obj.lambda = cParams.lam;
            obj.bL = cParams.bL;
            obj.bU = cParams.bU;
        end

        function computeNumericalDerivative(obj)
            numerical = bp_nhes(obj);
            numerical.compute();
            obj.numericalHessian = numerical.hessian;
        end

        function computeExactDerivative(obj)
            exact = bp_hes(obj);
            exact.compute();
            obj.exactHessian = exact.hess;
        end

        function checkDeviation(obj)
            deviation = max(max(abs(obj.numericalHessian - obj.exactHessian)));
            tol = 1e-3;
            if (deviation > tol)
                fprintf(1,'Hessian error %9.2e exceeds tolerance %9.2e\n',dev_2nd,tol) 
                obj.numericalHessian
                obj.exactHessian
                obj.okHessian = false;
            else
                obj.okHessian = true;
            end
        end
    end
end