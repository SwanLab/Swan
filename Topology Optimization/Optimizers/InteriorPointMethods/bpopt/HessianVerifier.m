classdef HessianVerifier < handle
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
        bp
    end

    methods (Access = public)
        function obj = HessianVerifier(cParams)
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
            obj.bp = cParams.bp;
        end

        function computeNumericalDerivative(obj)
            u.bp = obj.bp;
            u.x = obj.x;
            u.s = obj.s;
            u.lam = obj.lambda;
            u.bL = obj.bL;
            u.bU = obj.bU;
            numerical = NumericalHessianComputer(u);
            numerical.compute();
            obj.numericalHessian = numerical.hessian;
        end

        function computeExactDerivative(obj)
            u.x = obj.x;
            u.s = obj.s;
            u.lam = obj.lambda;
            u.bp = obj.bp;
            exact = HessianComputer(u);
            exact.create();
            obj.exactHessian = exact.hess;
        end

        function checkDeviation(obj)
            deviation = max(max(abs(obj.numericalHessian - obj.exactHessian)));
            tol = 1e-3;
            if (deviation > tol)
                fprintf(1,'Hessian error %9.2e exceeds tolerance %9.2e\n',deviation,tol) 
                obj.numericalHessian
                obj.exactHessian
                obj.okHessian = false;
            else
                obj.okHessian = true;
            end
        end
    end
end