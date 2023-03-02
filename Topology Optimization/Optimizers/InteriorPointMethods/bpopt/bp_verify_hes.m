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


function [ok_2nd] = bp_verify_hes(bp,x,s,lam,bL,bU)

% compute numerical 2nd deriv
h_num = bp_nhes(bp,x,s,lam,bL,bU);

% compute exact 2nd deriv
h_exact = bp_hes(bp,x,s,lam);

% evaluate deviation
dev_2nd = max(max(abs(h_num-h_exact)));

tol = 1e-3;
if (dev_2nd>tol),
    fprintf(1,'Hessian error %9.2e exceeds tolerance %9.2e\n',dev_2nd,tol) 
    h_num
    h_exact
    ok_2nd = false;
else
    ok_2nd = true;
end
