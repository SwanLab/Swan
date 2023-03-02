classdef bp_verify_objgrad < handle
    properties (Access = public)
        okFirstDerivative
    end
    properties (Access = private)
        numericalDerivative
        exactDerivative
        bp
        x 
        s 
    end

    methods (Access = public)
        function obj = bp_verify_objgrad(cParams)
            obj.init(cParams);
        end

        function check(obj)
            obj.computeNumericalDerivative();
            obj.computeExactDerivative();
            obj.checkDeviation();
        end
    end
    methods (Access = private)
        function init(obj,cParams)
            obj.bp = cParams.bp;
            obj.x = cParams.x;
            obj.s = cParams.s;
        end

        function computeNumericalDerivative(obj)
            num = bp_nobjgrad(obj);
            num.compute();
            obj.numericalDerivative = num.gradient;
        end

        function computeExactDerivative(obj)
            exact = bp_objgrad(obj);
            exact.compute();
            obj.exactDerivative = exact.objGradient;
        end

        function checkDeviation(obj)
            deviation = max(max(abs(g_num-g_exact)));
            tol = 1e-5;
            if (deviation > tol)
                fprintf(1,'Jacobian error %9.2e exceeds tolerance %9.2e\n',dev_objgrad,tol) 
                obj.numericalDerivative
                obj.exactDerivative
                obj.okFirstDerivative = false;
            else
                obj.okFirstDerivative = true;
            end
        end
    end
end

function [ok_1st] = bp_verify_objgrad(bp,x,s)

% compute numerical 1st deriv
g_num = bp_nobjgrad(bp,x,s);

% compute exact 1st deriv
g_exact = bp_objgrad(bp,x,s);

% evaluate deviation
dev_objgrad = max(max(abs(g_num-g_exact)));

tol = 1e-5;
if (dev_objgrad>tol),
    fprintf(1,'Jacobian error %9.2e exceeds tolerance %9.2e\n',dev_objgrad,tol) 
    g_num
    g_exact
    ok_1st = false;
else
    ok_1st = true;
end
