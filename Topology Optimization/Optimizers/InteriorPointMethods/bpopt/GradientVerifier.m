classdef GradientVerifier < handle
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
        function obj = GradientVerifier(cParams)
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
            u.x = obj.x;
            u.s = obj.s;
            u.bp = obj.bp;
            num = NumericalGradientComputer(u);
            num.compute();
            obj.numericalDerivative = num.gradient;
        end

        function computeExactDerivative(obj)
            u.x = obj.x;
            u.bp = obj.bp;
            u.s = obj.s;
            exact = GradientComputer(u);
            exact.create();
            obj.exactDerivative = exact.objGradient;
        end

        function checkDeviation(obj)
            deviation = max(max(abs(obj.numericalDerivative - obj.exactDerivative)));
            tol = 1e-5;
            if (deviation > tol)
                fprintf(1,'Jacobian error %9.2e exceeds tolerance %9.2e\n',deviation,tol) 
                obj.numericalDerivative
                obj.exactDerivative
                obj.okFirstDerivative = false;
            else
                obj.okFirstDerivative = true;
            end
        end
    end
end