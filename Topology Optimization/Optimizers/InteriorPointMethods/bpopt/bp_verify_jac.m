classdef bp_verify_jac < handle
    properties (Access = public)
        okJacobian
    end
    properties (Access = private)
        jacobianNum
        jacobianExact
        x 
        s 
        bL 
        bU
        bp
    end

    methods (Access = public)
        function obj = bp_verify_jac(cParams)
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
            obj.bL = cParams.bL;
            obj.bU = cParams.bU;
            obj.bp = cParams.bp;
        end

        function computeNumericalDerivative(obj)
            numerical = bp_njac(obj);
            numerical.compute();
            obj.jacobianNum = numerical.jacobian;
        end

        function computeExactDerivative(obj)
            exact = bp_jac(obj);
            exact.compute();
            obj.jacobianExact = exact.pd;
        end

        function checkDeviation(obj)
            deviation = max(max(abs(obj.jacobianNum - obj.jacobianExact)));
            tol = 1e-5;
            if (deviation > tol)
                fprintf(1,'Jacobian error %9.2e exceeds tolerance %9.2e\n',dev_1st,tol) 
                obj.jacobianNum
                obj.jacobianExact
                obj.okJacobian = false;
            else
                obj.okJacobian= true;
            end
        end
    end
end