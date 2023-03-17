classdef JacobianVerifier < handle
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
        function obj = JacobianVerifier(cParams)
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
            u.bp = obj.bp;
            u.x = obj.x;
            u.s = obj.s;
            u.bL = obj.bL;
            u.bU = obj.bU;
            numerical = NumericalJacobianComputer(u);
            numerical.compute();
            obj.jacobianNum = numerical.jacobian;
        end

        function computeExactDerivative(obj)
            u.bp = obj.bp;
            u.x = obj.x;
            u.bL = obj.bL;
            u.bU = obj.bU;
            exact = JacobianComputer(u);
            exact.compute();
            obj.jacobianExact = exact.pd;
        end

        function checkDeviation(obj)
            deviation = max(max(abs(obj.jacobianNum - obj.jacobianExact)));
            tol = 1e-5;
            if (deviation > tol)
                fprintf(1,'Jacobian error %9.2e exceeds tolerance %9.2e\n',deviation,tol) 
                obj.jacobianNum
                obj.jacobianExact
                obj.okJacobian = false;
            else
                obj.okJacobian= true;
            end
        end
    end
end