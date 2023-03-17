%% BPOPT Solver: Obtain Jacobian (1st derivatives) of Equations
classdef JacobianObtainer < handle
    properties (Access = public)
        pdComputation
    end
    properties (Access = private)
        xC
        bp
    end

    methods (Access = public)
        function obj = JacobianObtainer(cParams)
            obj.init(cParams);
        end

        function create(obj)
            obj.createCase();
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.xC = cParams.x;
            obj.bp = cParams.bp;
        end
        function createCase(obj)
            x = obj.xC;
            switch obj.bp.prob
            case(0)
                y = apm_details(obj.bp.server,obj.bp.app,x);
                pd = y.jac;
                obj.pdComputation = pd;
            case(1)
                pd(1,1) = x(2)*x(3)*x(4);
                pd(1,2) = x(1)*x(3)*x(4);
                pd(1,3) = x(1)*x(2)*x(4);
                pd(1,4) = x(1)*x(2)*x(3);
                pd(2,1) = 2*x(1);
                pd(2,2) = 2*x(2);
                pd(2,3) = 2*x(3);
                pd(2,4) = 2*x(4);
                obj.pdComputation = pd;
            case(2)
                n = size(x,2);
                pd = ones(1,n);
                obj.pdComputation = pd;
            case(3)
                %n = size(x,2);
                pd(1,1) = 0.1;
                pd(1,2) = -1;
                obj.pdComputation = pd;
            case(4)
                pd(1,1) = x(2);
                pd(1,2) = x(1);
                pd(2,1) = 2*x(1);
                pd(2,2) = 2*x(2);
                obj.pdComputation = pd;
            case(5)
                pd(1,1) = 0.1;
                pd(1,2) = -1.0;
                pd(2,1) = -10.0;
                pd(2,2) = 1.0;  
                obj.pdComputation = pd;      
            case(6)
                pd(1,1) = 2;
                pd(1,2) = 1;
                pd(2,1) = 1;
                pd(2,2) = 2;
                obj.pdComputation = pd;   
            case(7)
                pd(1,1) = 1;
                pd(1,2) = 1;
                obj.pdComputation = pd;
            case(8)
                pd(1,1) = -2 * x(1);
                pd(1,2) = -1;
                obj.pdComputation = pd;
            case(9)
                pd(1,1) = 2 * x(1);
                obj.pdComputation = pd;
            case(10)
                pd(1,1) = 0.1;
                pd(1,2) = -1.0;
                obj.pdComputation = pd;
            end
        end
    end
end