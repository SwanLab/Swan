classdef NumericalHessianComputer < handle
    properties (Access = public)
        hessian
    end

    properties (Access = private)
        dL
        dL1
        bp 
        x 
        s 
        lambda 
        bU
        bL
        %xL
        %xU
        xp
        sp
        dLBase
    end

    methods (Access = public)
        function obj = NumericalHessianComputer(cParams)
            obj.init(cParams);
        end
        function compute(obj)
            obj.computeLagrangianDerivativeBase();
            obj.computeHessian();
        end
    end
    methods (Access = private)
        function init(obj,cParams)
            obj.bp = cParams.bp;
            obj.x = cParams.x;
            obj.s = cParams.s;
            obj.lambda = cParams.lam;
            obj.bL = cParams.bL;
            obj.bU = cParams.bU;
            %obj.xL = cParams.xL;
            %obj.xU = cParams.xU;
        end

        function computeLagrangianDerivativeBase(obj)
            u.x = obj.x;
            u.bp = obj.bp;
            u.s = obj.s;
            u.lam = obj.lambda;
            u.bL = obj.bL;
            u.bU = obj.bU;
            dLB = LagrangianDerivativeComputer(u);
            dLB.compute();
            obj.dLBase = dLB.dLagrangian;
        end

        function computeLagrangianDerivative(obj)
            u.x = obj.x;
            u.bp = obj.bp;
            u.s = obj.s;
            u.lam = obj.lambda;
            u.bL = obj.bL;
            u.bU = obj.bU;
            derivL = LagrangianDerivativeComputer(u);
            derivL.compute();
            obj.dL1 = derivL.dLagrangian;
        end

        function computeHessian(obj)
            n = size(obj.x,2);
            m = size(obj.bL,2);
            ep = 1e-5;
            for i = 1:n
                obj.xp = obj.x;
                obj.xp(i) = obj.x(i) + ep;
                obj.computeLagrangianDerivative();
                h(:,i) = [(obj.dL1 - obj.dLBase)/ep]';
            end
            k = 0;
            for i = 1:m
                if (obj.bU(i) > obj.bL(i))
                k = k + 1;
                obj.sp = obj.s;
                obj.sp(i) = obj.s(i) + ep;
                obj.computeLagrangianDerivative();
                h(:,n + k) = [(obj.dL1 - obj.dLBase)/ep]';
                end
            end
            obj.hessian = h;
        end
    end
end