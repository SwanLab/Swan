classdef NumericalJacobianComputer < handle
    properties (Access = public)
        jacobian
    end

    properties (Access = private) 
        x 
        s 
        bL 
        bU
        bp
    end

    methods (Access = public)
        function obj = NumericalJacobianComputer(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.computeJacobian();
        end
    end
    methods (Access = private)
        function init(obj,cParams)
            obj.bp = cParams.bp;
            obj.x = cParams.x;
            obj.s = cParams.s;
            obj.bL = cParams.bL;
            obj.bU = cParams.bU;
        end

        function computeJacobian(obj)
            n = size(obj.x,2);
            m = size(obj.bL,2);
            ep = 1e-5;
            u.bp = obj.bp;
            u.x = obj.x;
            u.s = obj.s;
            u.bL = obj.bL;
            u.bU = obj.bU;
            residualBase = obj.computeResidual(u);
            for i = 1:n
                xp = u.x;
                xp(i) = u.x(i) + ep;
                u.x = xp;
                residual = obj.computeResidual(u);
                J(:,i) = [(residual - residualBase)/ep]';
            end
            k = 0;
            for i = 1:m
                if (obj.bU(i) > obj.bL(i))
                k = k + 1;
                sp = u.s;
                sp(i) = u.s(i) + ep;
                u.s = sp;
                residual = obj.computeResidual(u);
                J(:,n + k) = [(residual - residualBase)/ep]';
                end
            end
            obj.jacobian = J;
        end
    end
    methods (Static, Access = private)
        function residual = computeResidual(cParams)
            res = ResidualComputer(cParams);
            res.compute();
            residual = res.c;
        end
    end
end