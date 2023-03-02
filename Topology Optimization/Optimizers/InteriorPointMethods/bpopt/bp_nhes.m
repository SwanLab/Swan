classdef bp_nhes < handle
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
        xp
        sp
        dLBase
    end

    methods (Access = public)
        function obj = bp_nhes(cParams)
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
        end

        function computeLagrangianDerivativeBase(obj)
            dLB = bp_dL_dx(obj);
            dLB.compute();
            obj.dLBase = dLB.dL;
        end

        function computeLagrangianDerivative(obj)
            derivL = bp_dL_dx(obj);
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

function [h] = bp_nhes(bp,x,s,lam,bL,bU)

n = size(x,2);
m = size(bL,2);

% compute numerical 2nd deriv
dL_base = bp_dL_dx(bp,x,s,lam,bL,bU);
ep = 1e-5;
for i = 1:n,
    xp = x;
    xp(i) = x(i)+ep;
    dL1 = bp_dL_dx(bp,xp,s,lam,bL,bU);
    h(:,i) = [(dL1-dL_base)/ep]';
end
k = 0;
for i = 1:m,
    if (bU(i)>bL(i)),
       k = k + 1;
       sp = s;
       sp(i) = s(i)+ep;
       dL1 = bp_dL_dx(bp,x,sp,lam,bL,bU);
       h(:,n+k) = [(dL1-dL_base)/ep]';
    end
end
