classdef bp_njac < handle
    properties (Access = public)
        jacobian
    end

    properties (Access = private)
        residualBase
        residual
        bp 
        x 
        s 
        bL 
        bU
        xp
        sp
    end

    methods (Access = public)
        function obj = bp_njac(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.computeResidualBase();
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

        function computeResidualBase(obj)
            rBase = bp_res(obj);
            rBase.compute();
            obj.residualBase = rBase.c;
        end

        function computeResidual(obj)
            r = bp_res(obj);
            r.compute();
            obj.residual = r.c;
        end

        function computeJacobian(obj)
            n = size(obj.x,2);
            m = size(obj.bL,2);
            ep = 1e-5;
            for i = 1:n
                obj.xp = obj.x;
                obj.xp(i) = obj.x(i) + ep;
                obj.computeResidual();
                J(:,i) = [(obj.residual - obj.residualBase)/ep]';
            end
            k = 0;
            for i = 1:m
                if (obj.bU(i) > obj.bL(i))
                k = k + 1;
                obj.sp = obj.s;
                obj.sp(i) = obj.s(i) + ep;
                obj.computeResidual();
                J(:,n + k) = [(obj.residual - obj.residualBase)/ep]';
                end
            end
            obj.jacobian = J;
        end
    end
end

function [J] = bp_njac(bp,x,s,bL,bU)

% compute numerical 1st deriv
n = size(x,2);
m = size(bL,2);

% compute numerical 1st deriv
r_base = bp_res(bp,x,s,bL,bU);
ep = 1e-5;
for i = 1:n,
    xp = x;
    xp(i) = x(i)+ep;
    r = bp_res(bp,xp,s,bL,bU);
    J(:,i) = [(r-r_base)/ep]';
end
k = 0;
for i = 1:m,
    if (bU(i)>bL(i)),
       k = k + 1;
       sp = s;
       sp(i) = s(i)+ep;
       r = bp_res(bp,x,sp,bL,bU);
       J(:,n+k) = [(r-r_base)/ep]';
    end
end