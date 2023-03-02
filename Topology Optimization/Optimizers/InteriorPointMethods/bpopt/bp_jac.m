% rearrange Jacobian
classdef bp_jac < handle
    properties (Access = public)
        pd
    end
    properties (Access = private)
        pdComputation
        x 
        bL 
        bU
    end

    methods (Access = public)
        function obj = bp_jac(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.computeBaseJacobian();
            obj.computeJacobian();
        end
    end
    methods (Access = private)
        function init(obj,cParams)
            obj.x = cParams.x;
            obj.bL = cParams.bL;
            obj.bU = cParams.bU;
        end
        function computeBaseJacobian(obj)
            baseJac = bp_jac_stub(obj);
            baseJac.create();
            obj.pdComputation = baseJac.pdComputation;
        end

        function computeJacobian(obj)
            m = size(obj.pdComputation,1);
            n = size(obj.pdComputation,2);
            k = 0;
            for i = 1:m
                if(obj.bU(i) > obj.bL(i))
                    k = k + 1;
                    obj.pdComputation(i,n+k) = -1;
                end
            end
            obj.pd = obj.pdComputation;
        end
    end
end
function [pd] = bp_jac(bp,x,bL,bU)
pd = bp_jac_stub(bp,x);
m = size(pd,1);
n = size(pd,2);

% add slack variables for inequality constraints
k = 0;
for i = 1:m,
   if(bU(i)>bL(i)),
       k = k + 1;
       pd(i,n+k) = -1;
   end
end