% rearrange Jacobian
classdef JacobianComputer < handle
    properties (Access = public)
        pd
    end
    properties (Access = private)
        pdComputation
        x 
        bL 
        bU
        bp
    end

    methods (Access = public)
        function obj = JacobianComputer(cParams)
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
            obj.bp = cParams.bp;
            obj.bL = cParams.bL;
            obj.bU = cParams.bU;
        end
        
        function computeBaseJacobian(obj)
            u.bp = obj.bp;
            u.x = obj.x;
            baseJac = JacobianObtainer(u);
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