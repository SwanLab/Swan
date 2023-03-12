classdef bp_theta < handle
    properties (Access = public)
        theta
    end
    properties (Access = private)
        bp
        x 
        s 
        bL 
        bU 
        residual
    end

    methods (Access = public)
        function obj = bp_theta(cParams)
            obj.init(cParams);
        end
        
        function compute(obj)
            obj.computeResidual();
            obj.computeTheta();
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

        function computeResidual(obj)
            res = bp_res(obj);
            res.compute();
            obj.residual = res.c;
        end

        function computeTheta(obj)
            obj.theta = sum(abs(obj.residual));
        end
    end
end