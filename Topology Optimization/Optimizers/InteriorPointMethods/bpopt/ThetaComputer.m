classdef ThetaComputer < handle
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
        function obj = ThetaComputer(cParams)
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
            u.bp = obj.bp;
            u.x = obj.x;
            u.s = obj.s;
            u.bL = obj.bL;
            u.bU = obj.bU;
            res = ResidualComputer(u);
            res.compute();
            obj.residual = res.c;
        end

        function computeTheta(obj)
            obj.theta = sum(abs(obj.residual));
        end
    end
end