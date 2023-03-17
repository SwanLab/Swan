classdef LagrangianComputer < handle

    properties (Access = public)
        lagrangian
    end
    properties (Access = private)
        bp
        x
        s
        lam 
        bL 
        bU
        object
        residual
    end

    methods (Access = public)
        function obj = LagrangianComputer(cParams)
            obj.init(cParams)
        end
        function compute(obj)
            obj.computeObject();
            obj.computeResidual();
            obj.computeLagrangian();
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.bp = cParams.bp;
            obj.x = cParams.x;
            obj.s = cParams.s;
            obj.lam = cParams.lam;
            obj.bL = cParams.bL;
            obj.bU = cParams.bU;
        end

        function computeObject(obj)
            u.bp = obj.bp;
            u.x = obj.x;
            objective = ObjectiveFunctionComputer(u);
            objective.compute();
            obj.object = objective.objectiveFunc;
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

        function computeLagrangian(obj)
            obj.lagrangian = obj.object + obj.residual*obj.lam';
        end
    end
end