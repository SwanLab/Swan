classdef bp_lagrangian < handle

    properties (Access = public)
        lagrangian
        object
        residual
    end
    properties (Access = private)
        bp
        x
        s
        lam 
        bL 
        bU
    end

    methods (Access = public)
        function obj = bp_lagrangian(cParams)
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
            objective = bp_obj(u);
            objective.compute();
            obj.object = objective.Object;
        end

        function computeResidual(obj)
            u.bp = obj.bp;
            u.x = obj.x;
            u.s = obj.s;
            u.bL = obj.bL;
            u.bU = obj.bU;
            res = bp_res(u);
            res.compute();
            obj.residual = res.Residual;
        end

        function computeLagrangian(obj)
            obj.lagrangian = obj.object + obj.residual*obj.lam';
        end
    end
end