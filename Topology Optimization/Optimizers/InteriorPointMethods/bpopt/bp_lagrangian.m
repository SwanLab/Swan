classdef bp_lagrangian < handle

    properties (Access = public)
        Lagrangian
        Object
        Residual
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
        
        end

        function computeObject(obj)
            Objective = bp_obj(obj);
            Objective.compute();
            obj.Object = Objective.Object;
        end

        function computeResidual(obj)
            Res = bp_res(obj);
            Res.compute();
            obj.Residual = Res.Residual;
        end

        function computeLagrangian(obj)
            obj.Lagrangian = obj.Object + obj.Residual*lam';
        end
    end
end
function [L] = bp_lagrangian(bp,x,s,lam,bL,bU)

    L = bp_obj(bp,x) + bp_res(bp,x,s,bL,bU)*lam'; 
end
