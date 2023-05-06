classdef DualUpdater_IPM < handle

    properties (Access = private)
        dualVariable
        constraint
        nConstr
        constraintCase
        cost
        lam
        zLa
        zUa
        alphaDualMax
    end

    methods (Access = public)
        function obj = DualUpdater_IPM(cParams)
            obj.init(cParams);
        end

        function compute(obj,lz,uz)
            c   = obj.constraint.gradient';
            g = obj.cost.gradient;
            l = pinv(full(c*c'))*c*(lz'- uz'- g');
            obj.dualVariable.value = l';
        end

        function obj = update(obj,tau,g)
            obj.dualVariable.value = obj.dualVariable.value + tau * g';
        end
    end


    methods (Access = private)
        function init(obj,cParams)
            obj.dualVariable = cParams.dualVariable;
            obj.constraint     = cParams.constraint;
            obj.constraintCase = cParams.constraintCase;
            obj.nConstr        = cParams.constraint.nSF;
            obj.cost           = cParams.cost;
        end
    end
end