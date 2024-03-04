classdef DesignVariableAcademic < handle
    
    properties (Access = public)
        fun
        valueOld
    end
    
    methods (Access = public)
        function obj = DesignVariableAcademic(cParams)
            obj.init(cParams);
        end

        function update(obj,x)
            obj.fun.fValues = x;
        end

        function updateOld(obj)
            obj.valueOld = obj.fun.fValues;
        end

        function res = computeL2normIncrement(obj)
            incFun = obj.fun.fValues-obj.valueOld;
            nIncX  = norm(incFun);
            nX0    = norm(obj.valueOld);
            res    = nIncX/nX0;
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.fun.fValues = cParams.x0;
        end
       
    end
    
end