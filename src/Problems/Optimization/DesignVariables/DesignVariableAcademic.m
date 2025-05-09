classdef DesignVariableAcademic < handle & matlab.mixin.Copyable
    
    properties (Access = public)
        fun
        valuesOld
    end
    
    methods (Access = public)
        function obj = DesignVariableAcademic(cParams)
            obj.init(cParams);
        end

        function update(obj,x)
            obj.fun.fValues = x;
        end

        function updateOld(obj)
            obj.valuesOld = [obj.valuesOld,obj.fun.fValues];
        end

        function res = computeL2normIncrement(obj)
            incFun = obj.fun.fValues - obj.valuesOld(:,end);
            nIncX  = norm(incFun);
            nX0    = norm(obj.valuesOld(:,end));
            res    = nIncX/nX0;
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.fun.fValues = cParams.x0;
        end
       
    end
    
end