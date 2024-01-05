classdef DesignVariableAcademic < handle
    
    properties (Access = public)
        value
        valueOld
    end
    
    methods (Access = public)
        function obj = DesignVariableAcademic(cParams)
            obj.init(cParams);
        end

        function update(obj,x)
            obj.value = x;
        end

        function updateOld(obj)
            obj.valueOld = obj.value;
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.value = cParams.x0;
        end
       
    end
    
end