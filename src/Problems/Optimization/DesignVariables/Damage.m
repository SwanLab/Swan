classdef Damage < DesignVariable

    methods (Access = public)

        function obj = Damage(cParams)
            obj.nVariables = 1;
            obj.init(cParams);
        end

        function update(obj,value)
            obj.fun.setFValues(value);
        end

        function plot(obj)
            Plot(obj.fun);
        end
    
    end
    
end

