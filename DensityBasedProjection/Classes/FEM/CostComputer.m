classdef CostComputer < Cost
    properties (Access = public)
        cost
    end
    methods (Access = public)
        function compute(obj)
            obj.cost  = abs(obj.force'*obj.displacement);
        end
    end
    methods (Access = private)
    end 
end