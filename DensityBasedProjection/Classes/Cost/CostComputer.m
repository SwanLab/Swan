classdef CostComputer < handle
    properties (Access = public)
        cost
    end
    properties (Access = private)
        displacement 
        force
    end 
    methods (Access = public)
        function obj = CostComputer(cParams)
            obj.inputData(cParams)            
        end
        function compute(obj)
            obj.computeCost();
        end
    end
    methods (Access = private)
        function inputData(obj,cParams) 
            obj.displacement = cParams.displacement;
            obj.force = cParams.force;
        end      
        function computeCost(obj)
            obj.cost  = abs(obj.force'*obj.displacement);
        end 
    end 
end