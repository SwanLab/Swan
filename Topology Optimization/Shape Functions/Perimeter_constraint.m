classdef Perimeter_constraint < ShFunc_Perimeter
    
    properties
    end
    
    properties (Access = private)
        PerimeterTarget
    end
    
    methods
        
        function  obj = Perimeter_constraint(cParams)
            obj@ShFunc_Perimeter(cParams);
            obj.PerimeterTarget = cParams.Perimeter_target; %should be target parameter?
        end
        
       function computeCostAndGradient(obj)
           computeCostAndGradient@ShFunc_Perimeter(obj);
           obj.value = obj.value/obj.PerimeterTarget - 1;
           obj.gradient = obj.gradient/obj.PerimeterTarget;
        end
        
        
    end
    
end

