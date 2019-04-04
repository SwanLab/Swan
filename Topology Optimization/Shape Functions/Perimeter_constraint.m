classdef Perimeter_constraint < ShFunc_Perimeter
    
    properties
    end
    
    properties (Access = private)
        PerimeterTarget
    end
    
    methods
        
        function  obj = Perimeter_constraint(settings)
            obj@ShFunc_Perimeter(settings);
            obj.PerimeterTarget = settings.Perimeter_target; %should be target parameter?
        end
        
       function computeCostAndGradient(obj,x)
           computeCostAndGradient@ShFunc_Perimeter(obj,x);
           obj.value = obj.value/obj.PerimeterTarget - 1;
           obj.gradient = obj.gradient/obj.PerimeterTarget;
        end
        
        
    end
    
end

