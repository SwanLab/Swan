classdef Perimeter_constraint < ShFunc_Perimeter
    
    properties
        target
    end
    
    methods
        
        
        function  obj = Perimeter_constraint(settings)
            obj@ShFunc_Perimeter(settings);
        end
        
       function computef(obj,x)
           
           computef@ShFunc_Perimeter(obj,x);
           
           obj.value = obj.value/obj.Perimeter_target - 1;
           obj.gradient = obj.gradient/obj.Perimeter_target;
            
        end
        
        
    end
    
end

