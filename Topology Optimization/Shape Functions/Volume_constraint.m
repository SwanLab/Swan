classdef Volume_constraint < ShFunc_Volume
    
    properties
        target
    end
    
    methods
        
        
        function  obj = Volume_constraint(settings)
            obj@ShFunc_Volume(settings);
        end
        
       function computef(obj,x)
           
           computef@ShFunc_Volume(obj,x);
           
           obj.value = obj.value/obj.Vfrac - 1;
           obj.gradient = obj.gradient/obj.Vfrac;
            
        end
        
        
    end
    
end

