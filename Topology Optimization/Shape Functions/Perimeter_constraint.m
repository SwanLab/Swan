classdef Perimeter_constraint < ShFunc_Perimeter
   
    properties (Access = private)
        perimeterTarget
    end
    
    methods
        
        function  obj = Perimeter_constraint(cParams)
            obj@ShFunc_Perimeter(cParams);
            obj.perimeterTarget = cParams.perimeterTarget; %should be target parameter?
        end
        
       function computeCostAndGradient(obj)
           computeCostAndGradient@ShFunc_Perimeter(obj);
           obj.value = obj.value/obj.perimeterTarget - 1;
           obj.gradient = obj.gradient/obj.perimeterTarget;
        end
        
        
    end
    
end

