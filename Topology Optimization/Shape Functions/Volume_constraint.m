classdef Volume_constraint < ShFunc_Volume
    
    properties (Access = private)
        VolumeTarget
    end
    
    methods
        
        
        function  obj = Volume_constraint(settings)
            obj@ShFunc_Volume(settings);
        end
        
       function VolumeTarget = get.VolumeTarget(obj)
            VolumeTarget = obj.target_parameters.Vfrac;
        end
        
       function computeCostAndGradient(obj,x)
           computeCostAndGradient@ShFunc_Volume(obj,x);
           obj.value = obj.value/obj.VolumeTarget - 1;
           obj.gradient = obj.gradient/obj.VolumeTarget;
        end
        
        
    end
    
end

