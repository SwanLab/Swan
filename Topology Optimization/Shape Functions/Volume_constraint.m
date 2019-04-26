classdef Volume_constraint < ShFunc_Volume
    
    properties (Access = private)
        VolumeTarget
    end
    
    methods (Access = public)        
        
        function  obj = Volume_constraint(settings)
            obj@ShFunc_Volume(settings);
        end
        
       function VolumeTarget = getVolumeTarget(obj)
            VolumeTarget = obj.target_parameters.Vfrac;
        end
        
       function computeCostAndGradient(obj)
           computeCostAndGradient@ShFunc_Volume(obj);
           obj.value = obj.value/obj.getVolumeTarget() - 1;
           obj.gradient = obj.gradient/obj.getVolumeTarget();
        end
        
        
    end
    
end

