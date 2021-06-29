classdef Volume_constraint < ShFunc_Volume
    
    properties (Access = private)
        VolumeTarget
        volum
    end
    
    methods (Access = public)        
        
        function  obj = Volume_constraint(cParams)
            obj@ShFunc_Volume(cParams);
        end
        
       function VolumeTarget = getVolumeTarget(obj)
            VolumeTarget = obj.target_parameters.Vfrac;
        end
        
       function computeFunctionAndGradient(obj)
           computeFunctionAndGradient@ShFunc_Volume(obj);
           obj.volum = obj.value;
           obj.value = obj.value/obj.getVolumeTarget() - 1;
           obj.gradient = obj.gradient/obj.getVolumeTarget();
       end
        
        function v = getVariablesToPlot(obj)
            v{1} = obj.volum;
        end
        
        function t = getTitlesToPlot(obj)
            t{1} = 'Volum';
        end       
        
        
    end
    
end

