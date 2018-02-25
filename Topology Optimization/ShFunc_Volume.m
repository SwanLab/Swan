classdef ShFunc_Volume< Shape_Functional
    properties 
        Vfrac
    end
    methods 
        function obj = ShFunc_Volume(settings)
           obj@Shape_Functional(settings);
        end
        function Vfrac = get.Vfrac(obj)
            Vfrac = obj.target_parameters.Vfrac;
        end
        function computef(obj, x, ~, ~)
            mass = obj.filter.Msmooth;
            rho = obj.filter.getP0fromP1(x);     
            
            %compute volume
            geometric_volume = sum(mass(:));
           
            volume = sum(obj.filter.dvolu*rho);
            volume = volume/(geometric_volume*obj.Vfrac) - 1;
            
            %compute gradient
            gradient_volume = 1/(geometric_volume*obj.Vfrac);
            gradient_volume = gradient_volume*ones(size(obj.filter.connectivities,1),1);
            gradient_volume = obj.filter.getP1fromP0(gradient_volume);
            gradient_volume = mass*gradient_volume;
            
            obj.value = volume;
            obj.gradient = gradient_volume;
        end
    end
end
