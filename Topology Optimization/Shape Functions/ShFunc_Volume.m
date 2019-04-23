classdef ShFunc_Volume < Shape_Functional
    
    properties (Access = private)
       % Vfrac
        geometric_volume
    end
    
    methods (Access = public)
        function obj = ShFunc_Volume(settings)
            obj.init(settings);
            obj.geometric_volume = sum(obj.dvolu(:));
        end
        
        function computeCostAndGradient(obj,x)
            obj.computeCost(x)
            obj.computeGradient()
        end
        
        function computeCost(obj,x)
            xc(:,1) = x;
            density = obj.filter.getP0fromP1(xc);            
            obj.computeCostFromDensity(density);           
        end        
        
        function computeCostFromDensity(obj,dens)
            densV(:,1) = dens;
            volume = sum(sum(obj.dvolu,2)'*densV);
            volume = volume/(obj.geometric_volume);
            obj.value = volume;            
        end
        
    end
    
    methods (Access = private)

        function computeGradient(obj)
            gradient_volume = 1/(obj.geometric_volume);
            nnodes = size(obj.filter.diffReacProb.mesh.connec,1);
            gradient_volume = gradient_volume*ones(nnodes,size(obj.dvolu,2));
            gradient_volume = obj.filter.getP1fromP0(gradient_volume);
            gradient_volume = obj.Msmooth*gradient_volume;            
            obj.gradient = gradient_volume;            
        end
        
    end
end

