classdef ShFunc_Volume < ShapeFunctional
    
    properties (Access = private)
        geometricVolume
    end
    
    methods (Access = public)
        function obj = ShFunc_Volume(cParams)
            cParams.filterParams.quadratureOrder = 'CONSTANT';            
            obj.init(cParams);         
            obj.geometricVolume = sum(obj.dvolu(:));
        end
        
        function computeCostAndGradient(obj,x)
            obj.updateHomogenizedMaterialProperties(x);            
            obj.computeCost()
            obj.computeGradient()
        end
        
        function computeCost(obj)          
            density = obj.homogenizedVariablesComputer.rho;
            obj.computeCostFromDensity(density);           
        end        
        
        function computeCostFromDensity(obj,dens)
            densV(:,1) = dens;
            volume = sum(sum(obj.dvolu,2)'*densV);
            volume = volume/(obj.geometricVolume);
            obj.value = volume;            
        end
        
    end
    
    methods (Access = private)

        function computeGradient(obj)
            drho = obj.homogenizedVariablesComputer.drho;
            gradient = drho/(obj.geometricVolume);
            gradient = obj.filter.getP1fromP0(gradient);
            gradient = obj.Msmooth*gradient;            
            obj.gradient = gradient;            
        end
        
        function updateHomogenizedMaterialProperties(obj,x)
            rho = obj.filter.getP0fromP1(x);
            obj.homogenizedVariablesComputer.computeDensity(rho);
        end        
        
    end
end

