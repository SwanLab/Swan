classdef ShFunc_Volume < Shape_Functional
    
    properties (Access = private)
        geometricVolume
        homogenizedVariablesComputer        
    end
    
    methods (Access = public)
        function obj = ShFunc_Volume(cParams)
            cParams.filterParams.quadratureOrder = 'CONSTANT';            
            obj.init(cParams);         
            obj.createHomogenizedVariablesComputer(cParams); 
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
            ngaus = obj.elemGradientSize.ngaus;
            gradient = repmat(gradient,1,ngaus);
            gradient = obj.filter.getP1fromP0(gradient);
            gradient = obj.Msmooth*gradient;            
            obj.gradient = gradient;            
        end
        
        function createHomogenizedVariablesComputer(obj,cParams)
            cP = cParams.materialInterpolationParams;
            cP.nelem = obj.elemGradientSize.nelem;
            cP.ngaus = obj.elemGradientSize.ngaus;
            h = HomogenizedVarComputer.create(cP);
            obj.homogenizedVariablesComputer = h;
        end                
        
        function updateHomogenizedMaterialProperties(obj,x)
            x0 = obj.filter.getP0fromP1(x);
            obj.homogenizedVariablesComputer.computeMatProp(x0);
        end        
        
    end
end

