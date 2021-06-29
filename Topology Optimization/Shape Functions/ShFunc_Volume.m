classdef ShFunc_Volume < ShapeFunctional
    
    properties (Access = public)
        geometricVolume
    end
    
    methods (Access = public)
        
        function obj = ShFunc_Volume(cParams)
            cParams.filterParams.quadratureOrder = 'CONSTANT';            
            obj.init(cParams);         
            obj.geometricVolume = sum(obj.dvolu(:));
        end
        
        function computeFunctionAndGradient(obj)
            obj.nVariables = obj.designVariable.nVariables;            
            obj.updateHomogenizedMaterialProperties();            
            obj.computeFunction();
            obj.computeGradient();
        end
        
        function computeFunction(obj)          
            density = obj.homogenizedVariablesComputer.rho;
            obj.computeFunctionFromDensity(density);  
        end        
        
        function computeFunctionFromDensity(obj,dens)
            densV(:,1) = dens;
            volume = sum(sum(obj.dvolu,2)'*densV);
            volume = volume/(obj.geometricVolume);
            obj.value = volume;            
        end
        
    end
    
    methods (Access = private)

        function computeGradient(obj)
            drho = obj.homogenizedVariablesComputer.drho;
            g = drho/(obj.geometricVolume);
            gf = zeros(size(obj.Msmooth,1),obj.nVariables);
            for ivar = 1:obj.nVariables
                gs = squeeze(g(:,ivar));
                gf(:,ivar) = obj.filter.getP1fromP0(gs);
            end
           % g = obj.Msmooth*gf;
            g = gf;
            obj.gradient = g(:);  
            
            obj.gradient = ones(size(g(:)))/obj.geometricVolume;
        end
        
        function updateHomogenizedMaterialProperties(obj)
            nx = length(obj.designVariable.value)/obj.designVariable.nVariables;
            x  = obj.designVariable.value;
            xf = cell(obj.designVariable.nVariables,1);            
            for ivar = 1:obj.nVariables
                i0 = nx*(ivar-1) + 1;
                iF = nx*ivar;
                xs = x(i0:iF);
                xf{ivar} = obj.filter.getP0fromP1(xs);
            end             
            obj.homogenizedVariablesComputer.computeDensity(xf);
        end        
        
    end
end

