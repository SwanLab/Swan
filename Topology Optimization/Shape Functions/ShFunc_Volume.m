classdef ShFunc_Volume < ShapeFunctional
    
    properties (Access = public)
        geometricVolume
    end
    
    methods (Access = public)
        
        function obj = ShFunc_Volume(cParams)
            cParams.filterParams.quadratureOrder = 'LINEAR';
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
            densV = dens;
            volume = (sum(obj.dvolu,2)'*mean(densV,2));
            volume = volume/(obj.geometricVolume);
            obj.value = volume;
        end

        function v = getVariablesToPlot(obj)
            v{1} = obj.value*obj.value0;
        end

        function t = getTitlesToPlot(obj)
            t{1} = 'Volume non scaled';
        end
        
    end
    
    methods (Access = private)

        function computeGradient(obj)
            drho = obj.homogenizedVariablesComputer.drho;
            g     = drho;
            gf    = zeros(size(obj.Msmooth,1),obj.nVariables);
            q     = Quadrature.set(obj.designVariable.mesh.type);
            q.computeQuadrature('LINEAR');
            for ivar = 1:obj.nVariables
                gs           = g{ivar}/obj.geometricVolume;
                nelem        = size(gs,1);
                ngaus        = size(gs,2);
                s.fValues    = reshape(gs',[1,ngaus,nelem]);
                s.mesh       = obj.designVariable.mesh;
                s.quadrature = q;
                f            = FGaussDiscontinuousFunction(s);
                gradP1       = obj.gradientFilter.getP1Function(f,'LINEAR');
                gf(:,ivar)   = gradP1.fValues;
            end
            g = obj.Msmooth*gf;
            obj.gradient = g(:);
        end
        
        function updateHomogenizedMaterialProperties(obj)
            obj.designVariable.updateFunction();
            mesh      = obj.designVariable.mesh;
            f         = obj.designVariable.fun;
            fG        = obj.filter.getFGaussFunction(f,'QUADRATICMASS');
            xP0       = squeeze(fG.fValues);
            xf{1}     = reshape(xP0',[mesh.nelem,fG.quadrature.ngaus]);
            obj.homogenizedVariablesComputer.computeDensity(xf);
        end
        
    end
end

