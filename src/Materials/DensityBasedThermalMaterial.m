classdef DensityBasedThermalMaterial < handle
    
    properties (Access = private)
       density 
       thermalmaterialInterpolator
       dim
       mesh
    end
    
    methods (Access = public)
        
        function obj = DensityBasedThermalMaterial(cParams)
            obj.init(cParams)
        end
        
        function kappa = obtainConductivity(obj)
            s.operation = @(xV) obj.evaluate(xV);
            s.mesh      = obj.mesh;
            kappa = DomainFunction(s);
        end

        function dkappa = obtainConductivityDerivative(obj)
            mI  = obj.thermalmaterialInterpolator;
            rho = obj.density;
            dkappa = mI.computeConductivityDerivative(rho);
            s.operation = @(xV) obj.evaluateGradient(xV);
            s.mesh      = obj.mesh;
            dkappa = DomainFunction(s);
        end

        function setDesignVariable(obj,x)
            obj.density = x;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.density              = cParams.density;
            obj.materialInterpolator = cParams.materialInterpolator;
            obj.dim                  = cParams.dim;
            obj.mesh                 = cParams.mesh;
        end
        
        function kappa = evaluate(obj,xV)
            mI  = obj.materialInterpolator;
            rho = obj.density;
            [mu,kappa] = mI.computeConsitutiveTensor(rho);
            m = obj.createMaterial(mu,kappa);
            kappa = m.evaluate(xV);
        end
        
        function dC = evaluateGradient(obj,dmu,dkappa,xV)
            m  = obj.createMaterial(dmu,dkappa);
            dC = m.evaluate(xV);
        end

        function ndim = computeNdim(obj)
            switch obj.dim
                case '2D'
                  ndim = 2;
                case '3D'
                  ndim = 3;
            end
        end
        
    end
    
end