classdef DensityBasedMaterialAnisotropic < handle
    
    properties (Access = private)
       density 
       materialInterpolator
       dim
       mesh
    end
    
    methods (Access = public)
        
        function obj = DensityBasedMaterialAnisotropic(cParams)
            obj.init(cParams)
        end
        
        function C = obtainTensor(obj)
            s.operation = @(xV) obj.evaluate(xV);
            s.mesh      = obj.mesh;
            C = DomainFunction(s);
        end

        function dC = obtainTensorDerivative(obj)
            s.operation = @(xV) obj.evaluateGradient(xV);
            s.mesh      = obj.mesh;
            dC{1} = DomainFunction(s);
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
        
        
        function C = evaluate(obj,xV)
            mI  = obj.materialInterpolator;
            rho = obj.density;
            C = mI.computeConsitutiveTensor(rho,xV);
        end
        
        function dC = evaluateGradient(obj,xV)
            mI  = obj.materialInterpolator;
            rho = obj.density;
            dC = mI.computeConsitutiveTensorDerivative(rho,xV);
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