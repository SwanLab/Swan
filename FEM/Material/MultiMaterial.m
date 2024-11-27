classdef MultiMaterial < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
       density 
       materialInterpolator
       dim
    end
    
    methods (Access = public)
        
        function obj = MultiMaterial(cParams)
            obj.init(cParams)
        end
        
        function C = obtainTensor(obj)
            s.operation = @(xV) obj.evaluate(xV);
            C = DomainFunction(s);
        end

        function dC = obtainTensorDerivative(obj)
            mI  = obj.materialInterpolator;
            rho = obj.density;
            [dmu,dkappa] = mI.computeConsitutiveTensorDerivative(rho);
            for i = 1:size(dmu,1)
                for j = 1:size(dmu,1)
                    s.operation = @(xV) obj.evaluateGradient(dmu{i,j},dkappa{i,j},xV);
                    dC{i,j} = DomainFunction(s);
                end
            end
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
        end
        
        function m = createMaterial(obj,mu,kappa)
            s.type    = 'ISOTROPIC';
            s.ptype   = 'ELASTIC';
            s.ndim    = obj.computeNdim();
            s.shear   = mu;
            s.bulk    = kappa;
            m = Material.create(s);
        end
        
        function C = evaluate(obj,xV)
            mI  = obj.materialInterpolator;
            rho = obj.density;
            [mu,kappa] = mI.computeConsitutiveTensor(rho);
            m = obj.createMaterial(mu,kappa);
            C = m.evaluate(xV);
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