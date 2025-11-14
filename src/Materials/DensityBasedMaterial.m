classdef DensityBasedMaterial < Material
    
    properties (Access = private)
       density 
       materialInterpolator
       dim
    end
    
    methods (Access = public)
        
        function obj = DensityBasedMaterial(cParams)
            obj.init(cParams)
            obj.density              = cParams.density;
            obj.materialInterpolator = cParams.materialInterpolator;
            obj.dim                  = cParams.dim;
        end
        
        function C = obtainTensor(obj)
            s.operation = @(xV) obj.evaluateNew(xV);
            s.mesh      = obj.mesh;
            s.ndimf     = repmat(obj.mesh.ndim,1,4);
            C = DomainFunction(s);
        end

        function dC = obtainTensorDerivative(obj)
            mI  = obj.materialInterpolator;
            rho = obj.density;
            [dmu,dkappa] = mI.computeConsitutiveTensorDerivative(rho);
            n1 = size(dmu,1);
            n2 = size(dmu,2);
            dC = cell(n1,n2);
            for i = 1:n1
                for j = 1:n2
                    s.operation = @(xV) obj.evaluateGradient(dmu{i,j},dkappa{i,j},xV);
                    s.mesh      = obj.mesh;
                    s.ndimf     = repmat(obj.mesh.ndim,1,4);
                    dC{i,j} = DomainFunction(s);
                end
            end
        end

        function setDesignVariable(obj,x)
            obj.density = x;
        end
        
    end

    methods (Access =protected)

        function C = evaluateNew(obj,xV)
            mI  = obj.materialInterpolator;
            rho = obj.density;
            [mu,kappa] = mI.computeConsitutiveTensor(rho);
            m = obj.createMaterial(mu,kappa);
            C = m.evaluate(xV);
        end

    end
    
    methods (Access = private)
        
        function m = createMaterial(obj,mu,kappa)
            s.type    = 'ISOTROPIC';
            s.ptype   = 'ELASTIC';
            s.ndim    = obj.computeNdim();
            s.shear   = mu;
            s.bulk    = kappa;
            s.mesh    = obj.mesh;
            m = Material.create(s);
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