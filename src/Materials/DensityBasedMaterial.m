classdef DensityBasedMaterial < handle
    
    properties (Access = private)
       density 
       materialInterpolator
       dim
       mesh
       fibreOrientation
    end
    
    methods (Access = public)
        
        function obj = DensityBasedMaterial(cParams)
            obj.init(cParams)
        end
        
        function C = obtainTensor(obj)
            s.operation = @(xV) obj.evaluate(xV);
            s.mesh      = obj.mesh;
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
            obj.mesh                 = cParams.mesh;
            obj.fibreOrientation     = cParams.fibreOrientation
        end
        
        function m = createMaterial(obj,mu,kappa)
            s.type    = 'ANISOTROPIC';
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
            type = obj.fibreOrientation;
            C_voigt = Cvoigt.create(type);
            C = m.evaluate(xV,C_voigt);
            rhoEv = rho{1}.evaluate(xV);
            nGauss = size(rhoEv,2);
            nElem = size(rhoEv,3);
            rhoEv = reshape(rhoEv,[1 1 1 1 nGauss nElem]);
            C = (1-rhoEv.^3)*1e-3.*C + (rhoEv.^3).*C;
        end
        
        function dC = evaluateGradient(obj,dmu,dkappa,xV)
            mI  = obj.materialInterpolator;
            rho = obj.density;
            [mu,kappa] = mI.computeConsitutiveTensor(rho);
            m = obj.createMaterial(mu,kappa);
            type = obj.fibreOrientation;
            C_voigt = Cvoigt.create(type);
            C = m.evaluate(xV,C_voigt);
            rhoEv = rho{1}.evaluate(xV);
            nGauss = size(rhoEv,2);
            nElem = size(rhoEv,3);
            rhoEv = reshape(rhoEv,[1 1 1 1 nGauss nElem]);
            dC = (-3.*rhoEv.^2)*1e-3.*C + 3.*(rhoEv.^2).*C; 
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