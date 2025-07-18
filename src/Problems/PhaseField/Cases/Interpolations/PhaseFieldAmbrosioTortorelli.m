classdef PhaseFieldAmbrosioTortorelli < handle
    
   properties (Access = private)
        shear
        bulk
        pExp
   end

    methods (Access = public)
        function obj = PhaseFieldAmbrosioTortorelli(cParams)
            obj.init(cParams)
        end

        function [mu,kappa] = computeConstitutiveTensorParams(obj,phi)
            mu    = obj.computeMuFunction(phi);
            kappa = obj.computeKappaFunction(phi);
        end

        function [dmu,dkappa] = computeConstitutiveTensorDerivativeParams(obj,phi)
            dmu    = obj.computeMuDerivative(phi);
            dkappa = obj.computeKappaDerivative(phi);
        end

        function [ddmu,ddkappa] = computeConstitutiveTensorSecondDerivativeParams(obj,phi)
            ddmu    = obj.computeMuSecondDerivative(phi);
            ddkappa = obj.computeKappaSecondDerivative(phi);
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            if isfield(cParams,'young') && isfield(cParams,'poisson')
                ndim = cParams.ndim;
                E    = cParams.young;
                nu   = cParams.poisson;
                obj.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E,nu,ndim);
                obj.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E,nu);
            elseif isfield(cParams,'bulk') && isfield(cParams,'shear')
                obj.shear = cParams.shear;
                obj.bulk  = cParams.bulk;
            end
            obj.pExp  = 2;
        end

        function mu = computeMuFunction(obj,phi)
            mu0 = obj.shear;
            mu  = obj.interpolate(phi,mu0);
        end

        function k = computeKappaFunction(obj,phi)
            k0 = obj.bulk;
            k  = obj.interpolate(phi,k0);
        end

        function dmu = computeMuDerivative(obj,phi)
            mu0 = obj.shear;
            dmu = obj.derive(phi,mu0);
        end

        function dk = computeKappaDerivative(obj,phi)
            k0 = obj.bulk;
            dk = obj.derive(phi,k0);
        end

        function ddmu = computeMuSecondDerivative(obj,phi)
            mu0  = obj.shear;
            ddmu = obj.derive2(phi,mu0);
        end

        function ddk = computeKappaSecondDerivative(obj,phi)
            k0  = obj.bulk;
            ddk = obj.derive2(phi,k0);
        end


        function f = interpolate(obj,phi,f0)
            p = obj.pExp;
            f = ((1-phi.fun).^p).*f0;

            d = 1;
            cw = 8/3; E = 210; Gc = 5e-3; l0 = 0.1;
            sigCrit = 1.5;            
            c = -2*(Gc/(cw*l0))*E*(1/sigCrit)^2;
            b = -10;
            a = -2*(c+1)-b/3;
            f = ((b/6).*phi.fun.^3 + (a/2).*phi.fun.^2 + c.*phi.fun + d).*f0;

            % fun = @(x) (b/6)*x^3 + (a/2).*x^2 + c.*x + d;
            % dfun = diff(fun)
        end
        
        function f = derive(obj,phi,f0)
            p = obj.pExp;
            f = -p*((1-phi.fun).^(p-1)).*f0;

            d = 1;
            cw = 8/3; E = 210; Gc = 5e-3; l0 = 0.1;
            sigCrit = 1.5;
            c = -2*(Gc/(cw*l0))*E*(1/sigCrit)^2;
            b = -10;
            a = -2*(c+1)-b/3;
            f = ((b/2).*phi.fun.^2 + (a).*phi.fun + c).*f0;
        end

        function f = derive2(obj,phi,f0)
            p = obj.pExp;
            f = p*(p-1)*((1-phi.fun).^(p-2)).*f0;

            d = 1;
            cw = 8/3; E = 210; Gc = 5e-3; l0 = 0.1;
            sigCrit = 1.5;
            c = -2*(Gc/(cw*l0))*E*(1/sigCrit)^2;
            b = -10;
            a = -2*(c+1)-b/3;
            f = (b.*phi.fun + a).*f0;
        end
    end

end