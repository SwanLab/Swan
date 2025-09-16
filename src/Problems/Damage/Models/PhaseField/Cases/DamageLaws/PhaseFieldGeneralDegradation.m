classdef PhaseFieldGeneralDegradation < handle
    
   properties (Access = private)
        shear
        bulk
        pExp
   end

    methods (Access = public)
        function obj = PhaseFieldGeneralDegradation(cParams)
            obj.init(cParams)
        end

        function [mu,kappa] = computeConstitutiveTensor(obj,phi)
            mu    = obj.computeMuFunction(phi);
            kappa = obj.computeKappaFunction(phi);
        end

        function [dmu,dkappa] = computeConstitutiveTensorDerivative(obj,phi)
            dmu    = obj.computeMuDerivative(phi);
            dkappa = obj.computeKappaDerivative(phi);
        end

        function [ddmu,ddkappa] = computeConstitutiveTensorSecondDerivative(obj,phi)
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
            f = (1-phi.fun).*f0;
        end
        
        function f = derive(obj,phi,f0)
            f = 0.*phi.fun - f0;
        end

        function f = derive2(obj,phi,f0)
            f = 0.*phi.fun;
        end
    end

end