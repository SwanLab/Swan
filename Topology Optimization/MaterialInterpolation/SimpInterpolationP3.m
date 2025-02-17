classdef SimpInterpolationP3 < handle
    
   properties (Access = private)
        matA
        matB
        pExp
   end

    methods (Access = public)
        function obj = SimpInterpolationP3(cParams)
            obj.init(cParams)
        end

        function [mu,kappa] = computeConsitutiveTensor(obj,rho)
            mu    = obj.computeMuFunction(rho{1});
            kappa = obj.computeKappaFunction(rho{1});
        end

        function [dmu,dkappa] = computeConsitutiveTensorDerivative(obj,rho)
            dmu{1}    = obj.computeMuDerivative(rho{1});
            dkappa{1} = obj.computeKappaDerivative(rho{1});
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.matA = cParams.matA;
            obj.matB = cParams.matB;
            obj.pExp = 8;
        end

        function mu = computeMuFunction(obj,rho)
            mu0 = obj.matA.shear;
            mu1 = obj.matB.shear;
            mu  = obj.interpolate(rho,mu0,mu1);
        end

        function k = computeKappaFunction(obj,rho)
            k0 = obj.matA.bulk;
            k1 = obj.matB.bulk;
            k  = obj.interpolate(rho,k0,k1);
        end

        function dmu = computeMuDerivative(obj,rho)
            mu0  = obj.matA.shear;
            mu1  = obj.matB.shear;
            dmu  = obj.derive(rho,mu0,mu1);
        end

        function dk = computeKappaDerivative(obj,rho)
            k0  = obj.matA.bulk;
            k1  = obj.matB.bulk;
            dk  = obj.derive(rho,k0,k1);
        end

        function f = interpolate(obj,rho,f0,f1)
            p = obj.pExp;
            [drho0,drho1] = obj.computeDensities(rho);
            f = (drho0^p)*f0 + (drho1^p)*f1;
        end
        
        function f = derive(obj,rho,f0,f1)
            p = obj.pExp;
            [drho0,drho1] = obj.computeDensities(rho);
            f = -p*(drho0^(p-1))*f0 + p*(drho1^(p-1))*f1;
        end
    end

    methods (Static, Access = private)
        function [drho0,drho1] = computeDensities(rho)
            drho0 = 1 - rho;
            drho1 = rho;
        end
    end
end