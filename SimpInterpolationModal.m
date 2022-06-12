classdef SimpInterpolationModal < SimpInterpolation

    properties (Access = public)
        pFrac
    end

    methods (Access = public)

        function obj= SimpInterpolationModal(cParams)
            obj.init(cParams)
            obj.pExp = 3;
            obj.pFrac = 1/50;
            obj.computeSymbolicInterpolationFunctions();
        end

    end

    methods (Access = protected)

        function [k_low,k_high,dk_low,dk_high] = computeKappaSymbolicFunctionAndDerivative(obj)
            k0 = obj.matProp.kappa0;
            k1 = obj.matProp.kappa1;
            [k_low,k_high] = obj.interpolate(k0,k1);
            dk_low = diff(k_low);
            dk_high= diff(k_high);
        end

        function [mu_low,mu_high,dmu_low,dmu_high] = computeMuSymbolicFunctionAndDerivative(obj)
            mu0 = obj.matProp.mu0;
            mu1 = obj.matProp.mu1;
            [mu_low,mu_high] = obj.interpolate(mu0,mu1);
            dmu_low = diff(mu_low);
            dmu_high= diff(mu_high);
        end

        function [f_1,f_2] = interpolate(obj,f0,f1)
            p = obj.pExp;
            frac = obj.pFrac;
            [drho0,drho1] = obj.computeDensities();
            f_1 =  (drho1*frac)*f1; %(drho0*frac)*f0
            f_2 = (drho1^p)*f1; %  (drho0^p)*f0 + 
        end

        function [drho0,drho1] = computeDensities(obj)
            rho1 = obj.matProp.rho1;
            rho0 = obj.matProp.rho0;
            rho = sym('rho','real');
            inc  = rho1 - rho0;
            drho0 = (rho1 - rho)/inc;
            drho1 = (rho  - rho0)/inc;
        end

    end

end