classdef SimpInterpolationP3 < handle
    
   properties (Access = private)
        muFunc
        dmuFunc
        kappaFunc
        dkappaFunc
   end

   properties (Access = private)
        matA
        matB
        pExp
   end

    methods (Access = public)
        function obj = SimpInterpolationP3(cParams)
            obj.init(cParams)
            obj.computeSymbolicInterpolationFunctions();
        end

        function [mu,kappa] = computeConsitutiveTensor(obj,rho)
            mu      = CompositionFunction.create(obj.muFunc,rho);
            kappa   = CompositionFunction.create(obj.kappaFunc,rho);
        end

        function [dmu,dkappa] = computeConsitutiveTensorDerivative(obj,rho)
            dmu      = CompositionFunction.create(obj.dmuFunc,rho);
            dkappa   = CompositionFunction.create(obj.dkappaFunc,rho);
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.matA = cParams.matA;
            obj.matB = cParams.matB;
            obj.pExp = 3;
        end

        function computeSymbolicInterpolationFunctions(obj)
            [muS,dmuS,kS,dkS] = obj.computeSymbolicMuKappa();
            obj.muFunc        = matlabFunction(muS);
            obj.dmuFunc       = matlabFunction(dmuS);
            obj.kappaFunc     = matlabFunction(kS);
            obj.dkappaFunc    = matlabFunction(dkS);
        end

        function [muS,dmuS,kS,dkS] = computeSymbolicMuKappa(obj)
            [muS,dmuS] = obj.computeMuSymbolicFunctionAndDerivative();
            [kS,dkS]   = obj.computeKappaSymbolicFunctionAndDerivative();
        end

        function [mu,dmu] = computeMuSymbolicFunctionAndDerivative(obj)
            mu0 = obj.matA.shear;
            mu1 = obj.matB.shear;
            mu  = obj.interpolate(mu0,mu1);
            dmu = diff(mu);
        end

        function [k,dk] = computeKappaSymbolicFunctionAndDerivative(obj)
            k0 = obj.matA.bulk;
            k1 = obj.matB.bulk;
            k  = obj.interpolate(k0,k1);
            dk = diff(k);
        end

        function f = interpolate(obj,f0,f1)
            p = obj.pExp;
            [drho0,drho1] = obj.computeDensities();
            f = (drho0^p)*f0 + (drho1^p)*f1;
        end
    end

    methods (Static, Access = private)
        function [drho0,drho1] = computeDensities()
            rho   = sym('rho','real');
            drho0 = 1 - rho;
            drho1 = rho;
        end
    end
end