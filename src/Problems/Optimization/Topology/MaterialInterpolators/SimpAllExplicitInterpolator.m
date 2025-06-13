classdef SimpAllExplicitInterpolator < handle

   properties (Access = private)
        ndim
        matA
        matB
   end

    methods  (Access = public)
        function obj = SimpAllExplicitInterpolator(cParams)
            obj.init(cParams);
        end
        
        function [mu,kappa] = computeConsitutiveTensor(obj,rho)
            mu    = obj.computeMuFunction(rho(1)); %rho{1}
            kappa = obj.computeKappaFunction(rho(1));  %rho{1}
        end

        function [dmu,dkappa] = computeConsitutiveTensorDerivative(obj,rho)
            dmu{1}    = obj.computeMuDerivative(rho(1)); %rho{1}
            dkappa{1} = obj.computeKappaDerivative(rho(1)); %rho{1}
        end          
    end
     
    methods (Access = private)
        function init(obj,cParams)
            obj.matA = cParams.matA;
            obj.matB = cParams.matB;
            obj.computeNDim(cParams);
        end        

        function computeNDim(obj,cParams)
            switch cParams.dim
                case '2D'
                  obj.ndim = 2;
                case '3D'
                  obj.ndim = 3;
            end
        end                       
        
        function mu = computeMuFunction(obj,rho)
            m0   = obj.matA.shear;
            m1   = obj.matB.shear;
            k0   = obj.matA.bulk;
            k1   = obj.matB.bulk;
            eta0 = obj.computeEtaMu(m0,k0);
            eta1 = obj.computeEtaMu(m1,k1);
            c    = obj.computeCoeff(m0,m1,eta0,eta1);
            mu   = obj.computeRationalFunction(rho,c);
        end
        
        function kappa = computeKappaFunction(obj,rho)
            m0    = obj.matA.shear;
            m1    = obj.matB.shear;
            k0    = obj.matA.bulk;
            k1    = obj.matB.bulk;
            eta0  = obj.computeEtaKappa(m0);
            eta1  = obj.computeEtaKappa(m1);
            c     = obj.computeCoeff(k0,k1,eta0,eta1);
            kappa = obj.computeRationalFunction(rho,c);
        end

        function dmu = computeMuDerivative(obj,rho)
            m0    = obj.matA.shear;
            m1    = obj.matB.shear;
            k0    = obj.matA.bulk;
            k1    = obj.matB.bulk;
            eta0  = obj.computeEtaMu(m0,k0);
            eta1  = obj.computeEtaMu(m1,k1);
            c     = obj.computeCoeff(m0,m1,eta0,eta1);
            dmu   = obj.computeRationalDerivative(rho,c);
        end

        function dkappa = computeKappaDerivative(obj,rho)
            m0     = obj.matA.shear;
            m1     = obj.matB.shear;
            k0     = obj.matA.bulk;
            k1     = obj.matB.bulk;
            eta0   = obj.computeEtaKappa(m0);
            eta1   = obj.computeEtaKappa(m1);
            c      = obj.computeCoeff(k0,k1,eta0,eta1);
            dkappa = obj.computeRationalDerivative(rho,c);
        end

        function etaMu = computeEtaMu(obj,mu,kappa)
            N   = obj.ndim;
            num = -mu.*(4.*mu - kappa.*N^2 - 2.*mu.*N^2 + 2.*mu.*N);
            den = 2.*N.*(kappa + 2*mu);
            etaMu = num./den;
        end
        
        function etaKappa = computeEtaKappa(obj,mu)
            N   = obj.ndim;
            num = 2.*mu.*(N-1);
            den = N;
            etaKappa = num./den;
        end
    end
    
    methods (Access = private, Static)
        function c = computeCoeff(f0,f1,eta0,eta1)
            c.n01 =  -(f1 - f0).*(eta1 - eta0);
            c.n0  =  f0.*(f1 + eta0);
            c.n1  =  f1.*(f0 + eta1);
            c.d0  =  (f1 + eta0);
            c.d1  =  (f0 + eta1);
        end
        
        function f = computeRationalFunction(rho,s)
            rho = rho.fun.fValues;
            n01 = s.n01;
            n0  = s.n0;
            n1  = s.n1;
            d0  = s.d0;
            d1  = s.d1;
            num = n01.*(1-rho).*(rho) + n0.*(1-rho) + n1.*rho;
            den = d0.*(1-rho) + d1.*rho;
            f   = num./den;
        end

        function f = computeRationalDerivative(rho,s)
            n01    = s.n01;
            n0     = s.n0;
            n1     = s.n1;
            d0     = s.d0;
            d1     = s.d1;
            num    = n01.*(1-rho).*(rho) + n0.*(1-rho) + n1.*rho;
            den    = d0.*(1-rho) + d1.*rho;
            derNum = (n01.*(1-2.*rho)-n0+n1).*den - num.*(-d0+d1);
            f      = derNum ./ (den.^2);
        end
    end
end