classdef SimpAllInterpolationExplicit < MaterialInterpolation
     
    methods (Access = protected)
        
        function [mS,dmS] = computeMuSymbolicFunctionAndDerivative(obj)
            mS = obj.computeSymMu();
            dmS = diff(mS);
        end
        
        function [kS,dkS] = computeKappaSymbolicFunctionAndDerivative(obj)
            kS = obj.computeSymKappa();
            dkS = diff(kS);
        end        
        
        function mu = computeSymMu(obj)
            m0 = obj.matProp.mu0;
            m1 = obj.matProp.mu1;
            k0 = obj.matProp.kappa0;
            k1 = obj.matProp.kappa1;            
            eta0 = obj.computeEtaMu(m0,k0);
            eta1 = obj.computeEtaMu(m1,k1);
            c  = obj.computeCoeff(m0,m1,eta0,eta1);
            mu = obj.computeRationalFunction(c);            
        end
        
        function kappa = computeSymKappa(obj)
            m0 = obj.matProp.mu0;
            m1 = obj.matProp.mu1;
            k0 = obj.matProp.kappa0;
            k1 = obj.matProp.kappa1;            
            eta0 = obj.computeEtaKappa(m0);
            eta1 = obj.computeEtaKappa(m1);
            c     = obj.computeCoeff(k0,k1,eta0,eta1);
            kappa = obj.computeRationalFunction(c);            
        end
        
    end
    
    methods (Access = protected, Static)
        
        function eta = computeEtaMu(mu,kappa)
            eta = (kappa*mu)/(2*mu+kappa);
        end
        
        function eta = computeEtaKappa(mu)
            eta = mu;
        end
        
        function c = computeCoeff(f0,f1,eta0,eta1)
            c.n01 =  -(f1 - f0)*(eta1 - eta0);
            c.n0 =  f0*(f1 + eta0);
            c.n1 =  f1*(f0 + eta1);
            c.d0 =  (f1 + eta0);
            c.d1 =  (f0 + eta1);
        end
        
        function f = computeRationalFunction(s)
            rho = sym('rho','positive');            
            n01 = s.n01;
            n0  = s.n0;
            n1  = s.n1;
            d0  = s.d0;
            d1  = s.d1;
            num = n01*(1-rho)*(rho) + n0*(1-rho) + n1*rho;
            den = d0*(1-rho)+d1*rho;
            f = num/den;
        end
        
    end    
    
end