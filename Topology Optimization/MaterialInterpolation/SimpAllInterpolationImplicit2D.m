classdef SimpAllInterpolationImplicit2D < SimpAllInterpolationImplicit
    
   methods  (Access = public)
        
        function obj = SimpAllInterpolationImplicit2D(cParams)
            obj.init(cParams);
            obj.computeNstre();
            obj.computeSymbolicInterpolationFunctions();
            obj.dmu0 = obj.computeDmu0();
            obj.dmu1 = obj.computeDmu1();
            obj.dk0  = obj.computeDKappa0();
            obj.dk1  = obj.computeDKappa1();
        end

   end
        
    methods  (Access = protected)
        
        function [pMu,pKappa] = computePolarizationTensorAsMuKappa(eMatrix,eInclusion,nuMatrix,nuInclusion)
            coef = obj.compute2Dcoefficients(eMatrix,eInclusion,nuMatrix,nuInclusion);
            [p1,p2] = obj.computeP1P2(coef);
            [pMu,pKappa] = obj.computePkappaPmu(p1, p2);
        end
        
        function dmu0 = computeDmu0(obj)
            E1  = obj.matProp.E1;
            E0  = obj.matProp.E0;
            nu1 = obj.matProp.nu1;
            nu0 = obj.matProp.nu0;
            mu0 = obj.matProp.mu0;
            [pMu0,~] = obj.computePolarizationTensorAsMuKappa(E0,E1,nu0,nu1);
            dmu0     = -mu0*pMu0;
        end
        
        function dmu1 = computeDmu1(obj)
            E1  = obj.matProp.E1;
            E0  = obj.matProp.E0;
            nu1 = obj.matProp.nu1;
            nu0 = obj.matProp.nu0;
            mu1 = obj.matProp.mu1;
            [pMu1,~] = obj.computePolarizationTensorAsMuKappa(E1,E0,nu1,nu0);
            dmu1     = mu1*pMu1;
        end
        
        function dkappa0 = computeDKappa0(obj)
            E1     = obj.matProp.E1;
            E0     = obj.matProp.E0;
            nu1    = obj.matProp.nu1;
            nu0    = obj.matProp.nu0;
            kappa0 = obj.matProp.kappa0;
            [~,pKappa0] = obj.computePolarizationTensorAsMuKappa(E0,E1,nu0,nu1);
            dkappa0     = -kappa0*pKappa0;
        end
        
        function dkappa1 = computeDKappa1(obj)
            E1     = obj.matProp.E1;
            E0     = obj.matProp.E0;
            nu1    = obj.matProp.nu1;
            nu0    = obj.matProp.nu0;
            kappa1 = obj.matProp.kappa1;
            [~,pKappa1] = obj.computePolarizationTensorAsMuKappa(E1,E0,nu1,nu0);
            dkappa1     = kappa1*pKappa1;
        end

    end
    
    methods  (Access = protected, Static)
    
        function coef = compute2Dcoefficients(eMatrix,eInclusion,nuMatrix,nuInclusion)
            coef.a    = (1 + nuMatrix)/(1 - nuMatrix);
            coef.b    = (3 - nuMatrix)/(1 + nuMatrix);
            coef.gam  = eInclusion/eMatrix;
            coef.tau1 = (1 + nuInclusion)/(1 + nuMatrix);
            coef.tau2 = (1 - nuInclusion)/(1 - nuMatrix);
            coef.tau3 = (nuInclusion*(3*nuMatrix - 4) + 1)/(nuMatrix*(3*nuMatrix - 4) + 1);
        end
        
        function [p1,p2] = computeP(s)
            a    = s.a;
            b    = s.b;
            gam  = s.gam;
            tau1 = s.tau1;
            tau2 = s.tau2;
            tau3 = s.tau3;
            p1   = 1/(b*gam+tau1)*(1+b)*(tau1-gam);
            p2   = 0.5*(a-b)/(b*gam+tau1)*(gam*(gam-2*tau3)+tau1*tau2)/(a*gam+tau2);
        end
        
        function [pMu,pKappa] = computePKappaMu(p1, p2)
            pMu    = p1;
            pKappa = 2*p2 + p1;
        end
        
    end
  
end