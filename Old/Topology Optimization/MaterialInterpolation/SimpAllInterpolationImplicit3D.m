classdef SimpAllInterpolationImplicit3D < SimpAllInterpolationImplicit
    
   methods  (Access = public)
        
        function obj = SimpAllInterpolationImplicit3D(cParams)
            obj.init(cParams);
            obj.computeNstre();
            obj.dmu0 = obj.computeDmu0();
            obj.dmu1 = obj.computeDmu1();
            obj.dk0  = obj.computeDKappa0();
            obj.dk1  = obj.computeDKappa1();
            obj.computeSymbolicInterpolationFunctions();
        end

   end
        
    methods  (Access = protected)
        
        function [dMu,dKappa] = computePolarizationParametersAsMuKappa(obj, eMatrix,eInclusion,nuMatrix,nuInclusion)
            [m1,m2] = obj.compute3Dcoefficients(eMatrix,eInclusion,nuMatrix,nuInclusion);
            [dMu,dKappa] = obj.computeDkappaDmu(m1,m2);
        end
        
        function dmu0 = computeDmu0(obj)
            E1  = obj.matProp.E1;
            E0  = obj.matProp.E0;
            nu1 = obj.matProp.nu1;
            nu0 = obj.matProp.nu0;
            mu0 = obj.matProp.mu0;
            mu1 = obj.matProp.mu1;
            [dMu0,~] = obj.computePolarizationParametersAsMuKappa(E0,E1,nu0,nu1);
            qmu0 = dMu0/(mu0*(mu1-mu0));
            dmu0 = mu0*(mu1-mu0)*qmu0;
        end
        
        function dmu1 = computeDmu1(obj)
            E1  = obj.matProp.E1;
            E0  = obj.matProp.E0;
            nu1 = obj.matProp.nu1;
            nu0 = obj.matProp.nu0;
            mu0 = obj.matProp.mu0;
            mu1 = obj.matProp.mu1;
            [dMu1,~] = obj.computePolarizationParametersAsMuKappa(E1,E0,nu1,nu0);
            qmu1 = dMu1/(mu1*(mu0-mu1));
            dmu1 = mu1*(mu1-mu0)*qmu1;
        end
        
        function dkappa0 = computeDKappa0(obj)
            E1     = obj.matProp.E1;
            E0     = obj.matProp.E0;
            nu1    = obj.matProp.nu1;
            nu0    = obj.matProp.nu0;
            kappa0 = obj.matProp.kappa0;
            kappa1 = obj.matProp.kappa1;
            [~,dKappa0] = obj.computePolarizationParametersAsMuKappa(E0,E1,nu0,nu1);
            qKappa0 = dKappa0/(kappa0*(kappa1-kappa0));
            dkappa0 = kappa0*(kappa1-kappa0)*qKappa0;
        end
        
        function dkappa1 = computeDKappa1(obj)
            E1     = obj.matProp.E1;
            E0     = obj.matProp.E0;
            nu1    = obj.matProp.nu1;
            nu0    = obj.matProp.nu0;
            kappa0 = obj.matProp.kappa0;
            kappa1 = obj.matProp.kappa1;
            [~,dKappa1] = obj.computePolarizationParametersAsMuKappa(E1,E0,nu1,nu0);
            qKappa1 = dKappa1/(kappa1*(kappa0-kappa1));
            dkappa1 = kappa1*(kappa1-kappa0)*qKappa1;
        end

    end
    
    methods  (Access = protected, Static)
    
        function [m1,m2] = compute3Dcoefficients(eMatrix,eInclusion,nuMatrix,nuInclusion)
            mu      = eMatrix/(2*(1+nuMatrix));
            muI     = eInclusion/(2*(1+nuInclusion));
            mu2     = mu - muI;
            lambda  = eMatrix*nuMatrix/((1+nuMatrix)*(1-2*nuMatrix));
            lambdaI = eInclusion*nuInclusion/((1+nuInclusion)*(1-2*nuInclusion));
            lam2    = lambda - lambdaI;
            m1n     = 15*mu*mu2*(nuMatrix-1);
            m1d     = 15*mu*(1-nuMatrix)+2*mu2*(5*nuMatrix-4);
            m1      = m1n/m1d;
            m2n     = lam2*(15*mu*lambda*(1-nuMatrix) + 2*lambda*mu2*(5*nuMatrix-4)) - 2*mu2*(lambda*mu2-5*mu*nuMatrix*lam2);
            m2d     = 5*mu2*(3*mu*lambda*(1-nuMatrix)-3*mu*nuMatrix*lam2-lambda*mu2*(1-2*nuMatrix));
            m2      = m2n/m2d;
        end
        
        function [dMu,dKappa] = computeDkappaDmu(m1,m2)
            dMu    = m1;
            dKappa = m1*m2+2/3*m1;
        end
       
    end
  
end