classdef SimpAllInterpolationImplicit < MaterialInterpolation
    
    methods (Access = protected)
        
        function [mS,dmS] = computeMuSymbolicFunctionAndDerivative(obj)
            dmu0 = obj.computeDmu0();
            dmu1 = obj.computeDmu1();            
            s.f0  = obj.matProp.mu0;
            s.f1  = obj.matProp.mu1;
            s.df0 = dmu0;
            s.df1 = dmu1;               
            [mS,dmS] = obj.computeParameterInterpolationAndDeirvative(s);
        end
        
        function [kS,dkS] = computeKappaSymbolicFunctionAndDerivative(obj)
            dk0 = obj.computeDmu0();
            dk1 = obj.computeDmu1();                        
            s.f0  = obj.matProp.kappa0;
            s.f1  = obj.matProp.kappa1;
            s.df0 = dk0;
            s.df1 = dk1;            
            [kS,dkS] = obj.computeParameterInterpolationAndDeirvative(s);
        end
                
        function [f,df] = computeParameterInterpolationAndDeirvative(obj,s)
            c = obj.computeCoefficients(s);
            rho = sym('rho','positive');            
            fSym  = obj.rationalFunction(c,rho);
            dfSym = obj.rationalFunctionDerivative(c,rho);
            f  = simplify(fSym);
            df = simplify(dfSym);
        end
        
        function c = computeCoefficients(obj,s)
            f1  = s.f1;
            f0  = s.f0;
            df1 = s.df1;
            df0 = s.df0;
            c1 = sym('c1','real'); 
            c2 = sym('c2','real');            
            c3 = sym('c3','real');            
            c4 = sym('c4','real');            
            coef = [c1 c2 c3 c4];
            r1 = obj.matProp.rho1;
            r0 = obj.matProp.rho0;
            eq(1) = obj.rationalFunction(coef,r1) - f1;
            eq(2) = obj.rationalFunction(coef,r0) - f0;
            eq(3) = obj.rationalFunctionDerivative(coef,r1) - df1;
            eq(4) = obj.rationalFunctionDerivative(coef,r0) - df0;
            c = solve(eq,[c1,c2,c3,c4]);
            c = struct2cell(c);    
            c = [c{:}];
        end
        
        function dmu0 = computeDmu0(obj)
            E1  = obj.matProp.E1;
            E0  = obj.matProp.E0;            
            nu1 = obj.matProp.nu1;            
            nu0 = obj.matProp.nu0;                        
            mu0 = obj.matProp.mu0;
            [pMu0,~] = obj.computePmuPkappa(E0,E1,nu0,nu1);
            dmu0    = -mu0*pMu0;
        end
        
        function dmu1 = computeDmu1(obj)
            E1  = obj.matProp.E1;
            E0  = obj.matProp.E0;            
            nu1 = obj.matProp.nu1;            
            nu0 = obj.matProp.nu0;                        
            mu1 = obj.matProp.mu1;
            [pMu1,~] = obj.computePmuPkappa(E1,E0,nu1,nu0);
            dmu1    = mu1*pMu1;
        end
        
        function dkappa0 = computeDKappa0(obj)
            E1  = obj.matProp.E1;
            E0  = obj.matProp.E0;            
            nu1 = obj.matProp.nu1;            
            nu0 = obj.matProp.nu0;                        
            kappa0 = obj.matProp.kappa0;            
            [~,pKappa0] = obj.computePmuPkappa(E0,E1,nu0,nu1);
            dkappa0 = -kappa0*pKappa0;            
        end
        
        function dkappa1 = computeDKappa1(obj)
            E1  = obj.matProp.E1;
            E0  = obj.matProp.E0;            
            nu1 = obj.matProp.nu1;            
            nu0 = obj.matProp.nu0;                        
            kappa1 = obj.matProp.kappa1;            
            [~,pKappa1] = obj.computePmuPkappa(E1,E0,nu1,nu0);
            dkappa1 = kappa1*pKappa1;            
        end                
        
    end
    
    methods (Access = protected, Static)
        
        function r = rationalFunction(coef,rho)
            c1 = coef(1);
            c2 = coef(2);
            c3 = coef(3);
            c4 = coef(4);
            n = (c1*rho^2 + c2*rho + 1);
            d = (c4 + rho*c3);
            r = n/d;
        end
        
        function dr = rationalFunctionDerivative(coef,rho)
            c1 = coef(1);
            c2 = coef(2);
            c3 = coef(3);
            c4 = coef(4);
            n1  = c2 + 2*rho*c1;
            d1  = c4 + rho*c3;
            dr1 = n1/d1;
            n2 = -c3*(c1*rho^2 + c2*rho + 1);
            d2 = (c4 + rho*c3)^2;
            dr2 = n2/d2;
            dr = dr1 + dr2;
        end
        
        function [pMu,pKappa] = computePmuPkappa(E_matrix,E_inclusion,nu_matrix,nu_inclusion)       
            a = (1+nu_matrix)/(1-nu_matrix);
            b = (3-nu_matrix)/(1+nu_matrix);
            gam = E_inclusion/E_matrix;
            tau1 = (1+nu_inclusion)/(1+nu_matrix);
            tau2 = (1-nu_inclusion)/(1-nu_matrix);
            tau3 = (nu_inclusion*(3*nu_matrix-4)+1)/(nu_matrix*(3*nu_matrix-4)+1);            
            p1 = (1/(b*gam+tau1)*(1+b)*(tau1-gam));
            p2 = (1/(b*gam+tau1)*0.5*(a-b)*(gam*(gam-2*tau3)+tau1*tau2)/(a*gam+tau2));            
            pMu = p1;
            pKappa = 2*p2 + p1;
        end        
        
    end
    
end