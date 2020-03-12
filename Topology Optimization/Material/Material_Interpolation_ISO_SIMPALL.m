classdef Material_Interpolation_ISO_SIMPALL < Material_Interpolation
    
    methods (Access = protected)
        
        function [mS,dmS] = computeMuSymbolicFunctionAndDerivative(obj)
            s = obj.computeMuLimitParams();
            [f,df] = obj.computeParameterInterpolationAndDeirvative(s);
            mS   = f;
            dmS  = df;
        end
        
        function [kS,dkS] = computeKappaSymbolicFunctionAndDerivative(obj)
            s = obj.computeKappaLimitParams();
            [f,df] = obj.computeParameterInterpolationAndDeirvative(s);
            kS     = f;
            dkS    = df;
        end
        
        function s = computeMuLimitParams(obj)
            [mu0,mu1] = obj.computeMuLimits();
            [dmu0,dmu1] = obj.computeDmuLimits();
            s.f0  = mu0;
            s.f1  = mu1;
            s.df0 = dmu0;
            s.df1 = dmu1;
        end
        
        function s = computeKappaLimitParams(obj)
            [k0,k1]   = obj.computeKappaLimits();
            [dk0,dk1] = obj.computeDkappaLimits();
            s.f0  = k0;
            s.f1  = k1;
            s.df0 = dk0;
            s.df1 = dk1;
        end
        
        function [dmu0,dmu1] = computeDmuLimits(obj)
            rho1 = obj.rho1;
            rho0 = obj.rho0;
            E1   = obj.E1;
            E0   = obj.E0;
            nu1  = obj.nu1;
            nu0  = obj.nu0;
            dmu0 = -(2*E0*(E0 - E1 + E0*nu1 - E1*nu0))/((rho1 - rho0)*(nu0 + 1)^2*(E0 + 3*E1 + E0*nu1 - E1*nu0));
            dmu1 = -(2*E1*(E0 - E1 + E0*nu1 - E1*nu0))/((rho1 - rho0)*(nu1 + 1)^2*(3*E0 + E1 - E0*nu1 + E1*nu0));
        end
        
        function [dk0,dk1] = computeDkappaLimits(obj)
            rho1 = obj.rho1;
            rho0 = obj.rho0;
            E1   = obj.E1;
            E0   = obj.E0;
            nu1  = obj.nu1;
            nu0  = obj.nu0;
            dk0 = -(E0*(E0 - E1 - E0*nu1 + E1*nu0))/((rho1 - rho0)*(nu0 - 1)^2*(E0 + E1 - E0*nu1 + E1*nu0));
            dk1 = -(E1*(E0 - E1 - E0*nu1 + E1*nu0))/((rho1 - rho0)*(nu1 - 1)^2*(E0 + E1 + E0*nu1 - E1*nu0));
        end
        
        function [f,df] = computeParameterInterpolationAndDeirvative(obj,s)
            f1  = s.f1;
            f0  = s.f0;
            df1 = s.df1;
            df0 = s.df0;
            syms c1 c2 c3 c4
            syms rho
            coef = [c1 c2 c3 c4];
            r1 = obj.rho1;
            r0 = obj.rho0;
            eq(1) = obj.rationalFunction(coef,r1)  - f1;
            eq(2) = obj.rationalFunction(coef,r0) - f0;
            eq(3) = obj.rationalFunctionDerivative(coef,r1) - df1;
            eq(4) = obj.rationalFunctionDerivative(coef,r0) - df0;
            c     = struct2array(solve(eq,[c1,c2,c3,c4]));
            fSym  = obj.rationalFunction(c,rho);
            dfSym = obj.rationalFunctionDerivative(c,rho);
            f  = simplify(fSym);
            df = simplify(dfSym);
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
        
    end
    
end