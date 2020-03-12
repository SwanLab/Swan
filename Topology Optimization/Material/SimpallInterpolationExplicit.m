classdef SimpallInterpolationExplicit < Material_Interpolation
    
    methods (Access = protected)
        
        function [mu,kappa] = computeAlternative(obj)
            mu = @(E,nu) E/(2*(1+nu));
            kappa = @(E,nu) E/(2*(1-nu));
            %
            mu0 = mu(obj.E0,obj.nu0);
            mu1 = mu(obj.E1,obj.nu1);
            %
            kappa0 = kappa(obj.E0,obj.nu0);
            kappa1 = kappa(obj.E1,obj.nu1);
            %
            % Auxiliar material property
            eta_mu = @(mu,kappa) (kappa*mu)/(2*mu+kappa);
            eta_kappa = @(mu,kappa) mu;
            %
            eta_mu0 = eta_mu(mu0,kappa0);
            eta_mu1 = eta_mu(mu1,kappa1); 
            %
            eta_kappa0 = eta_kappa(mu0,kappa0);
            eta_kappa1 = eta_kappa(mu1,kappa1);

            % Coeficients (n= numerator, d = denominator) of the rational function
            n01 = @(f0,f1,eta0,eta1) -(f1 - f0)*(eta1 - eta0);
            n0 = @(f0,f1,eta0) f0*(f1 + eta0);
            n1 = @(f0,f1,eta1) f1*(f0 + eta1);
            d0 = @(f1,eta0) (f1 + eta0);
            d1 = @(f0,eta1) (f0 + eta1);
            %
            n01_mu    = n01(mu0,   mu1,   eta_mu0,   eta_mu1);
            n01_kappa = n01(kappa0,kappa1,eta_kappa0,eta_kappa1);
            %
            n0_mu    = n0(mu0,   mu1,   eta_mu0);
            n0_kappa = n0(kappa0,kappa1,eta_kappa0);
            %
            n1_mu    = n1(mu0,   mu1,   eta_mu1);
            n1_kappa = n1(kappa0,kappa1,eta_kappa1);
            %
            d0_mu    = d0(mu1,   eta_mu0);
            d0_kappa = d0(kappa1,eta_kappa0);
            %
            d1_mu    = d1(mu0,   eta_mu1);
            d1_kappa = d1(kappa0,eta_kappa1);
            %
            % Density function (symbolic in order for further differentiation)
            rho = sym('rho','positive');
            %
            % SIMP-ALL as a rational function
            f = @(n01,n0,n1,d0,d1,rho) ...
                (n01*(1-rho)*(rho) + n0*(1-rho) + n1*rho)/(d0*(1-rho)+d1*rho);
            %
            mu =    f(n01_mu,   n0_mu,   n1_mu,   d0_mu,   d1_mu,   rho);
            kappa = f(n01_kappa,n0_kappa,n1_kappa,d0_kappa,d1_kappa,rho);
            
        end
        
        function [mS,dmS] = computeMuSymbolicFunctionAndDerivative(obj)
            [mu,kappa] = obj.computeAlternative();            
            mS = mu;
            dmS = diff(mS);
        end
        
        function [kS,dkS] = computeKappaSymbolicFunctionAndDerivative(obj)
           [mu,kappa] = obj.computeAlternative();                          
            kS     = kappa;
            dkS    = diff(kS);
        end
    
    end
    
end