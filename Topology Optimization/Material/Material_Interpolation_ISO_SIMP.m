classdef Material_Interpolation_ISO_SIMP < Material_Interpolation
    
    properties (Access = protected)             
        pExp
    end
    
    methods (Access = protected)
        
        function [muS,dmuS,kS,dkS] = computeSymbolicMuKappa(obj)
            [muS,kS] = obj.computeMuAndKappaSym();
            dmuS = diff(muS);
            dkS  = diff(kS);
        end
        
        function [mu_sym,kappa_sym] = computeMuAndKappaSym(obj)
            rho1 = obj.rho1;
            rho0 = obj.rho0;
            E1   = obj.E1;
            E0   = obj.E0;
            nu1  = obj.nu1;
            nu0  = obj.nu0;
            p = obj.pExp;
            
            rho = sym('rho','real');
            mu_sym = (E1*(rho - rho0)^p)/((2*nu1 + 2)*(rho1 - rho0)^p) - (E0*((rho - rho0)^p/(rho1 - rho0)^p - 1))/(2*nu0 + 2);
            kappa_sym = -(((2*E0*((rho - rho0)^p/(rho1 - rho0)^p - 1))/(nu0 + 1) - (2*E1*(rho - rho0)^p)/...
                ((rho1 - rho0)^p*(nu1 + 1)))*((E0*((rho - rho0)^p/...
                (rho1 - rho0)^p - 1))/(2*(nu0 + 1)) - (E0*nu0*((rho - rho0)^p/...
                (rho1 - rho0)^p - 1))/(nu0^2 - 1) - (E1*(rho - rho0)^p)/(2*(rho1 - rho0)^p*(nu1 + 1)) ...
                + (E1*nu1*(rho - rho0)^p)/((nu1^2 - 1)*(rho1 - rho0)^p)))/(((2*((E0*nu0*((rho - rho0)^p/...
                (rho1 - rho0)^p - 1))/(nu0^2 - 1) - (E1*nu1*(rho - rho0)^p)/((nu1^2 - 1)*(rho1 - rho0)^p)))/...
                ((E0*((rho - rho0)^p/(rho1 - rho0)^p - 1))/(nu0 + 1) - (E0*nu0*((rho - rho0)^p/...
                (rho1 - rho0)^p - 1))/(nu0^2 - 1) - (E1*(rho - rho0)^p)/((rho1 - rho0)^p*(nu1 + 1)) + ...
                (E1*nu1*(rho - rho0)^p)/((nu1^2 - 1)*(rho1 - rho0)^p)) + 2)*((E0*((rho - rho0)^p/(rho1 - rho0)^p...
                - 1))/(nu0 + 1) - (E0*nu0*((rho - rho0)^p/(rho1 - rho0)^p - 1))/(nu0^2 - 1) - (E1*(rho - rho0)^p)/...
                ((rho1 - rho0)^p*(nu1 + 1)) + (E1*nu1*(rho - rho0)^p)/((nu1^2 - 1)*(rho1 - rho0)^p)));            
        end
        
    end
end
