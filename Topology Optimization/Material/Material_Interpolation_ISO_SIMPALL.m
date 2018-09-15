classdef Material_Interpolation_ISO_SIMPALL < Material_Interpolation
    properties
        mu_func
        kappa_func
        dmu_func
        dlam_func
        mu_sym
        kappa_sym        
    end
    methods
        function obj=Material_Interpolation_ISO_SIMPALL(MaterialParameters)
            obj.rho_plus=MaterialParameters.rho_plus;
            obj.rho_minus=MaterialParameters.rho_minus;
            obj.E_plus=MaterialParameters.E_plus;
            obj.E_minus=MaterialParameters.E_minus;
            obj.nu_plus=MaterialParameters.nu_plus;
            obj.nu_minus=MaterialParameters.nu_minus;
        end
        function computeSymProps(obj, rho)
            ngauss=length(rho(1,:));
            
            rho_plus=obj.rho_plus;
            rho_minus=obj.rho_minus;
            E_plus=obj.E_plus;
            E_minus=obj.E_minus;
            nu_plus=obj.nu_plus;
            nu_minus=obj.nu_minus;
            
            syms c1 c2 c3 c4
            syms gamm
            
            eq_mu = [ (c1*rho_plus^2 + c2*rho_plus + 1)/(c4 + rho_plus*c3) - E_plus/(2*nu_plus + 2);
                (c1*rho_minus^2 + c2*rho_minus + 1)/(c4 + rho_minus*c3) - E_minus/(2*nu_minus + 2);
                (c2 + 2*rho_plus*c1)/(c4 + rho_plus*c3) - (c3*(c1*rho_plus^2 + c2*rho_plus + 1))/(c4 + rho_plus*c3)^2 + (2*E_plus*(E_minus - E_plus + E_minus*nu_plus - E_plus*nu_minus))/((rho_plus - rho_minus)*(nu_plus + 1)^2*(3*E_minus + E_plus - E_minus*nu_plus + E_plus*nu_minus));
                (c2 + 2*rho_minus*c1)/(c4 + rho_minus*c3) - (c3*(c1*rho_minus^2 + c2*rho_minus + 1))/(c4 + rho_minus*c3)^2 + (2*E_minus*(E_minus - E_plus + E_minus*nu_plus - E_plus*nu_minus))/((rho_plus - rho_minus)*(nu_minus + 1)^2*(E_minus + 3*E_plus + E_minus*nu_plus - E_plus*nu_minus))];
            
            eq_kappa = [E_plus/(2*nu_plus - 2) + (c1*rho_plus^2 + c2*rho_plus + 1)/(c4 + rho_plus*c3);
                E_minus/(2*nu_minus - 2) + (c1*rho_minus^2 + c2*rho_minus + 1)/(c4 + rho_minus*c3);
                (c2 + 2*rho_plus*c1)/(c4 + rho_plus*c3) - (c3*(c1*rho_plus^2 + c2*rho_plus + 1))/(c4 + rho_plus*c3)^2 + (E_plus*(E_minus - E_plus - E_minus*nu_plus + E_plus*nu_minus))/((rho_plus - rho_minus)*(nu_plus - 1)^2*(E_minus + E_plus + E_minus*nu_plus - E_plus*nu_minus));
                (c2 + 2*rho_minus*c1)/(c4 + rho_minus*c3) - (c3*(c1*rho_minus^2 + c2*rho_minus + 1))/(c4 + rho_minus*c3)^2 + (E_minus*(E_minus - E_plus - E_minus*nu_plus + E_plus*nu_minus))/((rho_plus - rho_minus)*(nu_minus - 1)^2*(E_minus + E_plus - E_minus*nu_plus + E_plus*nu_minus))];
            
            coef_kappa = (struct2array(solve(eq_kappa,[c1,c2,c3,c4])));
            coef_mu = (struct2array(solve(eq_mu,[c1,c2,c3,c4])));
            
            
            obj.mu_sym=simplify((coef_mu(1)*gamm^2 + coef_mu(2)*gamm + 1)/(coef_mu(3)*gamm + coef_mu(4)));
            obj.kappa_sym=simplify((coef_kappa(1)*gamm^2 + coef_kappa(2)*gamm + 1)/(coef_kappa(3)*gamm + coef_kappa(4)));
            
        end
        
    end
end