classdef Interpolation_ISO_SIMPALL < Interpolation
    properties
        mu
        kappa
    end
    methods 
        function obj=Interpolation_ISO_SIMPALL(HSbounds)
            obj.gamma_plus=HSbounds.gamma_plus;
            obj.gamma_minus=HSbounds.gamma_minus;
            obj.E_plus=HSbounds.E_plus;
            obj.E_minus=HSbounds.E_minus;
            obj.nu_plus=HSbounds.nu_plus;
            obj.nu_minus=HSbounds.nu_minus;
        end
        function matProps=computeMatProp(obj, gamma)
            gamma_plus=obj.gamma_plus;
            gamma_minus=obj.gamma_minus;
            E_plus=obj.E_plus;
            E_minus=obj.E_minus;
            nu_plus=obj.nu_plus;
            nu_minus=obj.nu_minus;
            ngauss=length(gamma(1,:));
            syms c1 c2 c3 c4
            syms gamm
            
            eq_mu = [ (c1*gamma_plus^2 + c2*gamma_plus + 1)/(c4 + gamma_plus*c3) - E_plus/(2*nu_plus + 2);
                (c1*gamma_minus^2 + c2*gamma_minus + 1)/(c4 + gamma_minus*c3) - E_minus/(2*nu_minus + 2);
                (c2 + 2*gamma_plus*c1)/(c4 + gamma_plus*c3) - (c3*(c1*gamma_plus^2 + c2*gamma_plus + 1))/(c4 + gamma_plus*c3)^2 + (2*E_plus*(E_minus - E_plus + E_minus*nu_plus - E_plus*nu_minus))/((gamma_plus - gamma_minus)*(nu_plus + 1)^2*(3*E_minus + E_plus - E_minus*nu_plus + E_plus*nu_minus));
                (c2 + 2*gamma_minus*c1)/(c4 + gamma_minus*c3) - (c3*(c1*gamma_minus^2 + c2*gamma_minus + 1))/(c4 + gamma_minus*c3)^2 + (2*E_minus*(E_minus - E_plus + E_minus*nu_plus - E_plus*nu_minus))/((gamma_plus - gamma_minus)*(nu_minus + 1)^2*(E_minus + 3*E_plus + E_minus*nu_plus - E_plus*nu_minus))];

            eq_kappa = [E_plus/(2*nu_plus - 2) + (c1*gamma_plus^2 + c2*gamma_plus + 1)/(c4 + gamma_plus*c3);
                E_minus/(2*nu_minus - 2) + (c1*gamma_minus^2 + c2*gamma_minus + 1)/(c4 + gamma_minus*c3);
                (c2 + 2*gamma_plus*c1)/(c4 + gamma_plus*c3) - (c3*(c1*gamma_plus^2 + c2*gamma_plus + 1))/(c4 + gamma_plus*c3)^2 + (E_plus*(E_minus - E_plus - E_minus*nu_plus + E_plus*nu_minus))/((gamma_plus - gamma_minus)*(nu_plus - 1)^2*(E_minus + E_plus + E_minus*nu_plus - E_plus*nu_minus));
                (c2 + 2*gamma_minus*c1)/(c4 + gamma_minus*c3) - (c3*(c1*gamma_minus^2 + c2*gamma_minus + 1))/(c4 + gamma_minus*c3)^2 + (E_minus*(E_minus - E_plus - E_minus*nu_plus + E_plus*nu_minus))/((gamma_plus - gamma_minus)*(nu_minus - 1)^2*(E_minus + E_plus - E_minus*nu_plus + E_plus*nu_minus))];
                  
            coef_kappa = (struct2array(solve(eq_kappa,[c1,c2,c3,c4])));
            coef_mu = (struct2array(solve(eq_mu,[c1,c2,c3,c4])));
            
            mu_func = matlabFunction(simplify((coef_mu(1)*gamm^2 + coef_mu(2)*gamm + 1)/(coef_mu(3)*gamm + coef_mu(4))));
            kappa_func = matlabFunction(simplify((coef_kappa(1)*gamm^2 + coef_kappa(2)*gamm + 1)/(coef_kappa(3)*gamm + coef_kappa(4))));
            
            for igauss=1:ngauss
            mu(:,igauss) =mu_func(gamma(:,igauss));
            kappa(:,igauss) =kappa_func(gamma(:,igauss));
            end
            matProps=struct;
            matProps.mu=mu;
            matProps.kappa=kappa;
        end
    end
end