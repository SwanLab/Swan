classdef Interpolation_ISO_SIMPALL < Interpolation
    properties
    end
    methods
        function obj=Interpolation_ISO_SIMPALL(HSbounds)
            obj.rho_plus=HSbounds.rho_plus;
            obj.rho_minus=HSbounds.rho_minus;
            obj.E_plus=HSbounds.E_plus;
            obj.E_minus=HSbounds.E_minus;
            obj.nu_plus=HSbounds.nu_plus;
            obj.nu_minus=HSbounds.nu_minus;
        end
        function matProps=computeMatProp(obj, rho)
            rho_plus=obj.rho_plus;
            rho_minus=obj.rho_minus;
            E_plus=obj.E_plus;
            E_minus=obj.E_minus;
            nu_plus=obj.nu_plus;
            nu_minus=obj.nu_minus;
            ngauss=length(rho(1,:));
            mu=sparse(length(rho(:,1)),ngauss);
            kappa=sparse(length(rho(:,1)),ngauss);
            p1=sparse(length(rho(:,1)),ngauss);
            p2=sparse(length(rho(:,1)),ngauss);
            
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
            
            mu_sym=simplify((coef_mu(1)*gamm^2 + coef_mu(2)*gamm + 1)/(coef_mu(3)*gamm + coef_mu(4)));
            kappa_sym=simplify((coef_kappa(1)*gamm^2 + coef_kappa(2)*gamm + 1)/(coef_kappa(3)*gamm + coef_kappa(4)));
            
            dmu = diff(mu_sym);
            d = 2; % the dimension, 2 in 2D
            lam = kappa_sym - 2/d*mu_sym;
            dlam = diff(lam);
            
            mu_func = matlabFunction(mu_sym);
            kappa_func = matlabFunction(kappa_sym);
            dmu_func=matlabFunction(dmu);
            dlam_func=matlabFunction(dlam);     
            
            for igauss=1:ngauss
                mu(:,igauss)=mu_func(rho(:,igauss));
                kappa(:,igauss)=kappa_func(rho(:,igauss));
                dC(1,1,:,igauss)=2*dmu_func(rho(:,igauss))+dlam_func(rho(:,igauss));
                dC(1,2,:,igauss)=  dlam_func(rho(:,igauss));
                dC(1,3,:,igauss)= zeros(length(rho(:,igauss)),1);
                dC(2,1,:,igauss)= dlam_func(rho(:,igauss));
                dC(2,2,:,igauss)= 2*dmu_func(rho(:,igauss))+dlam_func(rho(:,igauss));
                dC(2,3,:,igauss)=  zeros(length(rho(:,igauss)),1);
                dC(3,1,:,igauss)=   zeros(length(rho(:,igauss)),1);
                dC(3,2,:,igauss)=  zeros(length(rho(:,igauss)),1) ;
                dC(3,3,:,igauss)=    dmu_func(rho(:,igauss));
            end
            
            matProps=struct;
            matProps.mu=mu;
            matProps.kappa=kappa;
            matProps.dC=dC;
        end
        
    end
end