classdef Interpolation_ISO_SIMP_P3 < Interpolation
    properties
        mu_func
        kappa_func
        dmu_func
        dlam_func
        p
    end
    methods
        function obj=Interpolation_ISO_SIMP_P3(Bounds)
            obj.rho_plus=Bounds.rho_plus;
            obj.rho_minus=Bounds.rho_minus;
            obj.E_plus=Bounds.E_plus;
            obj.E_minus=Bounds.E_minus;
            obj.nu_plus=Bounds.nu_plus;
            obj.nu_minus=Bounds.nu_minus;
            obj.p=3;
        end
        function matProps=computeMatProp(obj, rho)
            ngauss=length(rho(1,:));
            if isempty(obj.mu_func)
                rho_plus=obj.rho_plus;
                rho_minus=obj.rho_minus;
                E_plus=obj.E_plus;
                E_minus=obj.E_minus;
                nu_plus=obj.nu_plus;
                nu_minus=obj.nu_minus;
                
                mu=sparse(length(rho(:,1)),ngauss);
                kappa=sparse(length(rho(:,1)),ngauss);
                
                gamm = sym('gamm','real');
                mu_sym = (E_plus*(gamm - rho_minus)^obj.p)/((2*nu_plus + 2)*(rho_plus - rho_minus)^obj.p) - (E_minus*((gamm - rho_minus)^obj.p/(rho_plus - rho_minus)^obj.p - 1))/(2*nu_minus + 2);
                kappa_sym = -(((2*E_minus*((gamm - rho_minus)^obj.p/(rho_plus - rho_minus)^obj.p - 1))/(nu_minus + 1) - (2*E_plus*(gamm - rho_minus)^obj.p)/...
                    ((rho_plus - rho_minus)^obj.p*(nu_plus + 1)))*((E_minus*((gamm - rho_minus)^obj.p/...
                    (rho_plus - rho_minus)^obj.p - 1))/(2*(nu_minus + 1)) - (E_minus*nu_minus*((gamm - rho_minus)^obj.p/...
                    (rho_plus - rho_minus)^obj.p - 1))/(nu_minus^2 - 1) - (E_plus*(gamm - rho_minus)^obj.p)/(2*(rho_plus - rho_minus)^obj.p*(nu_plus + 1)) ...
                    + (E_plus*nu_plus*(gamm - rho_minus)^obj.p)/((nu_plus^2 - 1)*(rho_plus - rho_minus)^obj.p)))/(((2*((E_minus*nu_minus*((gamm - rho_minus)^obj.p/...
                    (rho_plus - rho_minus)^obj.p - 1))/(nu_minus^2 - 1) - (E_plus*nu_plus*(gamm - rho_minus)^obj.p)/((nu_plus^2 - 1)*(rho_plus - rho_minus)^obj.p)))/...
                    ((E_minus*((gamm - rho_minus)^obj.p/(rho_plus - rho_minus)^obj.p - 1))/(nu_minus + 1) - (E_minus*nu_minus*((gamm - rho_minus)^obj.p/...
                    (rho_plus - rho_minus)^obj.p - 1))/(nu_minus^2 - 1) - (E_plus*(gamm - rho_minus)^obj.p)/((rho_plus - rho_minus)^obj.p*(nu_plus + 1)) + ...
                    (E_plus*nu_plus*(gamm - rho_minus)^obj.p)/((nu_plus^2 - 1)*(rho_plus - rho_minus)^obj.p)) + 2)*((E_minus*((gamm - rho_minus)^obj.p/(rho_plus - rho_minus)^obj.p...
                    - 1))/(nu_minus + 1) - (E_minus*nu_minus*((gamm - rho_minus)^obj.p/(rho_plus - rho_minus)^obj.p - 1))/(nu_minus^2 - 1) - (E_plus*(gamm - rho_minus)^obj.p)/...
                    ((rho_plus - rho_minus)^obj.p*(nu_plus + 1)) + (E_plus*nu_plus*(gamm - rho_minus)^obj.p)/((nu_plus^2 - 1)*(rho_plus - rho_minus)^obj.p)));
                
                dmu = diff(mu_sym);
                d = 2; % the dimension, 2 in 2D
                lam = kappa_sym - 2/d*mu_sym;
                dlam = diff(lam);
                
                obj.mu_func = matlabFunction(mu_sym);
                obj.kappa_func = matlabFunction(kappa_sym);
                obj.dmu_func=matlabFunction(dmu);
                obj.dlam_func=matlabFunction(dlam);
            end
            for igauss=1:ngauss
                mu(:,igauss)= obj.mu_func(rho(:,igauss));
                kappa(:,igauss)= obj.kappa_func(rho(:,igauss));
                dC(1,1,:,igauss)= 2*obj.dmu_func(rho(:,igauss))+obj.dlam_func(rho(:,igauss));
                dC(1,2,:,igauss)= obj.dlam_func(rho(:,igauss));
                dC(1,3,:,igauss)= zeros(length(rho(:,igauss)),1);
                dC(2,1,:,igauss)= obj.dlam_func(rho(:,igauss));
                dC(2,2,:,igauss)= 2*obj.dmu_func(rho(:,igauss))+obj.dlam_func(rho(:,igauss));
                dC(2,3,:,igauss)= zeros(length(rho(:,igauss)),1);
                dC(3,1,:,igauss)= zeros(length(rho(:,igauss)),1);
                dC(3,2,:,igauss)= zeros(length(rho(:,igauss)),1) ;
                dC(3,3,:,igauss)= obj.dmu_func(rho(:,igauss));
            end
            
            matProps=struct;
            matProps.mu=mu;
            matProps.kappa=kappa;
            matProps.dC=dC;
        end
        
    end
end